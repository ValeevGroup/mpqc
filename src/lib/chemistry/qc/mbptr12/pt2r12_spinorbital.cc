//
// pt2r12_spinorbital.cc
//
// Copyright (C) 2009 Edward Valeev
//
// Author: Edward Valeev <evaleev@vt.edu>
// Maintainer: EV
//
// This file is part of the SC Toolkit.
//
// The SC Toolkit is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The SC Toolkit is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the SC Toolkit; see the file COPYING.LIB.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//
// The U.S. Government is granted a limited license as per AL 91-7.
//

#include <fstream>
#include <numeric>
#include <cassert>
#include <chemistry/qc/mbptr12/pt2r12.h>
#include <util/misc/print.h>
#include <chemistry/qc/wfn/orbitalspace_utils.h>
#include <math/scmat/local.h>
#include <math/scmat/svd.h>
#include <chemistry/qc/mbptr12/compute_tbint_tensor.h>
#include <chemistry/qc/mbptr12/creator.h>
#include <chemistry/qc/mbptr12/container.h>

using namespace std;
using namespace sc;

#include <chemistry/qc/mbptr12/pt2r12_utils.h>

////////////////////

static ClassDesc SpinOrbitalPT2R12_cd(typeid(SpinOrbitalPT2R12),
                                      "SpinOrbitalPT2R12",
                                      1,"public Wavefunction",
                                      0,
                                      create<SpinOrbitalPT2R12>,
                                      create<SpinOrbitalPT2R12>);

SpinOrbitalPT2R12::SpinOrbitalPT2R12(const Ref<KeyVal> &keyval) : Wavefunction(keyval), B_(), X_(), V_()
{
  ///* comment out the following to test cabs singles correction
  string nfzc_str = keyval->stringvalue("nfzc",KeyValValuestring("0"));
  if (nfzc_str == "auto")  nfzc_ = molecule()->n_core_electrons()/2;
  else if (nfzc_str == "no" || nfzc_str == "false") nfzc_ = 0;
  else nfzc_ = atoi(nfzc_str.c_str());

  pt2_correction_ = keyval->booleanvalue("pt2_correction", KeyValValueboolean(true));
  omit_uocc_ = keyval->booleanvalue("omit_uocc", KeyValValueboolean(false));
  cabs_singles_ = keyval->booleanvalue("cabs_singles", KeyValValueboolean(false));
  cabs_singles_h0_ = keyval->stringvalue("cabs_singles_h0", KeyValValuestring(string("dyall")));
  cabs_singles_coupling_ = keyval->booleanvalue("cabs_singles_coupling", KeyValValueboolean(true));
  rotate_core_ = keyval->booleanvalue("rotate_core", KeyValValueboolean(true));

  rdm2_ = require_dynamic_cast<RDM<Two>*>(
        keyval->describedclassvalue("rdm2").pointer(),
        "SpinOrbitalPT2R12::ctor\n"
        );
  if (rdm2_.null()) throw InputError("missing rdm2", __FILE__, __LINE__,
                                    "rdm2", 0,
                                    this->class_desc());
  rdm1_ = rdm2_->rdm_m_1();

  // if world not given, make this the center of a new World
  Ref<WavefunctionWorld> world; world << keyval->describedclassvalue("world", KeyValValueRefDescribedClass(0));
  if (world.null())
    world = new WavefunctionWorld(keyval);
  if (world.null())
    throw InputError("SpinOrbitalPT2R12 requires a WavefunctionWorld; input did not specify it, neither could it be constructed",
                     __FILE__, __LINE__, "world");
  if (world->wfn() == 0) world->set_wfn(this);

  const bool spin_restricted = true;  // always use spin-restricted spaces
  // if omit_uocc is true, need to make an empty virtual space
  Ref<OrbitalSpace> virspace = 0;
  if (omit_uocc_) {
    virspace = new EmptyOrbitalSpace("", "", basis(), integral(), OrbitalSpace::symmetry);
  }

  Ref<RefWavefunction> ref;
  {
    Ref<Wavefunction> reference;
    reference << keyval->describedclassvalue("reference");
    if (reference.nonnull()) {
      MPQC_ASSERT(reference == rdm2_->wfn());
      ref = RefWavefunctionFactory::make(world, reference, spin_restricted,
                                         nfzc_, 0, virspace);
    }
    else  {
      ref << keyval->describedclassvalue("refwfn").pointer();
      if (ref.null() && reference.null()) throw InputError("missing reference", __FILE__, __LINE__,
                                                           "reference", 0,this->class_desc());
    }
  }

  // some defaults need to be overridden for R12WavefunctionWorld
  // spinadapted should by default be false
  {
    Ref<KeyVal> r12world_keyval = keyval;
    if (keyval->exists("spinadapted") == false) {
      Ref<AssignedKeyVal> akeyval = new AssignedKeyVal;
      akeyval->assignboolean("spinadapted", false);
      r12world_keyval = new AggregateKeyVal(keyval, akeyval);
    }
    r12world_ = new R12WavefunctionWorld(r12world_keyval, ref);
  }
  r12eval_ = new R12IntEval(r12world_);

  debug_ = keyval->intvalue("debug", KeyValValueint(0));
  r12eval_->debug(debug_);
  // this may update the accuracy of reference_ object
  this->set_desired_value_accuracy(desired_value_accuracy());
}

SpinOrbitalPT2R12::SpinOrbitalPT2R12(StateIn &s) : Wavefunction(s) {
  rdm2_ << SavableState::restore_state(s);
  rdm1_ << SavableState::restore_state(s);
  r12world_ << SavableState::restore_state(s);
  r12eval_ << SavableState::restore_state(s);
  s.get(nfzc_);
  s.get(omit_uocc_);
  s.get(cabs_singles_);
  s.get(cabs_singles_coupling_);
  s.get(debug_);
}

SpinOrbitalPT2R12::~SpinOrbitalPT2R12() {
}

void SpinOrbitalPT2R12::save_data_state(StateOut &s) {
  Wavefunction::save_data_state(s);
  SavableState::save_state(rdm2_, s);
  SavableState::save_state(rdm1_, s);
  SavableState::save_state(r12world_, s);
  SavableState::save_state(r12eval_, s);
  s.put(nfzc_);
  s.put(omit_uocc_);
  s.put(cabs_singles_coupling_);
  s.put(debug_);
}

void
SpinOrbitalPT2R12::obsolete() {
  r12eval_->obsolete();
  rdm1_->obsolete();
  rdm2_->obsolete();
  r12world_->world()->obsolete();
  r12world_->obsolete();
  Wavefunction::obsolete();
}

void
SpinOrbitalPT2R12::set_desired_value_accuracy(double acc)
{
  Function::set_desired_value_accuracy(acc);
  Ref<RefWavefunction> refwfn = r12world()->refwfn();
  if (refwfn->desired_value_accuracy_set_to_default()) {
    // reference should be computed to higher accuracy
    const double ref_acc = acc * ref_to_pt2r12_acc();
    refwfn->set_desired_value_accuracy(ref_acc);
  }
}

RefSymmSCMatrix SpinOrbitalPT2R12::hcore_mo() {
  Ref<SCMatrixKit> localkit = new LocalSCMatrixKit;
  Ref<OrbitalSpace> space = r12eval_->orbs(Alpha);
  RefSCMatrix coeffs = space->coefs();
  RefSCDimension nmodim = coeffs.rowdim();
  RefSCDimension naodim = coeffs.coldim();
  int nmo = nmodim.n();
  int nao = naodim.n();
  RefSCMatrix coeffs_nb = localkit->matrix(naodim,nmodim);
  for(int i=0; i<nao; i++) {
    for(int j=0; j<nmo; j++) {
      coeffs_nb.set_element(i,j,coeffs.get_element(i,j));
    }
  }

  Ref<OneBodyInt> Hcore=integral()->hcore();
  const double *buffer_Hcore = Hcore->buffer();
  Ref<GaussianBasisSet> basis = this->basis();
  int nshell = basis->nshell();

  RefSymmSCMatrix Hcore_ao = localkit->symmmatrix(naodim);

  for(int P=0; P<nshell; P++) {
    int nump = basis->shell(P).nfunction();
    for(int Q=0; Q<nshell; Q++) {
      int numq = basis->shell(Q).nfunction();

      Hcore->compute_shell(P,Q);
      int index = 0;
      for(int p=0; p<nump; p++) {
        int op = basis->shell_to_function(P) + p;
        for(int q=0; q<numq; q++) {
          int oq = basis->shell_to_function(Q) + q;
          index = p * numq + q;

          Hcore_ao.set_element(op,oq,buffer_Hcore[index]);
        }
      }
    }
  }

  RefSymmSCMatrix Hcore_mo = localkit->symmmatrix(nmodim);
  Hcore_mo.assign(0.0);
  Hcore_mo.accumulate_transform(coeffs_nb,Hcore_ao,SCMatrix::TransposeTransform);
  //Hcore_mo.accumulate_transform(coeffs_nb,Hcore_ao);

  if(debug_>=DefaultPrintThresholds::mostN2) {
    Hcore_mo->print("total hcore_mo");
  }

  return(Hcore_mo);
}



RefSymmSCMatrix SpinOrbitalPT2R12::hcore_mo(SpinCase1 spin) {
  const Ref<SCMatrixKit> localkit = new LocalSCMatrixKit;
  const Ref<OrbitalSpace> space = r12eval_->orbs(spin);
  const int nmo = space->rank();
  const RefSCMatrix coeffs = space->coefs();
  const RefSCDimension nmodim = new SCDimension(nmo);
  const RefSCDimension naodim = coeffs.coldim();
  const int nao = naodim.n();
  RefSCMatrix coeffs_nb = localkit->matrix(naodim,nmodim);
  for(int i=0; i<nao; i++) {
    for(int j=0; j<nmo; j++) {
      coeffs_nb.set_element(i,j,coeffs.get_element(i,j));
    }
  }

  const Ref<OneBodyInt> Hcore=integral()->hcore();
  const double *buffer_Hcore = Hcore->buffer();
  const Ref<GaussianBasisSet> basis = this->basis();
  const int nshell = basis->nshell();

  RefSymmSCMatrix Hcore_ao = localkit->symmmatrix(naodim);

  for(int P=0; P<nshell; P++) {
    int nump = basis->shell(P).nfunction();
    for(int Q=0; Q<nshell; Q++) {
      int numq = basis->shell(Q).nfunction();

      Hcore->compute_shell(P,Q);
      int index = 0;
      for(int p=0; p<nump; p++) {
        int op = basis->shell_to_function(P) + p;
        for(int q=0; q<numq; q++) {
          int oq = basis->shell_to_function(Q) + q;
          index = p * numq + q;

          Hcore_ao.set_element(op,oq,buffer_Hcore[index]);
        }
      }
    }
  }

  RefSymmSCMatrix Hcore_mo = localkit->symmmatrix(nmodim);
  Hcore_mo.assign(0.0);
  Hcore_mo.accumulate_transform(coeffs_nb,Hcore_ao,SCMatrix::TransposeTransform);

  if(debug_>=DefaultPrintThresholds::mostN2) {
    Hcore_mo.print(prepend_spincase(spin,"hcore_mo").c_str());
  }

  return(Hcore_mo);
}

#if 0
RefSymmSCMatrix SpinOrbitalPT2R12::moints() {
  Ref<SCMatrixKit> localkit = new LocalSCMatrixKit;
  Ref<OrbitalSpace> space = r12eval_->orbs(Alpha);

  RefSCMatrix coeffs = space->coefs();

  RefSCDimension nmodim = coeffs.rowdim();
  RefSCDimension naodim = coeffs.coldim();
  int nmo = nmodim.n();
  int nao = naodim.n();
  RefSCDimension nmodim_nb = new SCDimension(nmo);
  RefSCDimension naodim_nb = new SCDimension(nao);

  RefSCMatrix coeffs_nb = localkit->matrix(naodim,nmodim);
  for(int i=0; i<nao; i++) {
    for(int j=0; j<nmo; j++) {
      coeffs_nb.set_element(i,j,coeffs.get_element(i,j));
    }
  }

  RefSCDimension triangdim = new SCDimension(triang_half_INDEX(nmo-1,nmo-1)+1);

  /// getting AO integrals
  RefSymmSCMatrix aoints = localkit->symmmatrix(triangdim);

  Ref<TwoBodyInt> twoint = integral()->electron_repulsion();
  const double *buffer = twoint->buffer();

  Ref<GaussianBasisSet> basis = this->basis();

  int nshell = basis->nshell();
  for(int P = 0; P < nshell; P++) {
    int nump = basis->shell(P).nfunction();

    for(int Q = 0; Q < nshell; Q++) {
      int numq = basis->shell(Q).nfunction();

      for(int R = 0; R < nshell; R++) {
        int numr = basis->shell(R).nfunction();

        for(int S = 0; S < nshell; S++) {
          int nums = basis->shell(S).nfunction();

          twoint->compute_shell(P,Q,R,S);

          int index = 0;
          for(int p=0; p < nump; p++) {
            int op = basis->shell_to_function(P)+p;

            for(int q = 0; q < numq; q++) {
              int oq = basis->shell_to_function(Q)+q;

              for(int r = 0; r < numr; r++) {
                int oor = basis->shell_to_function(R)+r;

                for(int s = 0; s < nums; s++,index++) {
                  int os = basis->shell_to_function(S)+s;

                  aoints.set_element(triang_half_INDEX(op,oq),triang_half_INDEX(oor,os),buffer[index]);

                }
              }
            }
          }

        }
      }
    }
  }

  twoint = 0;

  /// tranformation of rs indices to MOs.
  RefSCMatrix moints_pqRS = localkit->matrix(triangdim,triangdim); // moints in chemist's notation order
  RefSymmSCMatrix mat_rs = localkit->symmmatrix(naodim);
  RefSymmSCMatrix mat_RS = localkit->symmmatrix(nmodim);
  for(int p=0; p<nmo; p++) {
    for(int q=0; q<=p; q++) {
      int ind_pq = triang_half_INDEX(p,q);
      RefSCVector vec_rs = aoints.get_row(ind_pq);
      vector_to_symmmatrix(mat_rs,vec_rs);
      mat_RS.assign(0.0);
      mat_RS.accumulate_transform(coeffs_nb,mat_rs,SCMatrix::TransposeTransform);
      symmmatrix_to_vector(vec_rs,mat_RS);
      moints_pqRS.assign_row(vec_rs,ind_pq);
    }
  }

  aoints = RefSymmSCMatrix();

  /// transformation of pq indices to MOs
  RefSymmSCMatrix moints_PQRS = localkit->symmmatrix(triangdim);  // moints in chemist's notation order
  RefSymmSCMatrix mat_pq = localkit->symmmatrix(naodim);
  RefSymmSCMatrix mat_PQ = localkit->symmmatrix(nmodim);
  for(int R=0; R<nmo; R++) {
    for(int S=0; S<=R; S++) {
      int ind_RS = triang_half_INDEX(R,S);
      RefSCVector vec_pq = moints_pqRS.get_column(ind_RS);
      vector_to_symmmatrix(mat_pq,vec_pq);
      mat_PQ.assign(0.0);
      mat_PQ.accumulate_transform(coeffs_nb,mat_pq,SCMatrix::TransposeTransform);
      symmmatrix_to_vector(vec_pq,mat_PQ);
      moints_PQRS.assign_row(vec_pq,ind_RS);
    }
  }

  moints_pqRS = RefSCMatrix();

  return(moints_PQRS);
}
#endif

RefSCMatrix SpinOrbitalPT2R12::moints(SpinCase2 pairspin) {
  Ref<SCMatrixKit> localkit = new LocalSCMatrixKit;
  Ref<OrbitalSpace> space1;
  Ref<OrbitalSpace> space2;

  if(pairspin == AlphaBeta) {
    space1 = r12eval_->orbs(Alpha);
    space2 = r12eval_->orbs(Beta);
  }
  if(pairspin == AlphaAlpha) {
    space1 = r12eval_->orbs(Alpha);
    space2 = r12eval_->orbs(Alpha);
  }
  if(pairspin == BetaBeta) {
    space1 = r12eval_->orbs(Beta);
    space2 = r12eval_->orbs(Beta);
  }

  RefSCMatrix coeffs1 = space1->coefs();
  RefSCMatrix coeffs2 = space2->coefs();

  RefSCDimension nmodim = coeffs1.rowdim();
  RefSCDimension naodim = coeffs1.coldim();
  int nmo = nmodim.n();
  int nao = naodim.n();
  RefSCDimension nmodim_nb = new SCDimension(nmo);
  RefSCDimension naodim_nb = new SCDimension(nao);

  RefSCMatrix coeffs1_nb = localkit->matrix(naodim,nmodim);
  for(int i=0; i<nao; i++) {
    for(int j=0; j<nmo; j++) {
      coeffs1_nb.set_element(i,j,coeffs1.get_element(i,j));
    }
  }
  RefSCMatrix coeffs2_nb = localkit->matrix(naodim,nmodim);
  for(int i=0; i<nao; i++) {
    for(int j=0; j<nmo; j++) {
      coeffs2_nb.set_element(i,j,coeffs2.get_element(i,j));
    }
  }

  RefSCDimension triangdim = new SCDimension(triang_half_INDEX(nmo-1,nmo-1)+1);

  /// getting AO integrals
  RefSymmSCMatrix aoints = localkit->symmmatrix(triangdim);

  Ref<TwoBodyInt> twoint = integral()->electron_repulsion();
  const double *buffer = twoint->buffer();

  Ref<GaussianBasisSet> basis = this->basis();

  int nshell = basis->nshell();
  for(int P = 0; P < nshell; P++) {
    int nump = basis->shell(P).nfunction();

    for(int Q = 0; Q < nshell; Q++) {
      int numq = basis->shell(Q).nfunction();

      for(int R = 0; R < nshell; R++) {
        int numr = basis->shell(R).nfunction();

        for(int S = 0; S < nshell; S++) {
          int nums = basis->shell(S).nfunction();

          twoint->compute_shell(P,Q,R,S);

          int index = 0;
          for(int p=0; p < nump; p++) {
            int op = basis->shell_to_function(P)+p;

            for(int q = 0; q < numq; q++) {
              int oq = basis->shell_to_function(Q)+q;

              for(int r = 0; r < numr; r++) {
                int oor = basis->shell_to_function(R)+r;

                for(int s = 0; s < nums; s++,index++) {
                  int os = basis->shell_to_function(S)+s;


                  aoints.set_element(triang_half_INDEX(op,oq),triang_half_INDEX(oor,os),buffer[index]);

                }
              }
            }
          }

        }
      }
    }
  }

  twoint = 0;

  /// tranformation of rs indices to MOs.
  RefSCMatrix moints_pqRS = localkit->matrix(triangdim,triangdim); // moints in chemist's notation order
  RefSymmSCMatrix mat_rs = localkit->symmmatrix(naodim);
  RefSymmSCMatrix mat_RS = localkit->symmmatrix(nmodim);
  for(int p=0; p<nmo; p++) {
    for(int q=0; q<=p; q++) {
      int ind_pq = triang_half_INDEX(p,q);
      RefSCVector vec_rs = aoints.get_row(ind_pq);
      vector_to_symmmatrix(mat_rs,vec_rs);
      mat_RS.assign(0.0);
      mat_RS.accumulate_transform(coeffs2_nb,mat_rs,SCMatrix::TransposeTransform);
      symmmatrix_to_vector(vec_rs,mat_RS);
      moints_pqRS.assign_row(vec_rs,ind_pq);
    }
  }

  aoints = RefSymmSCMatrix();

  /// transformation of pq indices to MOs
  RefSCMatrix moints_PQRS = localkit->matrix(triangdim,triangdim);  // moints in chemist's notation order
  RefSymmSCMatrix mat_pq = localkit->symmmatrix(naodim);
  RefSymmSCMatrix mat_PQ = localkit->symmmatrix(nmodim);
  for(int R=0; R<nmo; R++) {
    for(int S=0; S<=R; S++) {
      int ind_RS = triang_half_INDEX(R,S);
      RefSCVector vec_pq = moints_pqRS.get_column(ind_RS);
      vector_to_symmmatrix(mat_pq,vec_pq);
      mat_PQ.assign(0.0);
      mat_PQ.accumulate_transform(coeffs1_nb,mat_pq,SCMatrix::TransposeTransform);
      symmmatrix_to_vector(vec_pq,mat_PQ);
      moints_PQRS.assign_row(vec_pq,ind_RS);
    }
  }

  moints_pqRS = RefSCMatrix();

  return(moints_PQRS);
}

RefSCMatrix SpinOrbitalPT2R12::g(SpinCase2 pairspin,
                      const Ref<OrbitalSpace>& space1,
                      const Ref<OrbitalSpace>& space2) {
  const Ref<SCMatrixKit> localkit = new LocalSCMatrixKit;
  Ref<SpinMOPairIter> PQ_iter = new SpinMOPairIter(space1->rank(),space2->rank(),pairspin);
  const int nmo1 = space1->rank();
  const int nmo2 = space2->rank();
  const int pairrank = PQ_iter->nij();
  const RefSCDimension pairdim = new SCDimension(pairrank);
  RefSCMatrix G = localkit->matrix(pairdim,pairdim);
  G.assign(0.0);

  // find equivalent spaces in the registry
  Ref<OrbitalSpaceRegistry> oreg = this->r12world()->world()->tfactory()->orbital_registry();
  bool registered_space1=false, registered_space2=false;
  if (!oreg->value_exists(space1)) {
    oreg->add(make_keyspace_pair(space1));
    registered_space1 = true;
  }
  if (!oreg->value_exists(space2)) {
    oreg->add(make_keyspace_pair(space2));
    registered_space2 = true;
  }
  const string key1 = oreg->key(space1);
  Ref<OrbitalSpace> s1 = oreg->value(key1);
  const string key2 = oreg->key(space2);
  Ref<OrbitalSpace> s2 = oreg->value(key2);

  const bool antisymmetrize = (pairspin==AlphaBeta) ? false : true;
  std::vector<string> tforms;
  {
    const string tform_key = ParsedTwoBodyFourCenterIntKey::key(s1->id(),s2->id(),
                                                                     s1->id(),s2->id(),
                                                                     string("ERI"),
                                                                     string(TwoBodyIntLayout::b1b2_k1k2));
    tforms.push_back(tform_key);
  }

  r12eval_->compute_tbint_tensor<ManyBodyTensors::I_to_T,false,false>(G,TwoBodyOper::eri,s1,s1,s2,s2,
      antisymmetrize,tforms);

  if (registered_space1) oreg->remove(space1->id());
  if (registered_space2) oreg->remove(space2->id());

  return(G);
}

RefSCMatrix SpinOrbitalPT2R12::g(SpinCase2 pairspin,
                      const Ref<OrbitalSpace>& bra1,
                      const Ref<OrbitalSpace>& bra2,
                      const Ref<OrbitalSpace>& ket1,
                      const Ref<OrbitalSpace>& ket2)
{
      const Ref<SCMatrixKit> localkit = new LocalSCMatrixKit;
      Ref<SpinMOPairIter> braiter12 = new SpinMOPairIter(bra1->rank(),bra2->rank(),pairspin);
      Ref<SpinMOPairIter> ketiter12 = new SpinMOPairIter(ket1->rank(),ket2->rank(),pairspin);
      RefSCMatrix G = localkit->matrix(new SCDimension(braiter12->nij()),
                                       new SCDimension(ketiter12->nij()));
      G.assign(0.0);

      const int nket1 = ket1->rank();
      const int nket2 = ket2->rank();

      const bool antisymmetrize = (pairspin != AlphaBeta);
      const bool bra1_eq_bra2 = (*bra1 == *bra2);
      const bool ket1_eq_ket2 = (*ket1 == *ket2);


      // find equivalent spaces in the registry
      Ref<OrbitalSpaceRegistry> oreg = this->r12world()->world()->tfactory()->orbital_registry();
      if (!oreg->value_exists(bra1) || !oreg->value_exists(bra2) || !oreg->value_exists(ket1) || !oreg->value_exists(ket2))
        throw ProgrammingError("SpinOrbitalPT2R12::g() -- spaces must be registered",__FILE__,__LINE__);


      // compute (bra1 ket1 | bra2 ket2)
      Ref<DistArray4> da4_b1k1_b2k2;
      {
        const string tform_b1k1_b2k2_key = ParsedTwoBodyFourCenterIntKey::key(bra1->id(),bra2->id(),
                                                                                   ket1->id(),ket2->id(),
                                                                                   string("ERI"),
                                                                                   string(TwoBodyIntLayout::b1b2_k1k2));
        Ref<TwoBodyMOIntsTransform> tform_b1k1_b2k2 = this->r12world()->world()->moints_runtime4()->get( tform_b1k1_b2k2_key );
        tform_b1k1_b2k2->compute();
        da4_b1k1_b2k2 = tform_b1k1_b2k2->ints_distarray4();
      }
      da4_b1k1_b2k2->activate();


      Ref<DistArray4> da4_b2k1_b1k2;
      const bool need_b2k1_b1k2 = !bra1_eq_bra2 && !ket1_eq_ket2 && antisymmetrize; // if(...) then essentially there is no interchangablility since the index pair are
      if (need_b2k1_b1k2)                                                           // in different spaces;
      { // compute (bra2 ket1 | bra1 ket2)
        const string tform_b2k1_b1k2_key = ParsedTwoBodyFourCenterIntKey::key(bra2->id(),bra1->id(),
                                                                                   ket1->id(),ket2->id(),
                                                                                   string("ERI"),
                                                                                   string(TwoBodyIntLayout::b1b2_k1k2));
        Ref<TwoBodyMOIntsTransform> tform_b2k1_b1k2 = this->r12world()->world()->moints_runtime4()->get( tform_b2k1_b1k2_key );
        tform_b2k1_b1k2->compute();
        da4_b2k1_b1k2 = tform_b2k1_b1k2->ints_distarray4();
        da4_b2k1_b1k2->activate();
      }

      for(braiter12->start(); *braiter12; braiter12->next())
      {
        const int b1 = braiter12->i();
        const int b2 = braiter12->j();
        const int b12 = braiter12->ij();

        const double* blk_b1b2 = da4_b1k1_b2k2->retrieve_pair_block(b1, b2, 0);

        const double* blk_b2b1 = 0;
        if (antisymmetrize && !ket1_eq_ket2 && bra1_eq_bra2)  // bra1 and bra2 are in the same space, but ket1 and ket2 not; we get b2b1 block from da4_b1k1_b2k2
          blk_b2b1 = da4_b1k1_b2k2->retrieve_pair_block(b2, b1, 0);
        else if (antisymmetrize && need_b2k1_b1k2)            // neither of the index pairs are in the same space; then get b2b1 block from da4_b2k1_b1k2
          blk_b2b1 = da4_b2k1_b1k2->retrieve_pair_block(b2, b1, 0);

        for(ketiter12->start(); *ketiter12; ketiter12->next())
        {
          const int k1 = ketiter12->i();
          const int k2 = ketiter12->j();
          const int k12 = ketiter12->ij();

          if (!antisymmetrize)
            G.set_element( b12, k12, blk_b1b2[k12] );
          else
          {
            const int k12_rect = k1 * nket2 + k2;
            if (blk_b2b1) {  // this branch is for the case when we need b2b1 block
              G.set_element( b12, k12, blk_b1b2[k12_rect] - blk_b2b1[k12_rect] );
            }
            else {          // the other case: ket1 and ket2 are in the same space; we can antisymmetrize G by subtracting k21 from k12
              const int k21_rect = k2 * nket1 + k1;
              G.set_element( b12, k12, blk_b1b2[k12_rect] - blk_b1b2[k21_rect] );
            }
          }
        }

        da4_b1k1_b2k2->release_pair_block(b1, b2, 0);
        if (blk_b2b1) {
          if (need_b2k1_b1k2)
            da4_b2k1_b1k2->release_pair_block(b2, b1, 0);
          else
            da4_b1k1_b2k2->release_pair_block(b2, b1, 0);
        }
      }
  return(G);
}

RefSCMatrix SpinOrbitalPT2R12::f(SpinCase1 spin) {
  Ref<OrbitalSpace> space = rdm1_->orbs(spin);
  Ref<OrbitalSpaceRegistry> oreg = r12world()->world()->tfactory()->orbital_registry();
  if (!oreg->value_exists(space)) {
    oreg->add(make_keyspace_pair(space));
  }
  const string key = oreg->key(space);
  space = oreg->value(key);
  RefSCMatrix fmat = r12eval_->fock(space,space,spin);

  return(fmat);
}


RefSCMatrix SpinOrbitalPT2R12::C(SpinCase2 S) {
  //return(r12energy_->C(S));
  Ref<LocalSCMatrixKit> local_matrix_kit = new LocalSCMatrixKit();
  RefSCMatrix Cmat = local_matrix_kit->matrix(r12eval_->dim_GG(S),r12eval_->dim_gg(S));
  if(S==AlphaBeta) {
    SpinMOPairIter OW_iter(r12eval_->GGspace(Alpha)->rank(), r12eval_->GGspace(Beta)->rank(), S );
    SpinMOPairIter PQ_iter(r12eval_->ggspace(Alpha)->rank(), r12eval_->ggspace(Beta)->rank(), S );
    CuspConsistentGeminalCoefficient coeff_gen(S);
    for(OW_iter.start(); int(OW_iter); OW_iter.next()) {
      for(PQ_iter.start(); int(PQ_iter); PQ_iter.next()) {
        unsigned int O = OW_iter.i();
        unsigned int W = OW_iter.j();
        unsigned int P = PQ_iter.i();
        unsigned int Q = PQ_iter.j();
        int OW = OW_iter.ij();
        int PQ = PQ_iter.ij();
        Cmat.set_element(OW,PQ,coeff_gen.C(O,W,P,Q));
      }
    }
  }
  else {
    SpinCase1 spin = (S==AlphaAlpha) ? Alpha : Beta;
    SpinMOPairIter OW_iter(r12eval_->GGspace(spin)->rank(), r12eval_->GGspace(spin)->rank(), S );
    SpinMOPairIter PQ_iter(r12eval_->ggspace(spin)->rank(), r12eval_->ggspace(spin)->rank(), S );
    CuspConsistentGeminalCoefficient coeff_gen(S);
    for(OW_iter.start(); int(OW_iter); OW_iter.next()) {
      for(PQ_iter.start(); int(PQ_iter); PQ_iter.next()) {
        unsigned int O = OW_iter.i();
        unsigned int W = OW_iter.j();
        unsigned int P = PQ_iter.i();
        unsigned int Q = PQ_iter.j();
        int OW = OW_iter.ij();
        int PQ = PQ_iter.ij();
        Cmat.set_element(OW,PQ,coeff_gen.C(O,W,P,Q));
      }
    }
  }
  return(Cmat);
}

RefSCMatrix SpinOrbitalPT2R12::transform_MO() //transformation matrix between occupied orbitals (this also defines the matrix dimension); row: new MO, column: old MO
{                                       // assume mo_density is of the same ordering as occ_sb()
  const bool debugprint = false;
  RefSymmSCMatrix mo_density =  rdm1(Alpha) + rdm1(Beta);//this will eventually read the checkpoint file. I assume they are of the dimension of occ orb space
//  mo_density.print(string("transform_MO: mo_density (occ)").c_str());
  Ref<OrbitalSpace> unscreen_occ_act = r12world()->refwfn()->occ_act_sb(Alpha);
  Ref<OrbitalSpace> occ = r12world()->refwfn()->occ_sb(Alpha);
  std::vector<int> map1 = map(*occ, *unscreen_occ_act);

  int num_occ_act = unscreen_occ_act->rank();
  int num_occ = occ->rank();
  int rdmdim = mo_density->n();
  MPQC_ASSERT(rdmdim == num_occ);
  std::vector<int> occ_act_inds;
  std::vector<int> occ_act_mask(num_occ, 0);
  for (int i = 0; i < num_occ_act; ++i)
  {
    if(map1[i] == -1)
      throw ProgrammingError("some orbital in occ_act not belong to OBS;", __FILE__,__LINE__);
    else
      {
        occ_act_inds.push_back(map1[i]);
        occ_act_mask[map1[i]] = 1;
      }
  }

  RefSCDimension mo_dim = occ->coefs()->coldim();
  RefSCMatrix TransformMat = occ->coefs()->kit()->matrix(mo_dim, mo_dim);
  TransformMat.assign(0.0);
  for (int kk = 0; kk < num_occ; ++kk)
  {
    if(occ_act_mask[kk] == 0) TransformMat.set_element(kk, kk, 1.0);
  } // don't transform unrelated orbitals
  if(debugprint)
    TransformMat.print(prepend_spincase(AlphaBeta, "transform_MO: TransformMat initial").c_str());



  Ref<LocalSCMatrixKit> local_kit = new LocalSCMatrixKit();

  //the following code is copied from refwfn.cc: PopulatedOrbitalSpace(...)
  std::vector<unsigned int> blocks = occ->block_sizes();
  int nblocks = occ->nblocks();
  int blockoffset = 0;
  for (int i = 0; i < nblocks; ++i)
  {
    std::vector<int> occ_act_orb_inds; //record the position of occ_act orbitals
    for (int j = 0; j < blocks[i]; ++j)
    {
      const int ind = blockoffset + j;
      if (occ_act_mask[ind] == 1) occ_act_orb_inds.push_back(ind);
    }
    blockoffset += blocks[i];
    const int num_occ_act = occ_act_orb_inds.size();
    SCDimension *dim = new SCDimension(num_occ_act);
    RefSCMatrix occ_act_blockmat = local_kit->matrix(dim, dim);
    for (int k1 = 0; k1 < num_occ_act; ++k1)
    {
      for (int k2 = 0; k2 < num_occ_act; ++k2)
      {
        const double element = mo_density->get_element(occ_act_orb_inds[k1], occ_act_orb_inds[k2]);
        occ_act_blockmat->set_element(k1,k2, element);
      }
    } // finish constructing a block of RDM matrix: occ_act
    RefSCMatrix UU = occ_act_blockmat->clone();
    RefSCMatrix VV = occ_act_blockmat->clone();
    RefDiagSCMatrix DD = local_kit->diagmatrix(dim); // the matrix is a postive-semidefinite matrix, do SVD
    UU.assign(0.0); VV.assign(0.0);DD.assign(0.0);
    occ_act_blockmat->svd_this(UU,DD,VV);
#if 0
    occ_act_blockmat.print(prepend_spincase(AlphaBeta, "transform_MO: occ_act_block").c_str());
    UU.print(string("transform_MO: UU").c_str());
    DD.print(string("transform_MO: DD").c_str());
    VV.print(string("transform_MO: VV").c_str());
#endif
    for (int i2 = 0; i2 < num_occ_act; ++i2)
    {
      for (int j2 = 0; j2 < num_occ_act; ++j2)
      {
        TransformMat->set_element(occ_act_orb_inds[i2], occ_act_orb_inds[j2], UU->get_element(j2, i2)); // in the final transform matrix, row:new, col:old
      }
    }// finish building MO transform matrix in the block
  }
#if 0
    TransformMat.print(string("transform_MO: final TransformMat").c_str());
#endif

   return TransformMat;
}

template<SpinOrbitalPT2R12::Tensor4_Permute HowPermute>
RefSCMatrix SpinOrbitalPT2R12::RefSCMAT4_permu(RefSCMatrix rdm2_4space_int,
                                      const Ref<OrbitalSpace> b1space,
                                      const Ref<OrbitalSpace> b2space,
                                      const Ref<OrbitalSpace> k1space,
                                      const Ref<OrbitalSpace> k2space)
{
  RefSCMatrix result;
  const unsigned int b1dim = b1space->rank();
  const unsigned int b2dim = b2space->rank();
  const unsigned int k1dim = k1space->rank();
  const unsigned int k2dim = k2space->rank();
  Ref<SpinMOPairIter> upp_pair, low_pair;
  if(HowPermute == Permute34)
  {
    result = rdm2_4space_int.clone();
    result.assign(0.0);
    upp_pair = new SpinMOPairIter(b1dim, b2dim, AlphaBeta); // the iterators are for the resultant matrices
    low_pair = new SpinMOPairIter(k2dim, k1dim, AlphaBeta);
  }
  else if(HowPermute == Permute23)
  {
    const unsigned int upp_dim = b1dim * k1dim;
    const unsigned int low_dim = b2dim * k2dim;
    RefSCDimension upp_scdim = new SCDimension(upp_dim);
    RefSCDimension low_scdim = new SCDimension(low_dim);
    result = rdm2_4space_int->kit()->matrix(upp_scdim, low_scdim);
    result.assign(0.0);
    upp_pair = new SpinMOPairIter(b1dim, k1dim, AlphaBeta);
    low_pair = new SpinMOPairIter(b2dim, k2dim, AlphaBeta);
  }
  else if(HowPermute == Permute14)
  {
      const unsigned int upp_dim = k2dim * b2dim;
      const unsigned int low_dim = k1dim * b1dim;
      RefSCDimension upp_scdim = new SCDimension(upp_dim);
      RefSCDimension low_scdim = new SCDimension(low_dim);
      result = rdm2_4space_int->kit()->matrix(upp_scdim, low_scdim);
      result.assign(0.0);
      upp_pair = new SpinMOPairIter(k2dim, b2dim, AlphaBeta);
      low_pair = new SpinMOPairIter(k1dim, b1dim, AlphaBeta);
    }


  double oneelement;
  for(upp_pair->start(); *upp_pair; upp_pair->next())   // add BetaAlpha/AlphaAlpha/BetaBeta contribution to sf_rdm
  {
    for(low_pair->start(); *low_pair; low_pair->next())
      {
         const int U1 = upp_pair->i();
         const int U2 = upp_pair->j();
         const int L1 = low_pair->i();
         const int L2 = low_pair->j();
         if(HowPermute == Permute34)
             oneelement = rdm2_4space_int.get_element(upp_pair->ij(), L2* k2dim + L1);
         else if(HowPermute == Permute23)
             oneelement = rdm2_4space_int.get_element(U1 * b2dim + L1, U2* k2dim + L2);
         else if(HowPermute == Permute14)
             oneelement = rdm2_4space_int.get_element(L2 * b2dim + U2, L1* k2dim + U1);
         else
             abort();
         result(upp_pair->ij(), low_pair->ij()) = oneelement;
      }
  }
  return result;
}








RefSymmSCMatrix SpinOrbitalPT2R12::rdm1(SpinCase1 spin)
{
  return convert_to_local_kit(rdm1_->scmat(spin));
}



RefSymmSCMatrix SpinOrbitalPT2R12::rdm2(SpinCase2 spin)
{
  // since LocalSCMatrixKit is used everywhere, convert to Local kit
  return convert_to_local_kit(rdm2_->scmat(spin));
}



void SpinOrbitalPT2R12::compute()
{
  r12world()->initialize();

  const double energy_ref = r12world()->refwfn()->energy();
  double energy_correction_r12 = 0.0;
  double energy_pt2r12[NSpinCases2];
  double B_so = 0.0;
  double V_so = 0.0;
  double X_so = 0.0; // these 3 varialbes store the final values for B, V, X
  double energy_pt2r12_sf = 0.0;
  double cabs_singles_e = 0.0;
  const bool spin_polarized = r12world()->refwfn()->spin_polarized();

  if (pt2_correction_)
  {
    // use spin-orbital version
    {
      for(int i=0; i<NSpinCases2; i++) // may comment out this part for pure cas
      {
        SpinCase2 pairspin = static_cast<SpinCase2>(i);
        double scale = 1.0;
        if (pairspin == BetaBeta && !spin_polarized) continue;
        switch (r12world()->r12tech()->ansatz()->projector())
        {
          case R12Technology::Projector_1:
            energy_pt2r12[i] = energy_PT2R12_projector1(pairspin);
            break;
          case R12Technology::Projector_2:
            energy_pt2r12[i] = energy_PT2R12_projector2(pairspin);
            break;
          default:
            abort();
        }
      }
    }

    {
      if (!spin_polarized)
      {
        energy_pt2r12[BetaBeta] = energy_pt2r12[AlphaAlpha];
        B_.push_back(B_[AlphaAlpha]);
        V_.push_back(V_[AlphaAlpha]);
        X_.push_back(X_[AlphaAlpha]);
      }
      for(int i=0; i<NSpinCases2; i++)
      {
           energy_correction_r12 +=  energy_pt2r12[i];
           B_so += B_[i];
           V_so += V_[i];
           X_so += X_[i];
      }
    }
  }


  if(cabs_singles_)
  {
    double alpha_corre = 0.0, beta_corre = 0.0, cabs_singles_corre = 0.0;
    if(cabs_singles_h0_ == string("fock"))
    {
      beta_corre =  this->cabs_singles_Fock(Beta);
      if (spin_polarized)
        alpha_corre = this->cabs_singles_Fock(Alpha);
      else
        alpha_corre = beta_corre;
      cabs_singles_e = alpha_corre + beta_corre;
      ExEnv::out0() << indent << scprintf("CABS singles (fock):                   %17.12lf",
                                        cabs_singles_e) << endl;

    }
    else if(cabs_singles_h0_ == string("dyall"))
    {
      cabs_singles_e = cabs_singles_Dyall();
      ExEnv::out0() << indent << scprintf("CABS singles (dyall):                  %17.12lf",
                                      cabs_singles_e) << endl;
    }
    else
      abort();
    ExEnv::out0() << indent << scprintf("RASSCF+CABS singles:                   %17.12lf",
                                              energy_ref + cabs_singles_e) << endl << endl;
  }

  const double energy = energy_ref + energy_correction_r12 + cabs_singles_e;

  ExEnv::out0() <<indent << scprintf("Reference energy [au]:                 %17.12lf",
                                     energy_ref) << endl;
  #if 1
    const double recomp_ref_energy = this->energy_recomputed_from_densities();
    ExEnv::out0() <<  std::endl << std::endl << indent << scprintf("Reference energy (%9s) [au]:     %17.12lf",
                                        (this->r12world()->world()->basis_df().null() ? "   recomp" : "recomp+DF"),
                                        recomp_ref_energy) << endl;
  #endif

  if(pt2_correction_)
  {
    if (r12world_->spinadapted() == false) // if true, use spinadapted algorithm, then spin component contributions are not separated.
    {
      ExEnv::out0() << indent << scprintf("Alpha-beta [2]_R12 energy [au]:        %17.12lf",
                                          energy_pt2r12[AlphaBeta]) << endl;
      ExEnv::out0() << indent << scprintf("Alpha-alpha [2]_R12 energy [au]:       %17.12lf",
                                          energy_pt2r12[AlphaAlpha]) << endl;
      if (spin_polarized)
      {
        ExEnv::out0() << indent << scprintf("Beta-beta [2]_R12 energy [au]:       %17.12lf",
                                            energy_pt2r12[BetaBeta]) << endl;
      }
      else // if not polarized, then BebaBeta contribution == AlaphAlpha contribution
      {
        ExEnv::out0() << indent << scprintf("Singlet [2]_R12 energy [au]:           %17.12lf",
                                            energy_pt2r12[AlphaBeta] - energy_pt2r12[AlphaAlpha]) << endl;
        ExEnv::out0() << indent << scprintf("Triplet [2]_R12 energy [au]:           %17.12lf",
                                            3.0*energy_pt2r12[AlphaAlpha]) << endl;
      }
    }
  }

  ExEnv::out0() << indent << scprintf("[2]_R12 energy [au]:                   %17.12lf",
                                      energy_correction_r12) << endl;
  ExEnv::out0() << indent << scprintf("Total [2]_R12 energy [au]:             %17.12lf",
                                      energy) << std::endl;
  set_energy(energy);
}

double SpinOrbitalPT2R12::compute_energy(const RefSCMatrix &hmat,
                              SpinCase2 pairspin,
                              bool print_pair_energies,
                              std::ostream& os) {
  SpinCase1 spin1 = case1(pairspin);
  SpinCase1 spin2 = case2(pairspin);
  Ref<OrbitalSpace> gg1space = r12eval_->ggspace(spin1);
  Ref<OrbitalSpace> gg2space = r12eval_->ggspace(spin2);
  SpinMOPairIter gg_iter(gg1space->rank(),gg2space->rank(),pairspin);
  double energy = 0.0;

  if (print_pair_energies) {
    os << indent << prepend_spincase(pairspin, "[2]_R12 pair energies:") << endl;
    os << indent << scprintf("    i       j        e (ij)   ") << endl;
    os << indent << scprintf("  -----   -----   ------------") << endl;
  }
  for(gg_iter.start(); int(gg_iter); gg_iter.next()) {
    int i = gg_iter.i();
    int j = gg_iter.j();
    int ij = gg_iter.ij();
    const double e_ij = hmat.get_element(ij,ij);
    if (print_pair_energies) {
      os << indent << scprintf("  %3d     %3d     %12.9lf",
                               i+1,j+1,e_ij) << endl;
    }
    energy+=hmat.get_element(ij,ij);
  }
  if (print_pair_energies)
    os << indent << endl;
#if 0
  os << indent << "energy from spin " << pairspin << " :" << energy << std::endl;
  os << indent << "trace: " << hmat.trace() << std::endl;
#endif
  return energy;
}


double SpinOrbitalPT2R12::magnetic_moment() const
{
  return r12world()->refwfn()->magnetic_moment();
}



double SpinOrbitalPT2R12::cabs_singles_Fock(SpinCase1 spin)
{
  # define printout false

  Ref<OrbitalSpace> activespace = this->r12world()->refwfn()->occ_act_sb();
  Ref<OrbitalSpace> pspace = rdm1_->orbs(spin);
  Ref<OrbitalSpace> vspace = this->r12world()->refwfn()->uocc_act_sb(spin);
  Ref<OrbitalSpace> cabsspace = this->r12world()->cabs_space(spin);

  Ref<OrbitalSpaceRegistry> oreg = this->r12world()->world()->tfactory()->orbital_registry();
  if (!oreg->value_exists(pspace))
  {
    oreg->add(make_keyspace_pair(pspace));
  }
  const string key = oreg->key(pspace);
  pspace = oreg->value(key);

  Ref<OrbitalSpace> all_virtual_space;
  if (cabs_singles_coupling_)
  {
    all_virtual_space = new OrbitalSpaceUnion("AA", "all virtuals", *vspace, *cabsspace, true);
    if (!oreg->value_exists(all_virtual_space))
      oreg->add(make_keyspace_pair(all_virtual_space));
    const string AAkey = oreg->key(all_virtual_space);
    all_virtual_space = oreg->value(AAkey);

    { // make sure that the AO space that supports all_virtual_space is known
      Ref<AOSpaceRegistry> aoreg = this->r12world()->world()->tfactory()->ao_registry();
      if (aoreg->key_exists(all_virtual_space->basis()) == false) {
        Ref<Integral> localints = this->integral()->clone();
        Ref<OrbitalSpace> mu = new AtomicOrbitalSpace("mu''", "CABS(AO)+VBS(AO)", all_virtual_space->basis(), localints);
        oreg->add(make_keyspace_pair(mu));
        aoreg->add(mu->basis(),mu);
      }
    }
  }

  Ref<OrbitalSpace> Aspace;
  if (cabs_singles_coupling_)
    Aspace = all_virtual_space;
  else
    Aspace = cabsspace;

  const unsigned int num_blocks = vspace->nblocks();// since the two spaces have the same blocksizes, we only need to define one
  const std::vector<unsigned int>& v_block_sizes = vspace->block_sizes();
  const std::vector<unsigned int>& cabs_block_sizes = cabsspace->block_sizes();
  const std::vector<unsigned int>& p_block_sizes = pspace->block_sizes(); // for debugging purposes
  const std::vector<unsigned int>& A_block_sizes = Aspace->block_sizes(); // for debugging purposes

  // fock matrices
  RefSCMatrix F_pA = r12eval_->fock(pspace,Aspace,spin);
  RefSCMatrix F_AA = r12eval_->fock(Aspace,Aspace,spin);
  RefSCMatrix F_pp = this->f(spin);
  RefSCMatrix F_pp_os = this->f(other(spin));

  // RDMs
  RefSymmSCMatrix gamma1 = rdm1_->scmat(spin);
  RefSymmSCMatrix gamma1_otherspin = rdm1_->scmat(other(spin));
  RefSymmSCMatrix gamma2_ss = this->rdm2( case12(spin,spin) ); //ss: same spin
  RefSymmSCMatrix gamma2_ba = this->rdm2( case12(spin,other(spin)) );//os: opposite spin

  // define H0 and necessary vectors
  const int no = pspace->rank();
  const int nX = Aspace->rank();
  const int nv = vspace->rank();
  const int n_cabs = cabsspace->rank();
  const int noX = no * nX;
  RefSCDimension dim_oX = new SCDimension(noX);
  RefSCDimension dim_o  = new SCDimension(no);
  RefSCDimension dim_X  = new SCDimension(nX);
  RefSymmSCMatrix H0 = gamma2_ss->kit()->symmmatrix(dim_oX); // H0 is the B matrix in our equaiton
  H0.assign(0.0);
  RefSymmSCMatrix Ixy = gamma2_ss->kit()->symmmatrix(dim_o);  // matrix to store intermediate quantities
  Ixy.assign(0.0);
  RefSCVector rhs_vector = gamma2_ss->kit()->vector(dim_oX);
  rhs_vector->assign(0.0);                                   // right hand vector


  if(!rotate_core_) // if forbit rotating core orbitals
  {
    std::vector<int> map_1_to_2 = map(*activespace, *pspace);
    for (int row_ind = 0; row_ind < no; ++row_ind)
    {
      if(map_1_to_2[row_ind] < 0)  // row_ind is a core-orbital index
      {
        for (int A1 = 0; A1 < nX; ++A1)
        {
          F_pA.set_element(row_ind, A1, 0.0);
        }
      }
    }
  }

  if (cabs_singles_coupling_) // when using the combined space, the Fock matrix component f^i_a must be zeroed since it belongs to zeroth order Hamiltonian;
   {
    unsigned int offset1 = 0;
    for (int block_counter1 = 0; block_counter1 < num_blocks; ++block_counter1)
    { for (int v_ind = 0; v_ind < v_block_sizes[block_counter1]; ++v_ind) // in each symmetry block, the OBS virtual indices are put before CABS indices
      { const unsigned int F_pA_v_ind =  offset1 + v_ind;
        for (int row_ind = 0; row_ind < no; ++row_ind)
          F_pA.set_element(row_ind, F_pA_v_ind, 0.0);
      }
      offset1 += A_block_sizes[block_counter1];
    }
   }

  // compute the Right-Hand vector
    RefSCMatrix GF = gamma1*F_pA;
    for(int x=0; x<no; ++x)
    {
      for(int B=0; B<nX; ++B)
      {
        rhs_vector->set_element(x * nX + B, -1.0*GF.get_element(x,B));
      }
    }

  //compute trace(f * gamma1) = f^p_q  gamma^q_p; this is an element needed to compute H0
  const double FtG = (F_pp * gamma1).trace() + (F_pp_os * gamma1_otherspin).trace();

  // compute H0; x1/x2 etc denote difference spins; H0 corresponds to the B matrix from the notes

  // first, calculate the intermediate I(x,y) = f^q_p * \gamma^{xp}_{yq};
  for(int x=0; x<no; ++x)
  {
    for(int y=0; y<no; ++y)
    {
      double Ixy_xy = 0.0;
          {
            for (int p = 0; p < no; ++p)
            {
              for (int q = 0; q < no; ++q)
              {// 2-rdm of different spins: arranged in the order of beta-alpha (alpha runs faster)
                const int x1p2 = (spin == Alpha)? (p*no + x):(x*no + p);
                const int y1q2 = (spin == Alpha)? (q*no + y):(y*no + q);
                Ixy_xy += F_pp_os(q, p) * gamma2_ba(x1p2, y1q2);
                if((x != p) && (y != q))  // contribution from gamma2_ss
                {
                  const int u_ss = antisym_pairindex(x,p);
                  const int l_ss = antisym_pairindex(y,q);
                  Ixy_xy += indexsizeorder_sign(x,p) * indexsizeorder_sign(y,q) * F_pp(q,p) * gamma2_ss(u_ss, l_ss);
                }
              }
            }
          }
      Ixy(x,y) = Ixy_xy;
    }
  }

  for(int x=0; x<no; ++x)
  {
    for(int y=0; y<no; ++y)
    {
      const double gamma_xy = gamma1(x,y);
      const double Ixy_xy = Ixy(x, y);
      for(int B=0; B<nX; ++B)
      {
        const int xB = x*nX + B; // calculate the row index
        for(int A=0; A<nX; ++A)
        {
          const int yA = y*nX + A; // calculate the column index
          double h0_xB_yA = 0.0;
          h0_xB_yA += gamma_xy * F_AA(A,B);
          if(A == B) // a Kronecker delta
          {
            h0_xB_yA += -1.0*gamma_xy * FtG + Ixy_xy;
          }
          H0(xB, yA) = h0_xB_yA;
        }
      }
    }
  }

// the old way, kept for testing
//  ExEnv::out0()  << indent << "old solver (comment out when done testing)" << std::endl;
//  H0.solve_lin(rhs_vector);
//  RefSCVector X = rhs_vector.copy();

  // more stable solver
  RefSCVector X = rhs_vector.clone();
  X.assign(0.0);
  const double rcond = lapack_linsolv_symmnondef(H0, X, rhs_vector);
  if (rcond < 1e-8)
    ExEnv::out0() << indent << "[1]_S wfn eqs rcond = " << std::setprecision(12) << rcond << std::endl;

  double E_spin = 0.0;
  for (int j = 0; j < no; ++j)
  {
    for (int A = 0; A < nX; ++A)
    {
      E_spin += GF.get_element(j,A) * X->get_element(j * nX + A);
    }
  }
  return E_spin;
}


double SpinOrbitalPT2R12::cabs_singles_Dyall()
{
  # define DEBUGG false

  const SpinCase1 spin = Alpha;
  Ref<OrbitalSpace> activespace = this->r12world()->refwfn()->occ_act_sb();
  Ref<OrbitalSpace> pspace = rdm1_->orbs(spin); // occupied orbs, determined by rdm1
  Ref<OrbitalSpace> vspace = this->r12world()->refwfn()->uocc_act_sb(spin);
  Ref<OrbitalSpace> cabsspace = this->r12world()->cabs_space(spin);

  Ref<OrbitalSpaceRegistry> oreg = this->r12world()->world()->tfactory()->orbital_registry();
  if (!oreg->value_exists(pspace))
    oreg->add(make_keyspace_pair(pspace));
  const string key = oreg->key(pspace);
  pspace = oreg->value(key);

  Ref<OrbitalSpace> all_virtual_space;
  if (cabs_singles_coupling_)
  {
    all_virtual_space = new OrbitalSpaceUnion("AA", "all virtuals", *vspace, *cabsspace, true);
    if (!oreg->value_exists(all_virtual_space)) oreg->add(make_keyspace_pair(all_virtual_space));
    const string AAkey = oreg->key(all_virtual_space);
    all_virtual_space = oreg->value(AAkey);

    { // make sure that the AO space that supports all_virtual_space is known
      Ref<AOSpaceRegistry> aoreg = this->r12world()->world()->tfactory()->ao_registry();
      if (aoreg->key_exists(all_virtual_space->basis()) == false) {
        Ref<Integral> localints = this->integral()->clone();
        Ref<OrbitalSpace> mu = new AtomicOrbitalSpace("mu''", "CABS(AO)+VBS(AO)", all_virtual_space->basis(), localints);
        oreg->add(make_keyspace_pair(mu));
        aoreg->add(mu->basis(),mu);
      }
    }
  }

  Ref<OrbitalSpace> Aspace;
  if (cabs_singles_coupling_)
    Aspace = all_virtual_space;
  else
    Aspace = cabsspace;

  const unsigned int num_blocks = vspace->nblocks();
  const std::vector<unsigned int>& p_block_sizes = pspace->block_sizes();
  const std::vector<unsigned int>& v_block_sizes = vspace->block_sizes();
  const std::vector<unsigned int>& cabs_block_sizes = cabsspace->block_sizes();
  const std::vector<unsigned int>& A_block_sizes = Aspace->block_sizes();

  // get the fock matrices
  RefSCMatrix F_pA_alpha  = r12eval_->fock(pspace,Aspace,spin);
  RefSCMatrix F_pA_beta   = r12eval_->fock(pspace,Aspace,other(spin));
  RefSCMatrix F_AA_alpha  = r12eval_->fock(Aspace,Aspace,spin);
  RefSCMatrix F_AA_beta   = r12eval_->fock(Aspace,Aspace,other(spin));
  RefSCMatrix F_pp_alpha  = this->f(spin);
  RefSCMatrix F_pp_beta   = this->f(other(spin));

  // RDMs
  RefSymmSCMatrix gamma1_alpha = rdm1_->scmat(spin);
  RefSymmSCMatrix gamma1_beta = rdm1_->scmat(other(spin));
  RefSymmSCMatrix gamma2_aa = this->rdm2(case12(spin,spin));
  RefSymmSCMatrix gamma2_bb = this->rdm2(case12(other(spin),other(spin)));
  RefSymmSCMatrix gamma2_ba = this->rdm2(case12(spin,other(spin)));

#if false
  gamma1_alpha.print("gamma1: alpha");
  gamma1_beta.print("gamma1: beta");
#endif

  // define H0 and necessary vectors
  const int no = pspace->rank();
  const int nv = vspace->rank();
  const int nX = Aspace->rank();
  const int n_cabs = cabsspace->rank();
  const int noX = no * nX;
  RefSCDimension dim_oX = new SCDimension(2 * noX);
  RefSCDimension dim_o  = new SCDimension(2 * no);
  RefSCDimension dim_X  = new SCDimension(2 * nX);
  RefSymmSCMatrix B = gamma2_aa->kit()->symmmatrix(dim_oX);
  B.assign(0.0);
  RefSCMatrix Ixy = gamma2_aa->kit()->matrix(dim_o, dim_o);
  Ixy.assign(0.0);
  RefSCVector rhs_vector = gamma2_aa->kit()->vector(dim_oX);
  rhs_vector->assign(0.0);

#if 0
  ExEnv::out0()  << "  primary, virtual and cabs space dimensions: " << no << ", " << nv << ", " << n_cabs << endl;
  ExEnv::out0()  << "  block dimensions of pspace, vspace, cabsspace: " << endl << "  ";
  { for (int block_counter = 0; block_counter < num_blocks; ++block_counter)
      ExEnv::out0() << scprintf("%5d", p_block_sizes[block_counter]);
    ExEnv::out0() << endl << "  ";
    for (int block_counter = 0; block_counter < num_blocks; ++block_counter)
      ExEnv::out0()  << scprintf("%5d", v_block_sizes[block_counter]);
    ExEnv::out0() << endl << "  ";
    for (int block_counter = 0; block_counter < num_blocks; ++block_counter)
      ExEnv::out0()  << scprintf("%5d", cabs_block_sizes[block_counter]);
    ExEnv::out0() << endl;}
#endif

  if (cabs_singles_coupling_) // when using the combined space, the Fock matrix component f^i_a must be zeroed since it belongs to zeroth order Hamiltonian;
   {
    unsigned int offset1 = 0;
    for (int block_counter1 = 0; block_counter1 < num_blocks; ++block_counter1)
    { // in each symmetry block, the OBS virtual indices are put before CABS indices
      for (int v_ind = 0; v_ind < v_block_sizes[block_counter1]; ++v_ind)
      {
        const unsigned int F_pA_v_ind =  offset1 + v_ind;
        for (int row_ind = 0; row_ind < no; ++row_ind)
        {
          F_pA_alpha.set_element(row_ind, F_pA_v_ind, 0.0);
          F_pA_beta.set_element(row_ind, F_pA_v_ind, 0.0);
        }
      }
      offset1 += A_block_sizes[block_counter1];
    }
   }

  std::vector<int> map_1_to_2 = map(*activespace, *pspace);
  if(!rotate_core_) // if forbit rotating core orbitals; default is 'rotate core orbs'
  {
    for (int row_ind = 0; row_ind < no; ++row_ind)
    {
      if(map_1_to_2[row_ind] < 0)  // row_ind is a core-orbital index
      {
        for (int A1 = 0; A1 < nX; ++A1)
        {
          F_pA_alpha.set_element(row_ind, A1, 0.0);
          F_pA_beta.set_element(row_ind, A1, 0.0);
        }
      }
    }
  }

  // compute the Right-Hand vector
  {
    for(int x1 = 0; x1<no; ++x1)
    {
      for(int B1=0; B1<nX; ++B1)
      {
        double x1B1_alpha = 0.0;
        double x1B1_beta = 0.0;
        for (int j1 = 0; j1 < no; ++j1)
        {
          x1B1_alpha -= gamma1_alpha(x1, j1) * F_pA_alpha(j1, B1);
          x1B1_beta  -= gamma1_beta(x1, j1) * F_pA_beta(j1, B1);
        }
        rhs_vector->set_element(x1 * nX + B1, x1B1_alpha);
        rhs_vector->set_element(x1 * nX + B1 + noX, x1B1_beta);
      }
    }
  }

  // compute H0; x1/x2 etc denote different spins; H0 corresponds to the B matrix from the notes

  { // calculate Ixy: the first two term;
    // the overall minus sign is taken into account when multiplied by delta to calculate B
    int x1, y1, i1,  j1, j2, k1, k2;
    RefSCMatrix g_pppp_ab   = this->g(case12(spin, other(spin)), pspace, pspace, pspace, pspace);
    for(x1=0; x1<no; ++x1)
    {
      for(y1=0; y1 < no; ++y1)
      {
        double I_aa = 0.0;
        double I_bb = 0.0;
        for (i1 = 0; i1 < no; ++i1)  // one contribution to Ixy: f^i1_y1 * gamma^x1_i1
        {
          I_aa += gamma1_alpha.get_element(x1, i1)*F_pp_alpha.get_element(i1,y1);
          I_bb += gamma1_beta.get_element(x1, i1)*F_pp_beta.get_element(i1,y1);
        }
        for (i1 = 0; i1 < no; ++i1)    // v^i1j2_y1k2 * (gamma^x1k2_i1j2 - gamma^x1_i1 * gamma^k2_j2)
        {
            for (j2 = 0; j2 < no; ++j2)
            {
                for (k2 = 0; k2 < no; ++k2)
                {
                    const double v_i1j2y1k2 = g_pppp_ab(i1* no + j2, y1 * no + k2);
//                        const double cumu_x1k2i1j2_aa = gamma2_ba.get_element(x1 * no + k2, i1 * no + j2)
//                            - gamma1_alpha.get_element(x1, i1) * gamma1_beta.get_element(k2, j2);
//                        const double cumu_x1k2i1j2_bb = gamma2_ba.get_element(k2 * no + x1, j2 * no + i1)
//                                                - gamma1_beta.get_element(x1, i1) * gamma1_alpha.get_element(k2, j2);
                    const double cumu_x1k2i1j2_aa = gamma2_ba.get_element(k2*no + x1, j2*no + i1)
                        - gamma1_alpha.get_element(x1, i1) * gamma1_beta.get_element(k2, j2);
                    const double cumu_x1k2i1j2_bb = gamma2_ba.get_element(x1*no + k2, i1*no + j2)
                                            - gamma1_beta.get_element(x1, i1) * gamma1_alpha.get_element(k2, j2);
                    I_aa +=  v_i1j2y1k2 * cumu_x1k2i1j2_aa;
                    I_bb +=  v_i1j2y1k2 * cumu_x1k2i1j2_bb;
                }
            }
        }
        Ixy(x1,y1) = I_aa;
        Ixy(no+x1, no+y1) = I_bb;
      }
    }
  }
//  Ixy.print("Ixy");

  {                                       // calculate Ixy: the last term
      int x1, y1, i1,  j1, j2, k1, k2;
      RefSCMatrix   g_pppp_aa   = this->g(case12(spin, spin), pspace, pspace, pspace, pspace);
      //  RefSCMatrix & g_pppp_bb   = g_pppp_aa;
      for(x1=0; x1<no; ++x1)
      {
          for(y1=0; y1 < no; ++y1)
          {
              double I_aa = 0.0;
              double I_bb = 0.0;
              for (i1 = 0; i1 < no; ++i1)
              {
                  for (j1 = 0; j1 < no; ++j1)
                  {
                      for (k1 = 0; k1 < no; ++k1) // v^i1j1_y1k1 * (0.5 * gamma^x1k1_i1j1 - gamma^x1_i1 * gamma^k1_j1)
                      {
                          if((i1 != j1) && (y1 != k1))
                          {
                            const int g_i1j1y1k1_upp_ind = antisym_pairindex(i1, j1);
                            const int g_i1j1y1k1_low_ind = antisym_pairindex(y1, k1);
                            const double g_i1j1y1k1 = g_pppp_aa.get_element(g_i1j1y1k1_upp_ind, g_i1j1y1k1_low_ind);
                            const int ga2_x1k1i1j1_upp_ind = antisym_pairindex(x1, k1);
                            const int ga2_x1k1i1j1_low_ind = g_i1j1y1k1_upp_ind;
                            double semi_cumu_aa = -1.0 * gamma1_alpha.get_element(x1, i1) * gamma1_alpha.get_element(k1, j1)
                                                  + gamma1_alpha.get_element(x1, j1) * gamma1_alpha.get_element(k1, i1);
                            double semi_cumu_bb = -1.0 * gamma1_beta.get_element(x1, i1)  * gamma1_beta.get_element(k1, j1)
                                                  +gamma1_beta.get_element(x1, j1) * gamma1_beta.get_element(k1, i1);
                            if(x1 != k1)
                            {
                              const double signn = indexsizeorder_sign(x1,k1) * indexsizeorder_sign(i1,j1);
                              semi_cumu_aa += signn * gamma2_aa.get_element(ga2_x1k1i1j1_upp_ind, ga2_x1k1i1j1_low_ind);
                              semi_cumu_bb += signn * gamma2_bb.get_element(ga2_x1k1i1j1_upp_ind, ga2_x1k1i1j1_low_ind);
                            }
                            I_aa +=  0.5 * indexsizeorder_sign(i1, j1) * indexsizeorder_sign(y1, k1) * g_i1j1y1k1 * semi_cumu_aa;
                            I_bb +=  0.5 * indexsizeorder_sign(i1, j1) * indexsizeorder_sign(y1, k1) * g_i1j1y1k1 * semi_cumu_bb;
                          }
                      }
                  }
              }
              Ixy(x1, y1) =  Ixy(x1, y1) + I_aa;
              Ixy(no+ x1, no + y1) =  Ixy(no+ x1, no + y1) + I_bb;
          }
      } // finsh calculating Ixy
  }


  // When dealing with unrelaxed core orbitals, Brillouin condition may not be satisfied;
  // thus we can not assume <| a^core_active H > = 0
  // in this case, our formulas still hold for (x, y) == (active, core);
  // since Ixy should be symmetric, we can calculate Ixy(active, core) for Ixy(core, active)
  // this works whether core orbital are relaxed or not
  for (int row = 0; row < no; ++row)
  {
    for (int col = 0; col < no; ++col) // only the lower (or upper) triangle, otherwise wrong
    {
      if(map_1_to_2[row] < 0 && map_1_to_2[col] > 0) // this means that row <--> core, col<-->active
      {
        Ixy(row, col) = Ixy(col, row);
        Ixy(no + row, no + col) = Ixy(no + col, no + row); // the corresponding beta-beta part
      }
    }
  }


  // explicitly symmetrize Ixy to counteract numerical inaccuracy
  for (int row = 0; row < 2 * no; ++row)
  {
    for (int col = 0; col < row; ++col) // only the lower (or upper) triangle, otherwise wrong
    {
      const double xy = Ixy(row, col);
      const double yx = Ixy(col, row);
      Ixy(row, col) = (xy + yx)/2.0;
      Ixy(col, row) = Ixy(row, col);
    }
  }



// calculate alpha-alpha and beta-beta portion of B: the first three terms
//  f^A1_B1 \gamma^x1_y1 - \delta^A1_B1 * I^x1_y1 + v^A1i2_B1j2 *(gamma^x1j2_y1i2 - \gamma^x1_y1 \gamma^j2_i2)
  {
    RefSCMatrix   g_ApAp_ab   = this->g(case12(spin, other(spin)), Aspace, pspace, Aspace, pspace);
    //  RefSCMatrix & g_pApA_ab   = g_ApAp_ab; // two-electron integrals are the same, as long as the indices are dealt with properly
    unsigned int x1, y1, B1, A1, i2, j2;
    for(x1 = 0; x1<no; ++x1)
    {
      for(y1=0; y1<no; ++y1)
      {
        const double gamma1_x1y1_alpha = gamma1_alpha(x1,y1);
        const double gamma1_x1y1_beta = gamma1_beta(x1, y1);
        const double I_x1y1_alpha = Ixy(x1, y1);
        const double I_x1y1_beta = Ixy(no + x1, no + y1);
        for(B1 = 0; B1 < nX; ++B1)
        {
          const int Baa_row_ind = x1*nX + B1;                            // the row index of  alpha-alpha portion of B matrix
          const int Bbb_row_ind = noX + Baa_row_ind;                     // the row index of  beta-beta portiono of B matrix
          for(A1 = 0; A1 < nX; ++A1)
          {
            double Baa = 0.0;                                            // initialize the B(x1B1, y1A1)
            double Bbb = 0.0;
            const int Baa_col_ind = y1*nX + A1;                          // the column index of alpha-alpha portion of B
            const int Bbb_col_ind = noX + Baa_col_ind;                   // corresponding beta-beta
            Baa += F_AA_alpha(A1, B1) * gamma1_x1y1_alpha;               // f^A1_B1 * gamma^x1_y1; THE MINUS SIGN OF THE INTERMEDIATE IXY IS TAKEN CARE OF HERE
            Bbb += F_AA_beta(A1, B1) * gamma1_x1y1_beta;
            if(A1 == B1)
            {
              Baa += -1.0 * I_x1y1_alpha;                                // -\delta^A1_B1 * I^x1_y1
              Bbb += -1.0 * I_x1y1_beta;
            }
            B(Baa_row_ind, Baa_col_ind) = Baa;
            B(Bbb_row_ind, Bbb_col_ind) = Bbb;
          }
        }
      }
    }
  }

//  B.solve_lin(rhs_vector);
  RefSCVector X = rhs_vector.clone();
  X.assign(0.0);
  const double rcond = lapack_linsolv_symmnondef(B, X, rhs_vector);
  //const double rcond = linsolv_symmnondef_cg(B, X, rhs_vector);
  if (rcond < 1e-8)
    ExEnv::out0() << indent << "[1]_S wfn eqs rcond = " << std::setprecision(12) << rcond << std::endl;

#if false
    B.print("H0");
#endif

  double E_cabs_singles = 0.0;
  for (int i = 0; i < no; ++i)
  {
    for (int j = 0; j < no; ++j)
    {
      const double gamma_ij_alpha = gamma1_alpha(i,j);
      const double gamma_ij_beta = gamma1_beta(i,j);
      for (int A = 0; A < nX; ++A)
      {
        E_cabs_singles += F_pA_alpha(i,A) * gamma_ij_alpha * X->get_element(j * nX + A);
        E_cabs_singles += F_pA_beta(i,A) * gamma_ij_beta * X->get_element(noX + j * nX + A);
      }
    }
  }

  return E_cabs_singles;
}

RefSymmSCMatrix SpinOrbitalPT2R12::density()
{
  throw FeatureNotImplemented("SpinOrbitalPT2R12::density() not yet implemented");
}

void SpinOrbitalPT2R12::brillouin_matrix() {
  RefSCMatrix fmat[NSpinCases1];
  RefSymmSCMatrix opdm[NSpinCases1];
  RefSymmSCMatrix tpcm[NSpinCases2];
  RefSCMatrix g[NSpinCases2];
  Ref<OrbitalSpace> pspace[NSpinCases1];  // all orbitals in OBS
  Ref<OrbitalSpace> mspace[NSpinCases1];  // all occupied orbitals (core + active)
  std::vector<unsigned int> m2p[NSpinCases1];

  Ref<OrbitalSpaceRegistry> oreg = this->r12world()->world()->tfactory()->orbital_registry();

  for(int s=0; s<NSpinCases1; s++) {
    SpinCase1 spin = static_cast<SpinCase1>(s);

    Ref<OrbitalSpace> ospace = rdm1_->orbs(spin);
    if (spin == Alpha) oreg->add(make_keyspace_pair(ospace)); // for some reason this space has not been added

    pspace[s] = this->r12world()->refwfn()->orbs_sb(spin);
    if (!oreg->value_exists(pspace[s])) {
      oreg->add(make_keyspace_pair(pspace[s]));
    }
    {
      const string key = oreg->key(pspace[s]);
      pspace[s] = oreg->value(key);
    }

    mspace[s] = this->r12world()->refwfn()->occ_sb(spin);
    if (!oreg->value_exists(mspace[s])) {
      oreg->add(make_keyspace_pair(mspace[s]));
    }
    {
      const string key = oreg->key(mspace[s]);
      mspace[s] = oreg->value(key);
    }

    m2p[s] = (*pspace[s] << *mspace[s]);

    fmat[s] = r12eval_->fock(pspace[s], pspace[s], spin);
    opdm[s] = rdm1(spin);
  }
  const int nmo = pspace[Alpha]->rank();
  const int nocc = mspace[Alpha]->rank();
  MPQC_ASSERT(mspace[Alpha]->rank() == mspace[Beta]->rank());
  MPQC_ASSERT(pspace[Alpha]->rank() == pspace[Beta]->rank());

  for(int i=0; i<NSpinCases2; i++) {
    SpinCase2 S = static_cast<SpinCase2>(i);
    const SpinCase1 spin1 = case1(S);
    const SpinCase1 spin2 = case2(S);
    const Ref<OrbitalSpace>& space1 = rdm2_->orbs(spin1);
    const Ref<OrbitalSpace>& space2 = rdm2_->orbs(spin2);

    tpcm[i] = lambda2(S);
    g[i] = this->g(S, pspace[spin1], pspace[spin2], space1, space2);
  }

  // precompute intermediates
  RefSCMatrix K[NSpinCases1];   // K = < a_p^q F_N > = \eta f \gamma = f \gamma - \gamma f \gamma
  for(int i=0; i<NSpinCases1; i++) {
    K[i] = fmat[i].kit()->matrix(fmat[i].rowdim(), fmat[i].coldim());   // K is a non-symmetric non-block-diagonal matrix
    K[i].assign(0.0);
    for(int y=0; y<nmo; y++) {
      for(int x=0; x<nocc; x++) {
        const int xx = m2p[i][x];
        for(int z=0; z<nocc; z++) {
          const int zz = m2p[i][z];
          K[i].accumulate_element(y,xx,fmat[i].get_element(y, zz) * opdm[i].get_element(z, x));
        }
      }
    }
#if 1
    for(int y=0; y<nocc; y++) {
      const int yy = m2p[i][y];
      for(int x=0; x<nocc; x++) {
        const int xx = m2p[i][x];
        for(int z1=0; z1<nocc; z1++) {
          const int zz1 = m2p[i][z1];
          for(int z2=0; z2<nocc; z2++) {
            const int zz2 = m2p[i][z2];
            K[i].accumulate_element(yy, xx, -opdm[i].get_element(y, z1) *
                                             fmat[i].get_element(zz1, zz2) *
                                             opdm[i].get_element(z2, x));
          }
        }
      }
    }
#endif
  }

  RefSCMatrix M[NSpinCases1];   // M = < a_p^q W_N >
  M[Alpha] = K[Alpha].clone(); M[Alpha].assign(0.0);
  M[Beta] = K[Beta].clone(); M[Beta].assign(0.0);
  for(int P=0; P<nocc; ++P) {
    for(int Q=0; Q<nmo; ++Q) {
      for(int R=0; R<nocc; ++R) {
        const int RR = m2p[Alpha][R];

        int PR_aa, sign_PR=+1;
        if (P > R) {
          PR_aa = P*(P-1)/2 + R;
        }
        else {
          PR_aa = R*(R-1)/2 + P;
          sign_PR = -1;
        }
        const int PR_ab = P * nocc + R;
        const int RP_ab = R * nocc + P;

        int QR_aa, sign_QR=+1;
        if (Q > RR) {
          QR_aa = Q*(Q-1)/2 + RR;
        }
        else {
          QR_aa = RR*(RR-1)/2 + Q;
          sign_QR = -1;
        }
        const int QR_ab = Q * nmo + RR;
        const int RQ_ab = RR * nmo + Q;

        for (int U = 0; U < nocc; ++U) {
          for (int V = 0; V < nocc; ++V) {
            int UV_aa, sign_UV = +1;
            if (U > V) {
              UV_aa = U * (U - 1) / 2 + V;
            } else {
              UV_aa = V * (V - 1) / 2 + U;
              sign_UV = -1;
            }
            const int UV_ab = U * nocc + V;
            const int VU_ab = V * nocc + U;

            double M_QP_a = 0.0;
            double M_QP_b = 0.0;
            if (P != R && Q != RR && U != V) {
              M_QP_a += 0.5 * sign_PR * sign_QR
                  * g[AlphaAlpha].get_element(QR_aa, UV_aa)
                  * tpcm[AlphaAlpha].get_element(UV_aa, PR_aa);
              M_QP_b += 0.5 * sign_PR * sign_QR
                  * g[BetaBeta].get_element(QR_aa, UV_aa)
                  * tpcm[BetaBeta].get_element(UV_aa, PR_aa);
            }
            M_QP_a += 0.5 * g[AlphaBeta].get_element(QR_ab, UV_ab)
                * tpcm[AlphaBeta].get_element(UV_ab, PR_ab);
            M_QP_a += 0.5 * g[AlphaBeta].get_element(QR_ab, VU_ab)
                * tpcm[AlphaBeta].get_element(VU_ab, PR_ab);
            M_QP_b += 0.5 * g[AlphaBeta].get_element(RQ_ab, VU_ab)
                * tpcm[AlphaBeta].get_element(VU_ab, RP_ab);
            M_QP_b += 0.5 * g[AlphaBeta].get_element(RQ_ab, UV_ab)
                * tpcm[AlphaBeta].get_element(UV_ab, RP_ab);
            M[Alpha].accumulate_element(Q, m2p[Alpha][P], M_QP_a);
            M[Beta].accumulate_element(Q, m2p[Beta][P], M_QP_b);
          }
        }
      }
    }
  }

  for(int s=0; s<NSpinCases1; ++s) {
    const SpinCase1 spin = static_cast<SpinCase1>(s);
    K[s].print(prepend_spincase(spin,"K = eta . f . gamma").c_str());
    M[s].print(prepend_spincase(spin,"M = g . lambda").c_str());
    (K[s] + M[s]).print(prepend_spincase(spin,"BC = K + M").c_str());
    fmat[s].print(prepend_spincase(spin,"f").c_str());
  }
}

void SpinOrbitalPT2R12::print(std::ostream & os) const
{
  os << indent << "SpinOrbitalPT2R12:" << endl;
  os << incindent;
  os << indent << "nfzc = " << nfzc_ << std::endl;
  os << indent << "omit_uocc = " << (omit_uocc_ ? "true" : "false") << std::endl;
  r12world()->print(os);
  Wavefunction::print(os);
  os << decindent;
}

RefSymmSCMatrix SpinOrbitalPT2R12::_rdm2_to_gg(SpinCase2 spin,
                                        RefSymmSCMatrix rdm)
{
  const SpinCase1 spin1 = case1(spin);
  const SpinCase1 spin2 = case2(spin);
  Ref<OrbitalSpace> orbs1 = rdm2_->orbs(spin1);
  Ref<OrbitalSpace> orbs2 = rdm2_->orbs(spin2);
  Ref<OrbitalSpace> gspace1 = r12eval_->ggspace(spin1);
  Ref<OrbitalSpace> gspace2 = r12eval_->ggspace(spin2);
  // if density is already in required spaces?
  if (*orbs1 == *gspace1 && *orbs2 == *gspace2)
    return rdm;

  Ref<LocalSCMatrixKit> local_kit = new LocalSCMatrixKit;
  RefSymmSCMatrix result = local_kit->symmmatrix(r12eval_->dim_gg(spin));
  result.assign(0.0);
  // it's possible for gspace to be a superset of orbs
  std::vector<int> map1 = map(*orbs1, *gspace1);
  std::vector<int> map2 = map(*orbs2, *gspace2);
  SpinMOPairIter UV_iter(gspace1->rank(),gspace2->rank(),spin);
  SpinMOPairIter PQ_iter(gspace1->rank(),gspace2->rank(),spin);
  const int nmo = orbs1->rank();

  for(PQ_iter.start(); int(PQ_iter); PQ_iter.next()) {
    const int P = PQ_iter.i();
    const int Q = PQ_iter.j();
    const int PQ = PQ_iter.ij();
    const int pp = map1[P];
    const int qq = map2[Q];
    if (pp == -1 || qq == -1) continue;   // skip if these indices are not in the source rdm2
    int pq;
    double pfac_pq = 1.0;
    switch(spin) {
      case AlphaBeta: pq = pp * nmo + qq; break;
      case AlphaAlpha:
      case BetaBeta:
        if (pp > qq) {
          pq = (pp * (pp-1)/2 + qq);
        }
        else {
          pq = (qq * (qq-1)/2 + pp);
          pfac_pq = -1.0;
        }
      default: MPQC_ASSERT(false);
    }

    for(UV_iter.start(); int(UV_iter); UV_iter.next()) {
      const int U = UV_iter.i();
      const int V = UV_iter.j();
      const int UV = UV_iter.ij();
      const int uu = map1[U];
      const int vv = map2[V];
      if (uu == -1 || vv == -1) continue;   // skip if these indices are not in the source rdm2
      int uv;
      double pfac_uv = 1.0;
      switch(spin) {
        case AlphaBeta: uv = uu * nmo + vv; break;
        case AlphaAlpha:
        case BetaBeta:
          if (uu > vv) {
            uv = (uu * (uu-1)/2 + vv);
          }
          else {
            uv = (vv * (vv-1)/2 + uu);
            pfac_uv = -1.0;
          }
        default: MPQC_ASSERT(false);
      }

      const double rdm_PQ_UV = pfac_pq * pfac_uv * rdm.get_element(pq, uv);
      result.set_element(PQ, UV, rdm_PQ_UV);
    }
  }

#if 0
  rdm.print(prepend_spincase(spin, "2-rdm (full)").c_str());
  result.print(prepend_spincase(spin, "2-rdm (occ)").c_str());
#endif

  return result;
}

RefSymmSCMatrix SpinOrbitalPT2R12::phi_gg(SpinCase2 spin)
{
  RefSymmSCMatrix phi = this->phi_cumulant(spin);
  return this->_rdm2_to_gg(spin, phi);
}

RefSymmSCMatrix SpinOrbitalPT2R12::lambda2_gg(SpinCase2 spin)
{
  RefSymmSCMatrix lambda = this->lambda2(spin);
  return this->_rdm2_to_gg(spin, lambda);
}

RefSymmSCMatrix SpinOrbitalPT2R12::rdm1_gg(SpinCase1 spin)
{
  RefSymmSCMatrix rdm = this->rdm1(spin);
  Ref<OrbitalSpace> orbs = rdm1_->orbs(spin);
  Ref<OrbitalSpace> gspace = r12eval_->ggspace(spin);
  // if density is already in required spaces?
  if (*orbs == *gspace)
    return rdm;

  Ref<LocalSCMatrixKit> local_kit = new LocalSCMatrixKit;
  RefSymmSCMatrix result = local_kit->symmmatrix(gspace->dim());
  result.assign(0.0);
  // it's possible for gspace to be a superset of orbs
  std::vector<int> omap = map(*orbs, *gspace);
  const int ng = gspace->rank();

  for(int R=0; R<ng; ++R)
    for(int C=0; C<=R; ++C) {
      const int rr = omap[R];
      const int cc = omap[C];
      if (rr == -1 || cc == -1) continue;
      const double rdm_R_C = rdm.get_element(rr, cc);
      result.set_element(R, C, rdm_R_C);
    }

  return result;
}

RefSymmSCMatrix SpinOrbitalPT2R12::rdm2_gg(SpinCase2 spin)
{
  RefSymmSCMatrix rdm = this->rdm2(spin);
  return this->_rdm2_to_gg(spin, rdm);
}

RefSymmSCMatrix SpinOrbitalPT2R12::lambda2(SpinCase2 spin)
{
  // since LocalSCMatrixKit is used everywhere, convert to Local kit
  return convert_to_local_kit(rdm2_->cumulant()->scmat(spin));
}

double SpinOrbitalPT2R12::energy_PT2R12_projector1(SpinCase2 pairspin) {
  const int nelectron = r12world()->refwfn()->nelectron();
  SpinCase1 spin1 = case1(pairspin);
  SpinCase1 spin2 = case2(pairspin);
  Ref<OrbitalSpace> gg1space = r12eval_->ggspace(spin1);
  Ref<OrbitalSpace> gg2space = r12eval_->ggspace(spin2);
  SpinMOPairIter gg_iter(gg1space->rank(),gg2space->rank(),pairspin);

  RefSymmSCMatrix tpdm = rdm2_gg(pairspin);
  RefSymmSCMatrix Phi = phi_gg(pairspin);
  RefSCMatrix VT = V_transformed_by_C(pairspin);
  RefSymmSCMatrix TXT = X_transformed_by_C(pairspin);
  RefSymmSCMatrix TBT = B_transformed_by_C(pairspin);

  ExEnv::out0() << "pairspin "
                << ((pairspin==AlphaBeta) ? "AlphaBeta" : ((pairspin==AlphaAlpha) ? "AlphaAlpha" : "BetaBeta")) << endl;

  // printing pair energies
  RefSCMatrix VT_t_tpdm = 2.0*VT*tpdm;
  RefSCMatrix TBT_t_tpdm = TBT*tpdm;
  RefSCMatrix TXT_t_Phi = TXT*Phi;
  RefSCMatrix HylleraasMatrix = VT_t_tpdm + TBT_t_tpdm - TXT_t_Phi;

  const double energy = this->compute_energy(HylleraasMatrix, pairspin);
  return(energy);
}

int SpinOrbitalPT2R12::nelectron()
{
  return r12world()->refwfn()->nelectron();
}

RefSCMatrix SpinOrbitalPT2R12::V_transformed_by_C(SpinCase2 pairspin) {
  RefSCMatrix T = C(pairspin);
  RefSCMatrix V = r12eval_->V(pairspin);
  RefSCMatrix Vtransposed = V.t();
  RefSCMatrix VT = Vtransposed*T;
  return(VT);
}

RefSymmSCMatrix SpinOrbitalPT2R12::X_transformed_by_C(SpinCase2 pairspin) {
  RefSCMatrix T = C(pairspin);
  RefSymmSCMatrix X = r12eval_->X(pairspin);
  RefSCDimension transformed_dim = T.coldim();
  RefSymmSCMatrix XT = T.kit()->symmmatrix(transformed_dim);
  XT.assign(0.0);
  XT.accumulate_transform(T,X,SCMatrix::TransposeTransform);

  return(XT);
}


RefSymmSCMatrix SpinOrbitalPT2R12::B_transformed_by_C(SpinCase2 pairspin) {
  RefSymmSCMatrix B = r12eval_->B(pairspin);
  RefSCMatrix T = C(pairspin);
  RefSCDimension transformed_dim = T.coldim();
  RefSymmSCMatrix BT = T.kit()->symmmatrix(transformed_dim);
  BT.assign(0.0);
  BT.accumulate_transform(T,B,SCMatrix::TransposeTransform);

  return(BT);
}

RefSymmSCMatrix SpinOrbitalPT2R12::phi_cumulant(SpinCase2 spin12) {
  RefSCMatrix fmat[NSpinCases1];
  RefSymmSCMatrix opdm[NSpinCases1];
  RefSymmSCMatrix tpcm[NSpinCases2];

  for(int s=0; s<NSpinCases1; s++) {
    SpinCase1 spin = static_cast<SpinCase1>(s);
    fmat[s] = f(spin);
    opdm[s] = rdm1(spin);
  }
  const int nmo = opdm[0].dim().n();

  // R12 intermediates, transformed with geminal amplitude matrix
  for(int i=0; i<NSpinCases2; i++) {
    SpinCase2 S = static_cast<SpinCase2>(i);
    tpcm[i] = lambda2(S);
  }

  // precompute intermediates
  RefSCMatrix K[NSpinCases1];   // \gamma f
  RefSCMatrix I[NSpinCases1];   // \gamma f \gamma
  RefSCMatrix M[NSpinCases1];   // trace(f lambda)
  for(int i=0; i<NSpinCases1; i++) {
    K[i] = fmat[i].clone();
    K[i].assign(0.0);
    for(int u=0; u<nmo; u++) {
      for(int z=0; z<nmo; z++) {
        for(int y=0; y<nmo; y++) {
          K[i].accumulate_element(u,z,opdm[i].get_element(u,y)*fmat[i].get_element(y,z));
        }
      }
    }
    I[i] = fmat[i].clone();
    I[i].assign(0.0);
    for(int q=0; q<nmo; q++) {
      for(int v=0; v<nmo; v++) {
        for(int z=0; z<nmo; z++) {
            I[i].accumulate_element(q,v,K[i].get_element(q,z)*opdm[i].get_element(z,v));
        }
      }
    }
  }
#if 0
  for(int s=0; s<NSpinCases1; ++s) {
    K[s] = opdm[s] * fmat[s];
    I[s] = K[s] * opdm[s];
  }
#endif

  M[Alpha] = fmat[Alpha].clone(); M[Alpha].assign(0.0);
  M[Beta] = fmat[Beta].clone(); M[Beta].assign(0.0);
  // each M has contributions from 2 pairspincases
  // AlphaAlpha contributes to Alpha only
  for(int P=0; P<nmo; ++P) {
    for(int Q=0; Q<nmo; ++Q) {
      int PQ_aa, sign_PQ=+1;
      if (P > Q) {
        PQ_aa = P*(P-1)/2 + Q;
      }
      else {
        PQ_aa = Q*(Q-1)/2 + P;
        sign_PQ = -1;
      }
      const int PQ_ab = P * nmo + Q;
      const int QP_ab = Q * nmo + P;
      for(int U=0; U<nmo; ++U) {
        for(int V=0; V<nmo; ++V) {
          int UV_aa, sign_UV=+1;
          if (U > V) {
            UV_aa = U*(U-1)/2 + V;
          }
          else {
            UV_aa = V*(V-1)/2 + U;
            sign_UV = -1;
          }
          const int UV_ab = U * nmo + V;
          const int VU_ab = V * nmo + U;

          double M_PU_a = 0.0;
          double M_PU_b = 0.0;
          if (P!=Q && U!=V) {
            M_PU_a += sign_PQ * sign_UV * fmat[Alpha].get_element(Q, V) * tpcm[AlphaAlpha].get_element(PQ_aa, UV_aa);
            M_PU_b += sign_PQ * sign_UV * fmat[Beta].get_element(Q, V) * tpcm[BetaBeta].get_element(PQ_aa, UV_aa);
          }
          M_PU_a += fmat[Beta].get_element(Q, V) * tpcm[AlphaBeta].get_element(PQ_ab, UV_ab);
          M_PU_b += fmat[Alpha].get_element(Q, V) * tpcm[AlphaBeta].get_element(QP_ab, VU_ab);
          M[Alpha].accumulate_element(P, U, M_PU_a);
          M[Beta].accumulate_element(P, U, M_PU_b);
        }
      }
    }
  }

  if (this->debug_ >= DefaultPrintThresholds::allN2) {
    for(int s=0; s<NSpinCases1; ++s) {
      const SpinCase1 spin = static_cast<SpinCase1>(s);
      K[s].print(prepend_spincase(spin,"K new").c_str());
      I[s].print(prepend_spincase(spin,"I new").c_str());
      M[s].print(prepend_spincase(spin,"M new").c_str());
      fmat[s].print(prepend_spincase(spin,"Fock matrix").c_str());
    }
  }

  // compute phi:
  // \phi^{u v}_{p q} = & P(p q) P(u v) (\gamma^u_p \gamma^{q_3}_q f^{q_2}_{q_3} \gamma^v_{q_2}
  //                                     + \frac{1}{2} \gamma^v_{q_2} f^{q_2}_{q_3} \lambda^{u q_3}_{p q}
  //                                     + \frac{1}{2} \gamma^{q_2}_p f^{q_3}_{q_2} \lambda^{u v}_{q_3 q}
  //                                     - \gamma^u_p f^{q_2}_{q_3} \lambda^{q_3 v}_{q_2 q})
  Ref<LocalSCMatrixKit> local_kit = new LocalSCMatrixKit();
  RefSymmSCMatrix phi = local_kit->symmmatrix(tpcm[spin12].dim());
  phi.assign(0.0);
  const SpinCase1 spin1 = case1(spin12);
  const SpinCase1 spin2 = case2(spin12);
  Ref<OrbitalSpace> orbs1 = rdm2_->orbs(spin1);
  Ref<OrbitalSpace> orbs2 = rdm2_->orbs(spin2);
  SpinMOPairIter UV_iter(orbs1->rank(),orbs2->rank(),spin12);
  SpinMOPairIter PQ_iter(orbs1->rank(),orbs2->rank(),spin12);

  if (spin12 == AlphaBeta) {
    for(PQ_iter.start(); int(PQ_iter); PQ_iter.next()) {
      const int P = PQ_iter.i();
      const int Q = PQ_iter.j();
      const int PQ = PQ_iter.ij();
      const int pp = P;
      const int qq = Q;
      const int pq = pp * nmo + qq;

      for(UV_iter.start(); int(UV_iter); UV_iter.next()) {
        const int U = UV_iter.i();
        const int V = UV_iter.j();
        const int UV = UV_iter.ij();
        const int uu = U;
        const int vv = V;
        const int uv = uu * nmo + vv;

        // first term is easy
        double phi_PQ_UV = (   I[Alpha].get_element(pp, uu) * opdm[Beta].get_element(qq, vv) +
                            opdm[Alpha].get_element(pp, uu) *    I[Beta].get_element(qq, vv)
                           );

        // second and third terms contain contributions from the cumulant of same spin
        for(int q3=0; q3<nmo; ++q3) {
          int uq3 = uu * nmo + q3;
          int q3v = q3 * nmo + vv;
          int pq3 = pp * nmo + q3;
          int q3q = q3 * nmo + qq;
          // +1 instead of +1/2 because there are 2 permutation operators!
          phi_PQ_UV += (tpcm[AlphaBeta].get_element(pq, uq3) * K[Beta].get_element(vv, q3) +
                        tpcm[AlphaBeta].get_element(pq, q3v) * K[Alpha].get_element(uu, q3) +
                        tpcm[AlphaBeta].get_element(pq3, uv) * K[Beta].get_element(qq, q3) +
                        tpcm[AlphaBeta].get_element(q3q, uv) * K[Alpha].get_element(pp, q3)
                       );
        }

        // fourth term is also easy
        phi_PQ_UV -= (   M[Alpha].get_element(pp, uu) * opdm[Beta].get_element(qq, vv) +
                      opdm[Alpha].get_element(pp, uu) *    M[Beta].get_element(qq, vv)
                     );

        phi(PQ,UV) = phi_PQ_UV;
      }
    }
  }
  else if (spin12 == AlphaAlpha || spin12 == BetaBeta) {
    const SpinCase1 spin = spin1;
    for(PQ_iter.start(); int(PQ_iter); PQ_iter.next()) {
      const int P = PQ_iter.i();
      const int Q = PQ_iter.j();
      const int PQ = PQ_iter.ij();
      const int pp = P;
      const int qq = Q;
      MPQC_ASSERT(pp >= qq);
      const int pq = pp * (pp - 1) / 2 + qq;

      for(UV_iter.start(); int(UV_iter); UV_iter.next()) {
        const int U = UV_iter.i();
        const int V = UV_iter.j();
        const int UV = UV_iter.ij();
        const int uu = U;
        const int vv = V;
        MPQC_ASSERT(uu >= vv);
        const int uv = uu * (uu - 1) / 2 + vv;

        // first term is easy
        double phi_PQ_UV = (   I[spin].get_element(pp, uu) * opdm[spin].get_element(qq, vv) +
                            opdm[spin].get_element(pp, uu) *    I[spin].get_element(qq, vv) -
                               I[spin].get_element(pp, vv) * opdm[spin].get_element(qq, uu) -
                            opdm[spin].get_element(pp, vv) *    I[spin].get_element(qq, uu)
                           );

        // second and third terms contain contributions from the cumulant of same spin
        // +1 instead of +1/2 because there are 2 permutation operators!
        // V
        for(int q3=0; q3<uu; ++q3) {
          const int uq3 = uu * (uu-1)/2 + q3;
          phi_PQ_UV += tpcm[spin12].get_element(pq, uq3) * K[spin].get_element(vv, q3);
        }
        for(int q3=uu+1; q3<nmo; ++q3) {
          const int q3u = q3 * (q3-1)/2 + uu;
          // - sign!
          phi_PQ_UV -= tpcm[spin12].get_element(pq, q3u) * K[spin].get_element(vv, q3);
        }
        // U
        for(int q3=0; q3<vv; ++q3) {
          const int vq3 = vv * (vv-1)/2 + q3;
          // - sign!
          phi_PQ_UV -= tpcm[spin12].get_element(pq, vq3) * K[spin].get_element(uu, q3);
        }
        for(int q3=vv+1; q3<nmo; ++q3) {
          const int q3v = q3 * (q3-1)/2 + vv;
          phi_PQ_UV += tpcm[spin12].get_element(pq, q3v) * K[spin].get_element(uu, q3);
        }
        // Q
        for(int q3=0; q3<pp; ++q3) {
          const int pq3 = pp * (pp-1)/2 + q3;
          phi_PQ_UV += tpcm[spin12].get_element(pq3, uv) * K[spin].get_element(qq, q3);
        }
        for(int q3=pp+1; q3<nmo; ++q3) {
          const int q3p = q3 * (q3-1)/2 + pp;
          // - sign!
          phi_PQ_UV -= tpcm[spin12].get_element(q3p, uv) * K[spin].get_element(qq, q3);
        }
        // P
        for(int q3=0; q3<qq; ++q3) {
          const int qq3 = qq * (qq-1)/2 + q3;
          // - sign!
          phi_PQ_UV -= tpcm[spin12].get_element(qq3, uv) * K[spin].get_element(pp, q3);
        }
        for(int q3=qq+1; q3<nmo; ++q3) {
          const int q3q = q3 * (q3-1)/2 + qq;
          // -1 instead of -1/2 because there are 2 permutation operators!
          phi_PQ_UV += tpcm[spin12].get_element(q3q, uv) * K[spin].get_element(pp, q3);
        }

        // fourth term is also easy
        phi_PQ_UV -= (    M[spin].get_element(pp, uu) * opdm[spin].get_element(qq, vv) +
                       opdm[spin].get_element(pp, uu) *    M[spin].get_element(qq, vv) -
                          M[spin].get_element(pp, vv) * opdm[spin].get_element(qq, uu) -
                       opdm[spin].get_element(pp, vv) *    M[spin].get_element(qq, uu)
                     );

        phi(PQ,UV) = phi_PQ_UV;
      }
    }
  }
  else { // invalid spin12
    MPQC_ASSERT(false);
  }

  if(debug_>=DefaultPrintThresholds::mostO4) {
    phi->print(prepend_spincase(spin12,"phi (new)").c_str());
  }

  return phi;
}



RefSCMatrix SpinOrbitalPT2R12::V_genref_projector2(SpinCase2 pairspin) {
  SpinCase1 spin1 = case1(pairspin);
  SpinCase1 spin2 = case2(pairspin);

  const Ref<OrbitalSpace>& GG1_space = r12eval_->GGspace(spin1);
  const Ref<OrbitalSpace>& GG2_space = r12eval_->GGspace(spin2);
  const Ref<OrbitalSpace>& gg1_space = r12eval_->ggspace(spin1);
  const Ref<OrbitalSpace>& gg2_space = r12eval_->ggspace(spin2);

  Ref<SCMatrixKit> localkit = new LocalSCMatrixKit;
  RefSCMatrix V_genref = localkit->matrix(r12eval_->dim_GG(pairspin),r12eval_->dim_gg(pairspin));
  V_genref.assign(0.0);

  const bool antisymmetrize = pairspin!=AlphaBeta;
  const bool part1_equiv_part2 = (GG1_space == GG2_space && gg1_space == gg2_space);

  const RefSCMatrix V_intermed = r12eval_->V(pairspin);
  const RefSymmSCMatrix tpdm = rdm2_gg(pairspin);
  V_genref.accumulate(V_intermed * tpdm);

#if 1
  V_intermed.print(prepend_spincase(pairspin, "V (in SpinOrbitalPT2R12::V_genref_projector2)").c_str());
  tpdm.print(prepend_spincase(pairspin, "gamma(in SpinOrbitalPT2R12::V_genref_projector2)").c_str());
  V_genref.print(prepend_spincase(pairspin, "V.gamma(in SpinOrbitalPT2R12::V_genref_projector2)").c_str());
#endif

  return(V_genref);
}

double SpinOrbitalPT2R12::energy_PT2R12_projector2(SpinCase2 pairspin) {
  SpinCase1 spin1 = case1(pairspin);
  SpinCase1 spin2 = case2(pairspin);
  Ref<OrbitalSpace> gg1space = r12eval_->ggspace(spin1);
  Ref<OrbitalSpace> gg2space = r12eval_->ggspace(spin2);


  RefSymmSCMatrix TBT = B_transformed_by_C(pairspin);
  RefSymmSCMatrix tpdm = rdm2_gg(pairspin);
  RefSCMatrix TBT_tpdm = TBT*tpdm;
  RefSCMatrix HylleraasMatrix = TBT_tpdm;
  B_.push_back(TBT_tpdm.trace());



  RefSCMatrix V_genref = V_genref_projector2(pairspin);
  RefSCMatrix T = C(pairspin);
  RefSCMatrix V_t_T = 2.0*V_genref.t()*T;
  HylleraasMatrix.accumulate(V_t_T);
  V_.push_back(V_t_T.trace());

  RefSymmSCMatrix TXT = X_transformed_by_C(pairspin);
  RefSymmSCMatrix Phi = phi_gg(pairspin);
  RefSCMatrix TXT_t_Phi = TXT*Phi; TXT_t_Phi.scale(-1.0);
  HylleraasMatrix.accumulate(TXT_t_Phi);
  X_.push_back(TXT_t_Phi.trace());

  if (this->debug_ >=  DefaultPrintThresholds::mostO4) {
    V_genref.print(prepend_spincase(pairspin,"Vg").c_str());
    V_t_T.print(prepend_spincase(pairspin,"gVT").c_str());

    TXT.print(prepend_spincase(pairspin,"TXT").c_str());
    TXT_t_Phi.print(prepend_spincase(pairspin,"-TXTf").c_str());

    //rdm2(pairspin).print(prepend_spincase(pairspin,"rdm2 (full)").c_str());
    tpdm.print(prepend_spincase(pairspin,"rdm2").c_str());
    rdm1(Alpha).print(prepend_spincase(Alpha,"rdm1 (full)").c_str());
    rdm1_gg(Alpha).print(prepend_spincase(Alpha,"rdm1").c_str());
    rdm1(Beta).print(prepend_spincase(Beta,"rdm1 (full)").c_str());
    rdm1_gg(Beta).print(prepend_spincase(Beta,"rdm1").c_str());
    TBT.print(prepend_spincase(pairspin,"TBT").c_str());
    TBT_tpdm.print(prepend_spincase(pairspin,"TBTg").c_str());
    const double E_V = V_t_T.trace();
    const double E_X = TXT_t_Phi.trace();
    const double E_B = TBT_tpdm.trace();
  }

  const double energy = this->compute_energy(HylleraasMatrix, pairspin);
  return(energy);
}
double
SpinOrbitalPT2R12::energy_recomputed_from_densities() {
  ExEnv::out0() << std::endl<< indent << "Entered energy_recomputed_from_densities" << std::endl;
  double twoparticle_energy[NSpinCases2];
  double oneparticle_energy[NSpinCases1];
  const int npurespincases2 = spin_polarized() ? 3 : 2;
  const int npurespincases1 = spin_polarized() ? 2 : 1;

  for(int spincase1=0; spincase1<npurespincases1; spincase1++) {
    const SpinCase1 spin = static_cast<SpinCase1>(spincase1);
    const RefSymmSCMatrix H = compute_obints<&Integral::hcore>(rdm1_->orbs(spin));
    RefSymmSCMatrix opdm = rdm1(spin);
    // H and opdm might use different SCMatrixKits -> copy H into another matrix with the same kit as opdm
    RefSymmSCMatrix hh = opdm.clone();
    hh->convert(H);
    hh->element_op(new SCElementScaleDiagonal(0.5));

    Ref<SCElementScalarProduct> trace_op = new SCElementScalarProduct;
    hh->element_op(trace_op, opdm);
    oneparticle_energy[spin] = trace_op->result() * 2.0;
  }

  for(int spincase2=0; spincase2<npurespincases2; spincase2++) {
    const SpinCase2 pairspin = static_cast<SpinCase2>(spincase2);
    const SpinCase1 spin1 = case1(pairspin);
    const SpinCase1 spin2 = case2(pairspin);
    const Ref<OrbitalSpace>& space1 = rdm2_->orbs(spin1);
    const Ref<OrbitalSpace>& space2 = rdm2_->orbs(spin2);

    const RefSymmSCMatrix tpdm = rdm2(pairspin);
    RefSCMatrix G = g(pairspin, space1, space2);
    // G and opdm might use different SCMatrixKits -> copy G into another matrix with the same kit as opdm
    RefSymmSCMatrix gg = tpdm.clone();  gg.assign(0.0);
    gg->accumulate_subblock(G.pointer(), 0, G.nrow()-1, 0, G.ncol() - 1);  // his automatically scales off diagonal elements by 2
                                                                           // no need to scale the diagonal when taking scalar product
    G = 0;

    Ref<SCElementScalarProduct> trace_op = new SCElementScalarProduct;
    gg->element_op(trace_op, tpdm);
    twoparticle_energy[pairspin] = trace_op->result();
  }

  if(!spin_polarized()) {
    twoparticle_energy[BetaBeta] = twoparticle_energy[AlphaAlpha];
    oneparticle_energy[Beta] = oneparticle_energy[Alpha];
  }

#define ENERGY_CONVENTIONAL_PRINT_CONTRIBUTIONS 1

  double energy_hcore = 0.0;
  for(int i=0; i<NSpinCases1; i++) {
    const SpinCase1 spin = static_cast<SpinCase1>(i);
#if ENERGY_CONVENTIONAL_PRINT_CONTRIBUTIONS
    ExEnv::out0() << ((spin==Alpha) ? "Alpha" : "Beta") << " hcore energy: "
                  << setprecision(12) << oneparticle_energy[i] << endl;
#endif
    energy_hcore += oneparticle_energy[i];
  }
#if ENERGY_CONVENTIONAL_PRINT_CONTRIBUTIONS
  ExEnv::out0() << "hcore energy: " << setprecision(12) << energy_hcore << endl;
#endif
  double energy_twoelec = 0.0;
  for(int i=0; i<NSpinCases2; i++) {
#if ENERGY_CONVENTIONAL_PRINT_CONTRIBUTIONS
    string pairspin_str = "";
    if(i==0) {
      pairspin_str = "AlphaBeta";
    }
    else if(i==1) {
      pairspin_str = "AlphaAlpha";
    }
    else if(i==2) {
      pairspin_str = "BetaBeta";
    }
    ExEnv::out0() << "pairspin " << pairspin_str << " energy: " << setprecision(12) << twoparticle_energy[i] << endl;
#endif
    energy_twoelec += twoparticle_energy[i];
  }
#if ENERGY_CONVENTIONAL_PRINT_CONTRIBUTIONS
  ExEnv::out0() << "two-electron energy: " << setprecision(12) << energy_twoelec << endl;
#endif
  double energy = energy_hcore + energy_twoelec;
  energy += this->r12world()->refwfn()->basis()->molecule()->nuclear_repulsion_energy();
  ExEnv::out0() << std::endl<< indent << "Exited energy_recomputed_from_densities" << std::endl << std::endl << std::endl;

  return(energy);
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:

//
// pt2r12.cc
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



#if defined(HAVE_MPQC3_RUNTIME)
#  include <TiledArray/algebra/conjgrad.h>
#endif

using namespace std;
using namespace sc;

#include <chemistry/qc/mbptr12/pt2r12_utils.h>

static ClassDesc PT2R12_cd(typeid(PT2R12),"PT2R12",
                           1,"public Wavefunction",
                           0,
                           create<PT2R12>,
                           create<PT2R12>);

PT2R12::PT2R12(const Ref<KeyVal> &keyval) : Wavefunction(keyval), B_(), X_(), V_()
{
  ///* comment out the following to test cabs singles correction
  string nfzc_str = keyval->stringvalue("nfzc",KeyValValuestring("0"));
  if (nfzc_str == "auto")  nfzc_ = molecule()->n_core_electrons()/2;
  else if (nfzc_str == "no" || nfzc_str == "false") nfzc_ = 0;
  else nfzc_ = atoi(nfzc_str.c_str());

  pt2_correction_ = keyval->booleanvalue("pt2_correction", KeyValValueboolean(true));
  omit_uocc_ = keyval->booleanvalue("omit_uocc", KeyValValueboolean(false));

#if defined(HAVE_MPQC3_RUNTIME)
  cabs_singles_ = keyval->booleanvalue("cabs_singles", KeyValValueboolean(false));
  cabs_singles_h0_ = keyval->stringvalue("cabs_singles_h0", KeyValValuestring(string("fock")));
  cabs_singles_coupling_ = keyval->booleanvalue("cabs_singles_coupling", KeyValValueboolean(true));
  use_mpqc3_ = keyval->booleanvalue("use_mpqc3",KeyValValueboolean(true));
#endif
  rotate_core_ = keyval->booleanvalue("rotate_core", KeyValValueboolean(true));

  rdm2_ = require_dynamic_cast<SpinFreeRDM<Two>*>(
        keyval->describedclassvalue("rdm2").pointer(),
        "PT2R12::PT2R12\n"
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
    throw InputError("PT2R12 requires a WavefunctionWorld; input did not specify it, neither could it be constructed",
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
    if (reference) {
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
      akeyval->assignboolean("spinadapted", true);
      r12world_keyval = new AggregateKeyVal(keyval, akeyval);
    }
    else //doubly check
    {
      bool adapted = keyval->booleanvalue("spinadapted");
      if(not adapted)
        throw InputError("spinadapted must be true for spin-free PT2R12 (the default is correct)",
                         __FILE__, __LINE__, "PT2R12");
    }
    r12world_ = new R12WavefunctionWorld(r12world_keyval, ref);
  }
  r12eval_ = new R12IntEval(r12world_);

  debug_ = keyval->intvalue("debug", KeyValValueint(0));
  r12eval_->debug(debug_);
  // this may update the accuracy of reference_ object
  this->set_desired_value_accuracy(desired_value_accuracy());
#if defined(HAVE_MPQC3_RUNTIME)
  bootup_mpqc3();
  CABS_Single_ = make_shared <CABS_Single> (srr12intrmds_);

#endif
}

PT2R12::PT2R12(StateIn &s) : Wavefunction(s) {
  rdm2_ << SavableState::restore_state(s);
  rdm1_ << SavableState::restore_state(s);
  r12world_ << SavableState::restore_state(s);
  r12eval_ << SavableState::restore_state(s);
  s.get(nfzc_);
  s.get(omit_uocc_);
#if defined(HAVE_MPQC3_RUNTIME)
  s.get(cabs_singles_);
  s.get(cabs_singles_coupling_);
#endif
  s.get(debug_);
}

PT2R12::~PT2R12() {
#if defined(HAVE_MPQC3_RUNTIME)
  shutdown_mpqc3();
#endif
}

void PT2R12::save_data_state(StateOut &s) {
  Wavefunction::save_data_state(s);
  SavableState::save_state(rdm2_, s);
  SavableState::save_state(rdm1_, s);
  SavableState::save_state(r12world_, s);
  SavableState::save_state(r12eval_, s);
  s.put(nfzc_);
  s.put(omit_uocc_);
  s.put(debug_);
}

void
PT2R12::obsolete() {
  r12eval_->obsolete();
  rdm1_->obsolete();
  rdm2_->obsolete();
  r12world_->world()->obsolete();
  r12world_->obsolete();
  Wavefunction::obsolete();
}

void
PT2R12::set_desired_value_accuracy(double acc)
{
  Function::set_desired_value_accuracy(acc);
  Ref<RefWavefunction> refwfn = r12world()->refwfn();
  if (refwfn->desired_value_accuracy_set_to_default()) {
    // reference should be computed to higher accuracy
    const double ref_acc = acc * ref_to_pt2r12_acc();
    refwfn->set_desired_value_accuracy(ref_acc);
  }
}

RefSymmSCMatrix PT2R12::hcore_mo() {
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

RefSCMatrix PT2R12::moints() {
  Ref<SCMatrixKit> localkit = new LocalSCMatrixKit;
  Ref<OrbitalSpace> space1;
  Ref<OrbitalSpace> space2;

  space1 = r12eval_->orbs(AnySpinCase1);
  space2 = r12eval_->orbs(AnySpinCase1);

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

RefSCMatrix PT2R12::C() {
  Ref<LocalSCMatrixKit> local_matrix_kit = new LocalSCMatrixKit();
  RefSCMatrix Cmat = local_matrix_kit->matrix(r12eval_->dim_GG(AlphaBeta),r12eval_->dim_gg(AlphaBeta));
  SpinMOPairIter OW_iter(r12eval_->GGspace(AnySpinCase1)->rank(), r12eval_->GGspace(AnySpinCase1)->rank(), AlphaBeta );
  SpinMOPairIter PQ_iter(r12eval_->ggspace(AnySpinCase1)->rank(), r12eval_->ggspace(AnySpinCase1)->rank(), AlphaBeta );
  CuspConsistentGeminalCoefficient coeff_gen(AlphaBeta);
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
  return(Cmat);
}

RefSCMatrix PT2R12::V_genref_projector2() {
  Ref<SCMatrixKit> localkit = new LocalSCMatrixKit;
  RefSCMatrix V_genref = localkit->matrix(r12eval_->dim_GG(AlphaBeta),r12eval_->dim_gg(AlphaBeta));
  V_genref.assign(0.0);

  Ref<OrbitalSpace> occspace = r12world()->refwfn()->occ_sb(Alpha);
  Ref<OrbitalSpace> gg_space = r12eval_->ggspace(Alpha);

//  ExEnv::out0() << "\n\n" << indent << "Started PT2R12::V_genref_projector2\n";
  RefSCMatrix V_intermed = r12eval_->V_genref_spinfree(occspace, occspace);
  V_intermed.scale(0.5);
  RefSCMatrix tpdm = rdm2_sf_4spaces(occspace, occspace, gg_space, gg_space);
  V_genref.accumulate(V_intermed * tpdm);
//  ExEnv::out0() << "\n\n" <<indent << "Exited PT2R12::V_genref_projector2\n\n";

#if 1
  ExEnv::out0() << __FILE__ << __LINE__ << "\n";
  V_intermed.print(prepend_spincase(AlphaBeta, "V(in PT2R12::V_genref_projector2)").c_str());
  tpdm.print(string("gamma(in PT2R12::V_genref_projector2)").c_str());
  V_genref.print(prepend_spincase(AlphaBeta, "V.gamma(in PT2R12::V_genref_projector2)").c_str());
#endif

  return(V_genref);
}

//
RefSCMatrix PT2R12::B_others() // the terms in B other than B' and X0
{
//  ExEnv::out0() << "\n\n" << indent << "Started PT2R12::B_others\n\n";
  Ref<OrbitalSpace> cabs = r12world_->cabs_space(Alpha);
  Ref<OrbitalSpace> occ_space = r12eval_->occ(Alpha);
  Ref<OrbitalSpace> GG_space = r12eval_->GGspace(Alpha);
  Ref<OrbitalSpace> gg_space = r12eval_->ggspace(Alpha);
//  Ref<OrbitalSpace> pdmspace = rdm1_->orbs(Alpha);
  const int occ_dim = occ_space->rank();
  const int gg_dim = gg_space->rank();

  Ref<DistArray4> wholeproduct; // for the whole contribution

  Ref<DistArray4> RFtimesT; // X^Ax_vw(R^A kappa_tu f^x_kappa * t^vw_tu),
  {

    Ref<OrbitalSpace> F_RI_occ = r12eval_->F_m_P(Alpha);
    R12TwoBodyIntKeyCreator IntCreator(r12eval_->moints_runtime4(),
                                       cabs, GG_space, F_RI_occ, GG_space,
                                       r12eval_->corrfactor(), true);
    std::vector<string> TensorString;
    fill_container(IntCreator, TensorString);
    Ref<TwoBodyMOIntsTransform> RFform = r12eval_->moints_runtime4()->get(TensorString[0]);
    RFform->compute();
    Ref<DistArray4> RF = RFform->ints_distarray4();
    const Ref<TwoBodyIntDescr> & intdescr = RFform->intdescr();
    TwoBodyOper::type f12int_type = r12eval_->corrfactor()->tbint_type_f12();
    const unsigned int f12int_index = intdescr->intset(f12int_type);

    RefSCMatrix T = C().t(); // inverse of T^GG_gg, in principle we should use the transpose matrix
    contract34_DA4_RefMat(RFtimesT, 1.0, RF, f12int_index, T, gg_dim, gg_dim);
  }


  Ref<DistArray4> RTgamma; // the terms in parentheses
  {
    Ref<DistArray4> RT; // X^Ay_rs = R^Ay pq * t^rs_pq
    {
      Ref<OrbitalSpace> cabs = r12world_->cabs_space(Alpha);
      Ref<OrbitalSpace> occ_space = r12eval_->occ(Alpha);
      const int occ_dim = occ_space->rank();
      R12TwoBodyIntKeyCreator IntCreator(r12eval_->moints_runtime4(),
                                         cabs, GG_space, occ_space, GG_space,
                                         r12eval_->corrfactor(), true);
      std::vector<string> TensorString;
      fill_container(IntCreator, TensorString);
      Ref<TwoBodyMOIntsTransform> Rform = r12eval_->moints_runtime4()->get(TensorString[0]);
      Rform->compute();
      Ref<DistArray4> R = Rform->ints_distarray4();
      const Ref<TwoBodyIntDescr> & intdescr = Rform->intdescr();
      TwoBodyOper::type f12int_type = r12eval_->corrfactor()->tbint_type_f12();
      const unsigned int f12int_index = intdescr->intset(f12int_type);
      RefSCMatrix T = C().t();
      contract34_DA4_RefMat(RT, 1.0, R, f12int_index, T, gg_dim, gg_dim);
    }

    RefSCMatrix onerdm = rdm1_gg_sf(); // G^s_v

    {
      Ref<DistArray4> I1; //the first two sets uses contract4
      contract4(RT, onerdm, I1);

      //the first set
      {
        Ref<DistArray4> RTgamma1;
        RefSCMatrix rdm2inter = rdm2_sf_4spaces_int(-0.5, 0.5, -0.5, occ_space, gg_space, gg_space, occ_space);
        rdm2inter = RefSCMAT4_permu<Permute23>(rdm2inter, occ_space, gg_space, gg_space, occ_space);
        rdm2inter = RefSCMAT4_permu<Permute34>(rdm2inter, occ_space, gg_space, gg_space, occ_space);
        Ref<DistArray4> pI1 = permute23(permute34(I1));
        contract34_DA4_RefMat(RTgamma1, 1.0, pI1, 0, rdm2inter, occ_dim, gg_dim);
        RTgamma1 = permute23(RTgamma1);
        RTgamma = RTgamma1;
//        ExEnv::out0() << "\n\n" << indent << indent << "Finished 1st group\n";
      }
      //the second set
      {
        Ref<DistArray4> RTgamma2;
        RefSCMatrix rdm2inter = rdm2_sf_4spaces_int(-0.5, 1.0, -0.25, occ_space, gg_space, occ_space, gg_space);
        rdm2inter = RefSCMAT4_permu<Permute34>(rdm2inter, occ_space, gg_space, occ_space, gg_space);
        rdm2inter = RefSCMAT4_permu<Permute23>(rdm2inter, occ_space, gg_space, gg_space, occ_space);
        rdm2inter = RefSCMAT4_permu<Permute34>(rdm2inter, occ_space, gg_space, gg_space, occ_space);
        Ref<DistArray4> I_Aw_yr = permute23(permute34(I1));
        contract34_DA4_RefMat(RTgamma2, 1.0, I_Aw_yr, 0, rdm2inter, occ_dim, gg_dim);
        RTgamma2 = permute34(permute23(RTgamma2));
        axpy(RTgamma2, 1.0, RTgamma, 1.0);
//        ExEnv::out0() << "\n\n" << indent << indent << "Finished 2nd group\n";
      }
    } // 1st and 2nd sets done

    {
       Ref<DistArray4> I2; //the first two sets uses contract4
       contract3(RT, onerdm, I2);

       //the third
       {
         Ref<DistArray4> RTgamma3;
         RefSCMatrix rdm2inter = rdm2_sf_4spaces_int(1.0, -1.0, 0.0, occ_space, gg_space, gg_space, occ_space);
         rdm2inter = RefSCMAT4_permu<Permute23>(rdm2inter, occ_space, gg_space, gg_space, occ_space);
         rdm2inter = RefSCMAT4_permu<Permute34>(rdm2inter, occ_space, gg_space, gg_space, occ_space);
         contract34_DA4_RefMat(RTgamma3, 1.0, permute23(I2), 0, rdm2inter, occ_dim, gg_dim);
         RTgamma3 = permute23(RTgamma3);
         axpy(RTgamma3, 1.0, RTgamma, 1.0);
//         ExEnv::out0() << "\n\n" << indent << indent << "Finished 3rd group\n";
       }
       //the fourth set
       {
         Ref<DistArray4> RTgamma4;
         RefSCMatrix rdm2inter = rdm2_sf_4spaces_int(-0.5, 0.5, 0.0, occ_space, gg_space, gg_space, occ_space);
         rdm2inter = RefSCMAT4_permu<Permute23>(rdm2inter, occ_space, gg_space, gg_space, occ_space);
         rdm2inter = RefSCMAT4_permu<Permute34>(rdm2inter, occ_space, gg_space, gg_space, occ_space);
         contract_DA4_RefMat_k2b2_34(RTgamma4, 1.0, I2, 0, rdm2inter, occ_dim, gg_dim); //RTgamma is initialized
         RTgamma4 = permute34(RTgamma4);
         axpy(RTgamma4, 1.0, RTgamma, 1.0);
//         ExEnv::out0() << "\n\n" << indent << indent << "Finished 4th group\n";
       }
     } // 3rd and 4th sets done

  } // RTgamma done
#if 0 // check RTgamma; it should be zero for CLHF
  RefSCMatrix RTgammaMat = copy_to_RefSCMat(RTgamma,0);
  RTgammaMat.print(prepend_spincase(AlphaBeta,"RTgammaMat").c_str());
#endif
  contract34(wholeproduct, 1.0, permute23(permute34(permute12(permute23(RFtimesT)))), 0,
                                permute23(permute34(permute12(permute23(RTgamma)))), 0); //(gg, gg)
//  ExEnv::out0() << "\n\n" << indent << "Finished PT2R12::B_others\n\n";
  RefSCMatrix TotalMat = copy_to_RefSCMat(wholeproduct,0);
  return TotalMat;
}



RefSCMatrix PT2R12::X_term_Gamma_F_T() {
  const Ref<OrbitalSpace> & gg_space = r12eval_->ggspace(Alpha);
  const Ref<OrbitalSpace> & GG_space = r12eval_->ggspace(Alpha);
  const int dimg = gg_space->dim()->n();
  const int dimG = GG_space->dim()->n();
  RefSCMatrix T = this->C(); // (dim_GG, dim_gg)

  RefSCMatrix F_gg  = r12eval_->fock(gg_space,gg_space,Alpha);
  RefSCMatrix tpdm =  rdm2_sf_4spaces(gg_space, gg_space, gg_space, gg_space);
  Ref<LocalSCMatrixKit> lmk = new LocalSCMatrixKit();
  RefSCMatrix GammaF = lmk->matrix(r12eval_->dim_gg(AlphaBeta), r12eval_->dim_gg(AlphaBeta));
  GammaF.assign(0.0);
  for (int r = 0; r < dimg; ++r)
  {
    for (int s = 0; s < dimg; ++s)
    {
      const int rowind = r*dimg + s;
      for (int v = 0; v < dimg; ++v)
      {
        for (int w = 0; w < dimg; ++w)
        {
            const int colind = v*dimg + w;
            double rsvw = 0.0;
            for (int x = 0; x < dimg; ++x)
            {
                rsvw += tpdm(rowind, v*dimg + x) * F_gg(x, w);
            }
            GammaF(rowind, colind) = rsvw;
        }
      }
    }
  }
  const bool debug_pp = false;
  if(debug_pp)
  {
    gg_space->basis()->print();
    gg_space->print_detail();
    F_gg.print("debug:F_gg in X_term_Gamma_F_T");
    tpdm.print("debug:tpdm in X_term_Gamma_F_T");
    T.print("debug:T in X_term_Gamma_F_T");
  }
  return GammaF*(T.t());//(gg, GG)   (Gamma^rs_vx * f^x_w) * t^vw_tu
}


double PT2R12::energy_PT2R12_projector2() {

  // 2*V*T constribution
  const bool print_all = true;
  if(print_all)
    ExEnv::out0() << std::endl << std::endl << indent << "Entered PT2R12::energy_PT2R12_projector2\n\n";

  if(print_all)
    ExEnv::out0() << "\n" << indent << "Started V_genref\n";
  RefSCMatrix V_genref = V_genref_projector2(); //(GG, gg)
  if(print_all)
    ExEnv::out0() << "\n\n" << indent << "Finished V_genref\n\n";

  RefSCMatrix T = C();   // C() is of dimension (GG, gg)
  RefSCMatrix V_t_T = 2.0*T.t()*V_genref;
  RefSCMatrix HylleraasMatrix = V_t_T.copy(); // (gg, gg)
  if(print_all)
    ExEnv::out0() << "\n\n" << indent << "Finished V_t_T\n\n";
  if (this->debug_ >=  DefaultPrintThresholds::mostO4 or print_all)
  {
    ExEnv::out0() << __FILE__ << __LINE__ << "\n";
    T.print(string("T").c_str());
    V_genref.print(string("V_genref").c_str());
    HylleraasMatrix.print(string("Hy:+V_t_T").c_str());
    ExEnv::out0() << "E(V_t_T) = " << V_t_T.trace() << std::endl;
  }

  // X contributions
  RefSCMatrix TGFT = X_term_Gamma_F_T(); //(GG, GG)
  if(print_all)
    ExEnv::out0() << "\n\n" << indent << "Finished TGFT\n\n";
  TGFT.scale(-1.0);
  RefSCMatrix Xpart = TGFT * r12eval_->X() * T;
  const bool debug_pp = true;
  if(debug_pp)
  {
    TGFT.print("debug:TGFT");
    r12eval_->X().print("debug:r12eval X");
    T.print("debug:T");
    Xpart.print("Hy:+TGfXT");
  }
  if(print_all)
    ExEnv::out0() << "\n\n" << indent << "Finished Xpart\n\n";
  HylleraasMatrix.accumulate(Xpart);//(gg, gg)
  if (this->debug_ >=  DefaultPrintThresholds::mostO4 or print_all)
  {
    Xpart.print(string("Xpart").c_str());
    HylleraasMatrix.print(string("Hy:+X").c_str());
    ExEnv::out0() << "E(TGfXT) = " << Xpart.trace() << std::endl;
  }


  // B' contribution
  Ref<OrbitalSpace> gg_space = r12eval_->ggspace(Alpha);
  RefSCMatrix TBTG =  T.t() * r12eval_->B() * T  * rdm2_sf_4spaces(gg_space, gg_space, gg_space, gg_space);//(gg,gg)
  if(print_all)
    ExEnv::out0() << "\n\n" << indent << "Finished TBTG\n\n";
  TBTG.scale(0.5);
  HylleraasMatrix.accumulate(TBTG);
  if (this->debug_ >=  DefaultPrintThresholds::mostO4 or print_all)
  {
    TBTG.print(string("TBTG").c_str());
    HylleraasMatrix.print(string("Hy:+TBTG").c_str());
    ExEnv::out0() << "E(TBTG) = " << TBTG.trace() << std::endl;
  }


  // the last messy term
  RefSCMatrix others = B_others(); //(gg, gg)
  HylleraasMatrix.accumulate(others);
  if (this->debug_ >=  DefaultPrintThresholds::mostO4 or print_all)
  {
      others.print(string("others").c_str());
      HylleraasMatrix.print(prepend_spincase(AlphaBeta,"Hy:+others").c_str());
      ExEnv::out0() << "E(others) = " << others.trace() << std::endl;
  }

  if(print_all)
    ExEnv::out0() << std::endl << std::endl << indent << "Exited PT2R12::energy_PT2R12_projector2\n\n";

  bool print_component = false;
  if(print_component)
  {
    const double E_V_t_T = V_t_T.trace();
    const double E_Xpart = Xpart.trace();
    const double E_TBTG = TBTG.trace();
    const double E_others = others.trace();
    const double Btotal =  E_Xpart + E_TBTG + E_others;
    const double E_total = E_V_t_T + Btotal;
    const double VBratio = -E_V_t_T/(2*Btotal);
    const double R12min = - E_V_t_T*E_V_t_T/(4 * Btotal);
    const double deviation = 1- E_total/R12min;
    ExEnv::out0() << std::endl << std::endl;
#define only_include_V false //for debugging
#if only_include_V
      HylleraasMatrix.assign(0.0);
      HylleraasMatrix.accumulate(V_t_T);
      ExEnv::out0() << indent << "Now only V term computed; others zeroed" << std::endl << std::endl;
#endif
    ExEnv::out0() << indent << "individual contributions::" << std::endl;
    ExEnv::out0() << indent << scprintf("V:                        %17.12lf", E_V_t_T) << std::endl;
    ExEnv::out0() << indent << scprintf("B'(0):                    %17.12lf", E_TBTG) << std::endl;
    ExEnv::out0() << indent << scprintf("B'(X):                    %17.12lf", E_Xpart) << std::endl;
    ExEnv::out0() << indent << scprintf("B remain:                 %17.12lf", E_others) << std::endl;
    ExEnv::out0() << indent << scprintf("VBratio remain:           %17.12lf", VBratio) << std::endl << std::endl;
    ExEnv::out0() << indent << scprintf("quadratic min:            %17.12lf", R12min) << std::endl;
    ExEnv::out0() << indent << scprintf("Hylleraas:                %17.12lf", E_total) << std::endl;
    ExEnv::out0() << indent << scprintf("deviation percentage:     %17.12lf", deviation) << std::endl << std::endl << std::endl;
  }

  const double energy = compute_energy(HylleraasMatrix);
  return energy;
}

#if defined(HAVE_MPQC3_RUNTIME)
  void PT2R12::bootup_mpqc3() {
    MPQC_ASSERT(not srr12intrmds_);
    srr12intrmds_ = make_shared<SingleReference_R12Intermediates<double>>(madness::World::get_default(),
        this->r12world());
    srr12intrmds_->set_rdm2(this->rdm2_);
  }

  void PT2R12::shutdown_mpqc3() {
    srr12intrmds_ = 0;
  }

#endif

std::pair<double,double>
PT2R12::energy_PT2R12_projector2_mpqc3() {

#if defined(HAVE_MPQC3_RUNTIME)

 // bootup_mpqc3();

  // see J. Chem. Phys. 135, 214105 (2011) for eqns.

  typedef SingleReference_R12Intermediates<double>::TArray4 TArray4;
  typedef SingleReference_R12Intermediates<double>::TArray2 TArray2;

  const bool print_all = true;
  if(print_all)
    ExEnv::out0() << std::endl << std::endl << indent << "Entered PT2R12::energy_PT2R12_projector2_mpqc3\n\n";

  TArray4 Tg_ij_kl; Tg_ij_kl("i,j,k,l") = _Tg("<i j|Tg|k l>");

  double VT2 = 0.0;
  {
    auto V_ij_mn = V_sf(true);

    // extra factor of 1/2 relative to Eq. (11), but gets cancelled by a factor of 2 as in 2 . V . T
    TArray4 Vg_ij_kl;
    Vg_ij_kl("i1,i2,k1,k2") = 0.5 * _4("<i1 i2|gamma|m1 m2>") * V_ij_mn("k1,k2,m1,m2");

    // cancellation of the previous 1/2 by this 2 to yield Eq. (11)
    VT2 = 2.0 * dot(Vg_ij_kl("i,j,k,l"), Tg_ij_kl("i,j,k,l"));
  }
  madness::World::get_default().gop.fence();
  ExEnv::out0() << indent << "VT2=" << VT2 << std::endl;

  double X = 0.0;
  {
    auto X_ij_kl = X_sf(true);
    TArray4 rdm2_F;
    rdm2_F("i1,i2,j1,j2") = _4("<i1 i2|gamma|j1 m3>") * _2("<m3|F|j2>");
    TArray4 TXT;
    TXT("i1,i2,l1,l2") = Tg_ij_kl("i1,i2,j1,j2") * X_ij_kl("j1,j2,k1,k2") * Tg_ij_kl("k1,k2,l1,l2");
    X = -dot(TXT("i1,i2,j1,j2"), rdm2_F("i1,i2,j1,j2"));
  }
  madness::World::get_default().gop.fence();
  ExEnv::out0() << indent << "X=" << X << std::endl;

  double B0 = 0.0;
  {
    auto B_ij_kl = B_sf(true);
    TArray4 TBT;
    TBT("i1,i2,l1,l2") = Tg_ij_kl("i1,i2,j1,j2") * B_ij_kl("j1,j2,k1,k2") * Tg_ij_kl("k1,k2,l1,l2");
    // extra 1/2 relative to Eq. (12), but B was scaled by factor of 2 relative to that Eq.
    B0 = 0.5 * dot(TBT("i1,i2,j1,j2"), _4("<i1 i2|gamma|j1 j2>"));
  }
  madness::World::get_default().gop.fence();
  ExEnv::out0() << indent << "B0=" << B0 << std::endl;

  double Delta = 0.0;
  {
    TArray4 Trf; Trf("i1,k,a',m") = Tg_ij_kl("i1,k,j1,j2") * _4("<j1 j2|r|a' m_F(p')>");
    TArray4 Tr ;  Tr("l,i2,a',n") = Tg_ij_kl("l,i2,j1,j2") * _4("<j1 j2|r|a' n>");

    TArray2 rdm1_oo;  rdm1_oo("m,n") = _2("<m|gamma|n>");
    TArray2 rdm1_oa;  rdm1_oa("m,i") = _2("<m|gamma|i>");
    TArray2 rdm1_ao;  rdm1_ao("i,m") = _2("<i|gamma|m>");
    TArray2 rdm1_aa;  rdm1_aa("i,j") = _2("<i|gamma|j>");

    {
      TArray4 lambda_1;
      lambda_1("m,k,l,n") = 0.5 * (rdm1_oa("m,k") * rdm1_ao("l,n")
                             - rdm1_oo("m,n") * rdm1_aa("k,l")
                             - _4("<m l|gamma|k n>")
                            );
      TArray4 Trf_gamma_Tr_1;
      Trf_gamma_Tr_1("m,k,l,n") = Tr("l,i2,a',n") * (rdm1_aa("i2,i1") * Trf("i1,k,a',m"));
      Delta += dot(Trf_gamma_Tr_1("m,k,l,n"), lambda_1("m,k,l,n"));
    }
    madness::World::get_default().gop.fence();

    {
      TArray4 lambda_2;
      lambda_2("m,n,l,k")  = (rdm1_oo("m,n") * rdm1_aa("l,k")
                             - 0.25 * rdm1_oa("m,k") * rdm1_ao("l,n")
                             - 0.5 * _4("<m l|gamma|n k>")
                            );
      TArray4 Trf_gamma_Tr_2;
      Trf_gamma_Tr_2("m,n,l,k") = Tr("l,i2,a',n") * (rdm1_aa("i2,i1") * Trf("k,i1,a',m"));
      Delta += dot(Trf_gamma_Tr_2("m,n,l,k"), lambda_2("m,n,l,k"));
    }
    madness::World::get_default().gop.fence();

    {
      TArray4 lambda_3;
      lambda_3("m,l,k,n") = (_4("<m l|gamma|k n>")
                             - rdm1_oa("m,k") * rdm1_ao("l,n")
                            );
      {
        TArray4 Trf_gamma_Tr_3;
        Trf_gamma_Tr_3("m,l,k,n") = Tr("i2,l,a',n") * (rdm1_aa("i2,i1") * Trf("i1,k,a',m"));
        Delta += dot(Trf_gamma_Tr_3("m,l,k,n"), lambda_3("m,l,k,n"));
      }
      madness::World::get_default().gop.fence();

      {
        // lambda_4 = -0.5 lambda_3
        TArray4 Trf_gamma_Tr_4;
        Trf_gamma_Tr_4("m,l,k,n") = Tr("i2,l,a',n") * (rdm1_aa("i2,i1") * Trf("k,i1,a',m"));
        Delta += -0.5 * dot(Trf_gamma_Tr_4("m,l,k,n"), lambda_3("m,l,k,n"));
      }
    }
    madness::World::get_default().gop.fence();

  }
  std::cout << indent << "Delta=" << Delta << std::endl;

  double eref_recomp = 0.0;
  {
    eref_recomp = dot(_2("<m1|h|n1>"), _2("<m1|gamma|n1>")) +
        0.5 * dot(_4("<m1 m2|g|n1 n2>"), _4("<m1 m2|gamma|n1 n2>"));
  }
  eref_recomp += r12world()->refwfn()->basis()->molecule()->nuclear_repulsion_energy();
  madness::World::get_default().gop.fence();

 // shutdown_mpqc3();

  return std::make_pair(VT2 + X + B0 + Delta, eref_recomp);
#else
  throw ProgrammingError("PT2R12::energy_PT2R12_projector2_mpqc3() called but MPQC3 runtime is not available",
                         __FILE__, __LINE__);
  return std::make_pair(0.0, 0.0); // unreachable
#endif
}

double PT2R12::compute_energy(const RefSCMatrix &hmat,
                              bool print_pair_energies,
                              std::ostream& os) {
  Ref<OrbitalSpace> ggspace = r12eval_->ggspace(AnySpinCase1);
  SpinMOPairIter gg_iter(ggspace->rank(),ggspace->rank(),AlphaBeta);
  double energy = 0.0;

  if (print_pair_energies) {
    os << indent << "[2]_R12 pair energies:" << endl;
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
  return energy;
}


RefSCMatrix PT2R12::transform_MO() //transformation matrix between occupied orbitals (this also defines the matrix dimension); row: new MO, column: old MO
{                                       // assume mo_density is of the same ordering as occ_sb()
  const bool debugprint = false;
  RefSymmSCMatrix mo_density =  rdm1();//this will eventually read the checkpoint file. I assume they are of the dimension of occ orb space
  Ref<OrbitalSpace> unscreen_occ_act = r12world()->refwfn()->occ_act_sb();
  Ref<OrbitalSpace> occ = r12world()->refwfn()->occ_sb();
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

RefSymmSCMatrix PT2R12::rdm1_sf_transform()
{
  static bool printed = false;
  RefSymmSCMatrix sf_opdm = rdm1();//converted to local
  RefSCMatrix transMO_nonlocal = this->transform_MO();
  MPQC_ASSERT((sf_opdm->n() == transMO_nonlocal->coldim()->n()) and (sf_opdm->n()==r12eval_->occ(AnySpinCase1)->rank()));//should be the dimension of occ space
  RefSCMatrix transMO = convert_RefSC_to_local_kit(transMO_nonlocal);//sf_opdm is converted to local, so transMO needs to do so too; otherwise aborts.
#if 0
    sf_opdm.print(prepend_spincase(AlphaBeta, "rdm1_sf: rdm before transformation").c_str());
    transMO.print(prepend_spincase(AlphaBeta, "rdm1_sf: transMO").c_str());
#endif
  RefSCMatrix res = (transMO*sf_opdm)*(transMO.t());
  RefSymmSCMatrix final = sf_opdm->clone();
  final->assign(0.0);
  final.copyRefSCMatrix(res);// convert a symmscmatrix

  //print occ numbers in rotates basis
  if(not printed)
  {
    std::vector<double> vec_RDM;
    RefDiagSCMatrix RDMdiag = transMO->kit()->diagmatrix(final->dim());//dimenstion (obs, obs)
    for (int x = 0; x < final->n(); ++x)
     {
       vec_RDM.push_back(final.get_element(x,x));
     }
     std::sort(vec_RDM.begin(), vec_RDM.end(), std::greater_equal<double>());
     for (int y = 0; y < final->n(); ++y)
     {
       RDMdiag->set_element(y, vec_RDM[y]);
     }
     RefDiagSCMatrix RDMratio = RDMdiag.clone();
     RDMratio.assign(0.0);
     const double tot = std::accumulate(vec_RDM.begin(), vec_RDM.end(), 0.0);
     std::vector<double>::iterator itt = vec_RDM.begin();
     itt++;
     for (int y = 0; y < final->n(); ++y, ++itt)
     {
       double sofar = std::accumulate(vec_RDM.begin(), itt, 0.0);
       RDMratio->set_element(y, sofar/tot);
     }
     RDMdiag.print(string("PT2R12 1-RDM: occ number").c_str());
     RDMratio.print(string("PT2R12 1-RDM: accumulated occ number percentage").c_str());
     printed = true;
   }
  return final;
}

RefSymmSCMatrix PT2R12::rdm1_sf()
{
  static bool printed = false;
  RefSymmSCMatrix sf_opdm = rdm1();//converted to local
  {
    return(sf_opdm);
  }
}

RefSCMatrix PT2R12::rdm1_gg_sf()
{
      Ref<OrbitalSpace> ggspace = r12eval_->ggspace(Alpha);
      return rdm1_sf_2spaces(ggspace, ggspace);
}

RefSCMatrix PT2R12::rdm1_sf_2spaces(const Ref<OrbitalSpace> b1space, const Ref<OrbitalSpace> k1space)
{
  const unsigned int b1dim = b1space->rank();
  const unsigned int k1dim = k1space->rank();
  RefSCDimension upp_scdim = new SCDimension(b1dim);
  RefSCDimension low_scdim = new SCDimension(k1dim);

  RefSCMatrix sfrdm1 = rdm1_sf().convert2RefSCMat();
  RefSCMatrix result = sfrdm1->kit()->matrix(upp_scdim, low_scdim);

  // 1. to avoid potential problems and for cleanness, all orb spaces should be accessed through r12eval_->r12world_->ref_->spinspaces_ (or screened_spinspaces_);
  // 2. take care of screening
  Ref<OrbitalSpace> pdmspace;
  pdmspace = get_r12eval()->r12world()->refwfn()->occ_sb();

  const unsigned int n_pdmspace = pdmspace->rank();

  std::vector<int> b1map = map(*pdmspace, *b1space);
  std::vector<int> k1map = map(*pdmspace, *k1space);

  for (int nb1 = 0; nb1 < b1dim; ++nb1)
  {
    if(b1map[nb1] == -1)
      throw ProgrammingError("some orbital in b1space not belong to pdmspace; unexpected.", __FILE__,__LINE__);
  }
  for (int nk1 = 0; nk1 < k1dim; ++nk1)
  {
    if(k1map[nk1] == -1)
      throw ProgrammingError("some orbital in k1space not belong to pdmspace; unexpected.", __FILE__,__LINE__);
  }

  for (int nb1 = 0; nb1 < b1dim; ++nb1)
  {
    for (int nk1 = 0; nk1 < k1dim; ++nk1)
    {
      const double rdmelement = sfrdm1.get_element(b1map[nb1], k1map[nk1]);
      result(nb1, nk1) = rdmelement;
    }
  }
//  result.print(string("rdm1_sf_2spaces").c_str());
  return result;
}

RefSymmSCMatrix PT2R12::rdm2_sf()
{
  RefSymmSCMatrix sf_rdm = rdm2();

    return sf_rdm;
}


RefSCMatrix PT2R12::rdm2_sf_4spaces(const Ref<OrbitalSpace> b1space, const Ref<OrbitalSpace> b2space,const Ref<OrbitalSpace> k1space, const Ref<OrbitalSpace> k2space)
{
  const unsigned int b1dim = b1space->rank();
  const unsigned int b2dim = b2space->rank();
  const unsigned int k1dim = k1space->rank();
  const unsigned int k2dim = k2space->rank();
  const unsigned int upp_dim = b1dim * b2dim;
  const unsigned int low_dim = k1dim * k2dim;
  RefSCDimension upp_scdim = new SCDimension(upp_dim);
  RefSCDimension low_scdim = new SCDimension(low_dim);

  RefSCMatrix sf2rdm = rdm2_sf().convert2RefSCMat();
  RefSCMatrix result = sf2rdm->kit()->matrix(upp_scdim, low_scdim);
  result.assign(0.0);

  // 1. to avoid potential problems and for cleanness, all orb spaces should be accessed through r12eval_->r12world_->ref_->spinspaces_ (or screened_spinspaces_);
  // 2. take care of screening
  Ref<OrbitalSpace> pdmspace = get_r12eval()->r12world()->refwfn()->occ_sb();

  const unsigned int n_pdmspace = pdmspace->rank();

#if 0
  if(fabs(r12world_->refwfn()->occ_thres()) > PT2R12::zero_occupancy())
  {
    // first test 1rdm, which shall be diagonal in our test case
    RefSymmSCMatrix ABX = rdm1_sf();
    ExEnv::out0() << __FILE__ << ": " << __LINE__ << "\n";
    ExEnv::out0() << "ordm dimension: " << ABX->n() << "\n";
    ABX.print(prepend_spincase(AlphaBeta, "ordm").c_str());

    RefSCMatrix old_mocoef2 = r12world()->refwfn()->orbs_sb()->coefs();
    ExEnv::out0() << __FILE__ << ": " << __LINE__ << "\n";
    ExEnv::out0() << "print old symmetry-block Mo coefs of dimension: " << old_mocoef2->ncol() << " (AO dimensin: )" << old_mocoef2->nrow() << "\n";
    old_mocoef2.print(prepend_spincase(AlphaBeta, "symm old mo").c_str());

    RefSCMatrix new_mocoef2 = r12world()->refwfn()->orbs_sb()->coefs();
    ExEnv::out0() << __FILE__ << ": " << __LINE__ << "\n";
    ExEnv::out0() << "print new symmetry-block Mo coefs of dimension: " << new_mocoef2->ncol() << " (AO dimensin: )" << new_mocoef2->nrow() << "\n";
    new_mocoef2.print(prepend_spincase(AlphaBeta, "symm new mo").c_str());

    ExEnv::out0() << __FILE__ << ": " << __LINE__ << "\n";
    ExEnv::out0() << "print b1space of dimension: " << b1space->coefs()->ncol() << " (AO dimensin: )" << b1space->coefs()->nrow() << "\n";
    b1space->coefs().print(prepend_spincase(AlphaBeta, "b1space").c_str());
    ExEnv::out0() << "\n\nprint pdmspace of dimension: " << pdmspace->coefs()->ncol() << " (AO dimensin: )" << pdmspace->coefs()->nrow() << "\n";
    pdmspace->coefs().print(prepend_spincase(AlphaBeta, "pdmspace").c_str());
  }
#endif
  std::vector<int> b1map = map(*pdmspace, *b1space);
  std::vector<int> b2map = map(*pdmspace, *b2space);
  std::vector<int> k1map = map(*pdmspace, *k1space);
  std::vector<int> k2map = map(*pdmspace, *k2space);

  for (int nb1 = 0; nb1 < b1dim; ++nb1)
  {
    if(b1map[nb1] == -1)
      {abort(); // debug now
      throw ProgrammingError("some orbital in b1space not belong to pdmspace; unexpected.", __FILE__,__LINE__);}
  }
  for (int nb2 = 0; nb2 < b2dim; ++nb2)
  {
    if(b2map[nb2] == -1)
      {abort();throw ProgrammingError("some orbital in b2space not belong to pdmspace; unexpected.", __FILE__,__LINE__);}
  }
  for (int nk1 = 0; nk1 < k1dim; ++nk1)
  {
    if(k1map[nk1] == -1)
      {abort();throw ProgrammingError("some orbital in k1space not belong to pdmspace; unexpected.", __FILE__,__LINE__);}
  }
  for (int nk2 = 0; nk2 < k2dim; ++nk2)
  {
    if(k2map[nk2] == -1)
      {abort();throw ProgrammingError("some orbital in k2space not belong to pdmspace; unexpected.", __FILE__,__LINE__);}
  }

  Ref<SpinMOPairIter> upp_pair = new SpinMOPairIter(b1dim, b2dim, AlphaBeta);
  Ref<SpinMOPairIter> low_pair = new SpinMOPairIter(k1dim, k2dim, AlphaBeta);
  for(upp_pair->start(); *upp_pair; upp_pair->next())   // add BetaAlpha/AlphaAlpha/BetaBeta contribution to sf_rdm
  {
    for(low_pair->start(); *low_pair; low_pair->next())
     {
         const int nb1 = upp_pair->i();
         const int nb2 = upp_pair->j();
         const int nk1 = low_pair->i();
         const int nk2 = low_pair->j();
         double rdmelement = sf2rdm.get_element(b1map[nb1] * n_pdmspace + b2map[nb2], k1map[nk1] * n_pdmspace + k2map[nk2]);
         result(upp_pair->ij(), low_pair->ij()) = rdmelement;
     }
  }

  return result;
}



//'int' means intermediate
RefSCMatrix PT2R12::rdm2_sf_4spaces_int(const double a, const double b, double const c,
                                                  const Ref<OrbitalSpace> b1space,
                                                  const Ref<OrbitalSpace> b2space,
                                                  const Ref<OrbitalSpace> k1space,
                                                  const Ref<OrbitalSpace> k2space)
{
  const unsigned int nb1 = b1space->rank();
  const unsigned int nb2 = b2space->rank();
  const unsigned int nk1 = k1space->rank();
  const unsigned int nk2 = k2space->rank();

  RefSCMatrix rdm2 = rdm2_sf_4spaces(b1space, b2space, k1space, k2space);
  RefSCMatrix result = rdm2.clone();
  result.assign(0.0);
  rdm2.scale(a);
  RefSCMatrix opdm1a = rdm1_sf_2spaces(b1space, k1space);
  RefSCMatrix opdm1b = rdm1_sf_2spaces(b2space, k2space);
  RefSCMatrix opdm2a = rdm1_sf_2spaces(b1space, k2space);
  RefSCMatrix opdm2b = rdm1_sf_2spaces(b2space, k1space);
  Ref<OrbitalSpace> occ_space = r12eval_->ggspace(Alpha);
//  ExEnv::out0() << indent << occ_space->rank() << std::endl;
  Ref<SpinMOPairIter> upp_pair = new SpinMOPairIter(b1space->rank(), b2space->rank(), AlphaBeta);
  Ref<SpinMOPairIter> low_pair = new SpinMOPairIter(k1space->rank(), k2space->rank(), AlphaBeta);
  for(upp_pair->start(); *upp_pair; upp_pair->next())
  {
    for(low_pair->start(); *low_pair; low_pair->next())
    {
      const unsigned int b1 = upp_pair->i();
      const unsigned int b2 = upp_pair->j();
      const unsigned int k1 = low_pair->i();
      const unsigned int k2 = low_pair->j();
      const double element = rdm2.get_element(upp_pair->ij(), low_pair->ij())
               + b * opdm1a.get_element(b1,k1) * opdm1b.get_element(b2,k2)
               + c * opdm2a.get_element(b1,k2) * opdm2b.get_element(b2,k1);
      result.set_element(upp_pair->ij(), low_pair->ij(), element);
    }
  }
  return result;
}


template<PT2R12::Tensor4_Permute HowPermute>
RefSCMatrix PT2R12::RefSCMAT4_permu(RefSCMatrix rdm2_4space_int,
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

RefSymmSCMatrix PT2R12::rdm1()
{
  return convert_to_local_kit(rdm1_->scmat());
}



RefSymmSCMatrix PT2R12::rdm2()
{
  // since LocalSCMatrixKit is used everywhere, convert to Local kit
  return convert_to_local_kit(rdm2_->scmat());
}



void PT2R12::compute()
{
  r12world()->initialize();
  const bool debug_printing = false;
  if(debug_printing)
  {
//    basis->print();
    rdm1_sf().print("debug: rdm1_sf");
    //rdm2_sf().print("debug: rdm2_sf");
    r12eval_->ggspace(Alpha)->print_detail();
    r12eval_->occ(Alpha)->print_detail();
    r12eval_->vir(Alpha)->print_detail();
    r12eval_->orbs(Alpha)->print_detail();
    r12world()->cabs_space(Alpha)->print_detail();
    r12eval_->ordm(Alpha)->print("ordm-alpha");
    r12eval_->ordm()->print("ordm");
  }
  r12world()->cabs_space(Alpha); // make sure that CABS spaces are computed

  const double energy_ref = r12world()->refwfn()->energy();
  double energy_correction_r12 = 0.0;
  double energy_pt2r12;
  double B_so = 0.0;
  double V_so = 0.0;
  double X_so = 0.0; // these 3 varialbes store the final values for B, V, X
  double energy_pt2r12_sf = 0.0;
  double cabs_singles_e = 0.0;
  const bool spin_polarized = r12world()->refwfn()->spin_polarized();

#if 0
  if(calc_davidson_)
  {
    if(not r12world()->refwfn()->force_rasscf())
      throw ProgrammingError("To calc Davidson correction, must set force_correlate_rasscf_ to true to enable proper initialization of conventional virtual orbital space.", __FILE__,__LINE__);
    Ref<OrbitalSpace> conv_vir = r12world()->refwfn()->conv_uocc_sb();
    Ref<OrbitalSpace> conv_occ = r12world()->refwfn()->conv_occ_sb();
    RefSCMatrix opdm_vv = rdm1_sf_2spaces(conv_vir, conv_vir);
    RefSCMatrix tpdm_vvvv = rdm2_sf_4spaces(conv_vir, conv_vir, conv_vir, conv_vir);
    const double vir_percentage = opdm_vv->trace()-0.5 * tpdm_vvvv->trace();
    ExEnv::out0() << std::endl << std::endl << indent <<scprintf("unoccupied occ num:                     %17.12lf",vir_percentage) << "\n";
    const double davidson = 1/(1 - vir_percentage) - 1;
    ExEnv::out0()  << indent << scprintf("Davidson correction coef:              %17.12lf",
                                            davidson) << std::endl << std::endl;
  }
#endif

  double recomp_ref_energy = 0.0;
  if (pt2_correction_)
  {
    MPQC_ASSERT(r12world()->r12tech()->ansatz()->projector() == R12Technology::Projector_2);

#if defined(HAVE_MPQC3_RUNTIME)
    if (use_mpqc3_) {
      std::pair<double,double> e = energy_PT2R12_projector2_mpqc3();
      energy_pt2r12_sf = e.first;
      recomp_ref_energy = e.second;
    } else
#endif
    {
      energy_pt2r12_sf = energy_PT2R12_projector2();
      recomp_ref_energy = this->energy_recomputed_from_densities();
    }
    energy_correction_r12 = energy_pt2r12_sf;
  }

#if defined(HAVE_MPQC3_RUNTIME)
  if(cabs_singles_ && use_mpqc3_)
  {
    cabs_singles_e = CABS_Single_->compute(cabs_singles_h0_);
  }
#endif

  const double energy = energy_ref + energy_correction_r12 + cabs_singles_e;



    ExEnv::out0() <<endl << indent << scprintf("Reference energy [au]:                 %17.12lf",
                                       energy_ref) << std::endl << std::endl;
#if defined(HAVE_MPQC3_RUNTIME)
    if(cabs_singles_)
    {
      std::string es = "CABS singles(" + cabs_singles_h0_ + ")";
      const unsigned int LL = std::string("Reference energy [au]:                 ").size();
      es.resize(LL, ' ');
      ExEnv::out0() << indent << scprintf((es + "%17.12lf").c_str(),  cabs_singles_e) << endl;
      ExEnv::out0() << indent << scprintf("RASSCF+CABS singles:                   %17.12lf",
                                                energy_ref + cabs_singles_e) << endl << endl;
    }
#endif

    ExEnv::out0() << std::endl << std::endl << indent << scprintf("Reference energy (%9s) [au]:     %17.12lf",
                                        (this->r12world()->world()->basis_df().null() ? "   recomp" : "recomp+DF"),
                                        recomp_ref_energy) << endl;

  if(pt2_correction_)
  {
    ExEnv::out0() << indent << scprintf("[2]_R12 energy [au]:                   %17.12lf",
                                        energy_correction_r12) << endl;
    ExEnv::out0() << indent << scprintf("Total [2]_R12 energy [au]:             %17.12lf",
                                        energy) << std::endl;
  }

  set_energy(energy);
}

double PT2R12::magnetic_moment() const
{
  return r12world()->refwfn()->magnetic_moment();
}
/*
double PT2R12::cabs_singles_Complete()
{
# define DEBUGG false

  const SpinCase1 spin = Alpha;
  Ref<OrbitalSpace> pspace = this->r12world()->refwfn()->occ_sb();
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

  // block size
  const unsigned int num_blocks = vspace->nblocks();
  const std::vector<unsigned int>& p_block_sizes = pspace->block_sizes();
  const std::vector<unsigned int>& v_block_sizes = vspace->block_sizes();
  const std::vector<unsigned int>& cabs_block_sizes = cabsspace->block_sizes();
  const std::vector<unsigned int>& A_block_sizes = Aspace->block_sizes();

  Ref<LocalSCMatrixKit> local_kit = new LocalSCMatrixKit;

  // dimension
  const int nA = Aspace->rank();
  const int ni = pspace->rank();
  RefSCDimension dimAA = new SCDimension(nA*nA);
  RefSCDimension dimii = new SCDimension(ni*ni);
  RefSCDimension dimiA = new SCDimension(ni*nA);
  RefSCDimension dimH = new SCDimension(ni*nA+1);
  RefSCDimension dimi = new SCDimension(ni);
  RefSCDimension dimA = new SCDimension(nA);
  RefSCVector vec_AA = local_kit->vector(dimAA);
  RefSCVector vec_ii = local_kit->vector(dimii);

  // get one-electron operator
  RefSCMatrix hcore_ii_block = r12eval_->fock(pspace, pspace, spin, 0.0, 0.0);
  RefSCMatrix hcore_AA_block = r12eval_->fock(Aspace, Aspace, spin, 0.0, 0.0);
  RefSCMatrix hcore_iA_block = r12eval_->fock(pspace, Aspace, spin, 0.0, 0.0);
  RefSymmSCMatrix hcore_AA = local_kit->symmmatrix(dimA);
  RefSymmSCMatrix hcore_ii = local_kit->symmmatrix(dimi);
  RefSCMatrix hcore_iA = local_kit->matrix(dimi, dimA);
 for (int aa = 0; aa < nA; ++aa) // can't use accumulate
  {
    for (int bb = 0; bb < nA; ++bb)
    {
      hcore_AA->set_element(aa,bb, hcore_AA_block->get_element(aa,bb));
    }
  }
  for (int ii = 0; ii < ni; ++ii)
  {
    for (int jj = 0; jj < ni; ++jj)
    {
      hcore_ii->set_element(ii,jj, hcore_ii_block->get_element(ii,jj));
    }
  }
  for (int ii = 0; ii < ni; ++ii)
  {
    for (int aa = 0; aa < nA; ++aa)
    {
      hcore_iA->set_element(ii,aa, hcore_iA_block->get_element(ii,aa));
    }
  }
#if DEBUGG
  hcore_AA.print(string("core AA").c_str());
  hcore_ii.print(string("core ii").c_str());
  hcore_iA.print(string("core iA").c_str());
#endif
  RefDiagSCMatrix delta_AA = hcore_AA.kit()->diagmatrix(hcore_AA.dim());
  delta_AA->assign(0.0);
  for (int i = 0; i < nA; ++i)
  {
    delta_AA->set_element(i, 1.0);
  }

  RefSCMatrix gamma1 = rdm1_sf_2spaces(pspace, pspace);
  gamma1.print(string("CABS singles: gamma1").c_str());
  RefSCMatrix gamma2 = rdm2_sf_4spaces(pspace, pspace, pspace, pspace);
  RefSCMatrix B_bar = local_kit->matrix(dimAA, dimii); // intermediate mat
//  RefSCMatrix B = local_kit->matrix(dimiA, dimiA);
  B_bar->assign(0.0);
  RefSCMatrix b_bar = local_kit->matrix(dimi, dimA);
  b_bar->assign(0.0);
  RefSCVector b = local_kit->vector(dimiA); // RHS
  RefSCVector b0 = local_kit->vector(dimH); // RHS


  // Compute B
  // Term1: h(alpha, beta)* gamma(i,j)
  {
    matrix_to_vector(vec_AA, hcore_AA);
    matrix_to_vector(vec_ii, gamma1);
    B_bar->accumulate_outer_product(vec_AA, vec_ii);
  }

  // Term2: -delta(alpha, beta)*(g(im,kl) * gamma(jm,kl) )
  {
    RefSCMatrix gbar_imkl = g(pspace, pspace, pspace, pspace);
    RefSCMatrix g_i_mkl = RefSCMAT_combine234(gbar_imkl, ni, ni, ni, ni);
    RefSCMatrix dbar_mkl_j = RefSCMAT_combine234(gamma2, ni, ni, ni, ni).t();
    RefSCMatrix gd = g_i_mkl * dbar_mkl_j;
    gd->scale(-1.0);
    matrix_to_vector(vec_AA, delta_AA);
    matrix_to_vector(vec_ii, gd);
    B_bar->accumulate_outer_product(vec_AA, vec_ii);
  }

  // Term3: -delta(alpha, beta)*( h(i,k) * gamma(k,j) )
  {
    RefSCMatrix hd = hcore_ii*gamma1;
    hd->scale(-1.0);
    matrix_to_vector(vec_AA, delta_AA);
    matrix_to_vector(vec_ii,hd);
    B_bar->accumulate_outer_product(vec_AA, vec_ii);
  }

  // Term4: +g(alpha k, beta l) *gamma(ki, lj)
  {
    RefSCMatrix gg = g(Aspace, pspace, Aspace, pspace);
    RefSCMatrix permu_gg = RefSCMAT4_permu<Permute23>(gg, Aspace, pspace, Aspace, pspace);
    RefSCMatrix permu_d = RefSCMAT4_permu<Permute23>(gamma2, pspace, pspace, pspace, pspace);
    B_bar->accumulate(permu_gg*permu_d);
  }

  // Term5: +g(alpha k, l beta) *gamma(ki, jl)
  {
    RefSCMatrix gg = g(Aspace, pspace, pspace, Aspace);
    RefSCMatrix permu_gg1 = RefSCMAT4_permu<Permute34>(gg, Aspace, pspace, pspace, Aspace);
    RefSCMatrix permu_gg2 = RefSCMAT4_permu<Permute23>(permu_gg1, Aspace, pspace, Aspace, pspace);
    RefSCMatrix permu_d1 = RefSCMAT4_permu<Permute34>(gamma2, pspace, pspace, pspace, pspace);
    RefSCMatrix permu_d2 = RefSCMAT4_permu<Permute23>(gamma2, pspace, pspace, pspace, pspace);
    B_bar->accumulate(permu_gg2*permu_d2);
  } // finish computing B

  //compute b_bar
  // - gamma(j,k) * h(k, beta) - gamma(jm,kl)*g(beta m, kl)
  {
     b_bar->accumulate(gamma1* hcore_iA);
#if DEBUGG
     gamma1.print(string("gamma1").c_str());
     hcore_iA.print(string("hcore iA").c_str());
     b_bar.print(string("b_bar term1").c_str());
#endif
     RefSCMatrix dd = RefSCMAT_combine234(gamma2, ni, ni, ni, ni);
     RefSCMatrix gg1 = g(Aspace, pspace, pspace, pspace);
     RefSCMatrix gg2 = RefSCMAT_combine234(gg1, nA, ni, ni, ni).t();
     b_bar->accumulate(dd*gg2);
#if DEBUGG
     b_bar.print(string("b_bar +term2").c_str());
#endif
     b_bar->scale(-1.0);
  }
#if DEBUGG
  b_bar.print(string("b_bar before zero").c_str());
#endif

  if (cabs_singles_coupling_) // zero Fock matrix component f^i_a
  {
    unsigned int offset1 = 0;
    int block_counter1, row_ind, v_ind;
    for (block_counter1 = 0; block_counter1 < num_blocks; ++block_counter1)
    {
      for (v_ind = 0; v_ind < v_block_sizes[block_counter1]; ++v_ind)
      {
        const unsigned int b_v_ind =  offset1 + v_ind;
        for (row_ind = 0; row_ind < ni; ++row_ind)
          b_bar.set_element(row_ind, b_v_ind, 0.0);
      }
      offset1 += A_block_sizes[block_counter1];
    }
  } // finish RHS construction

#if DEBUGG
  b_bar.print(string("b_bar after zero").c_str());
#endif
  matrix_to_vector(b, b_bar);

  RefSCMatrix B = RefSCMAT4_permu<Permute14>(B_bar, Aspace, Aspace, pspace, pspace);
  // symmetrize B
  B.accumulate(B.t());
  B.scale(0.5);

  double E;
  if(cabs_singles_h0_ == string("CI"))
  {
    RefSCMatrix H = local_kit->matrix(dimH, dimH);
    H.assign(0.0);
    H.assign_subblock(B, 1, ni*nA, 1, ni*nA);
#if false
    B.print(string("B").c_str());
    H.print(string("H").c_str());
#endif
    for(int i = 1; i < dimH.n(); ++i)
    {
      b0.set_element(i, b.get_element(i-1));
    }
    H.assign_column(b0, 0);
    H.assign_row(b0,0);
#if false
    H.print(string("H with b").c_str());
#endif
    RefSCMatrix U = H.clone();
    RefSCMatrix V = H.clone();
    RefDiagSCMatrix S = local_kit->diagmatrix(dimH);
    H.svd(U,S,V);
//    S.print(string("SVD").c_str());
    E = -1.0 * S.get_element(ni*nA); // seems the roots are sorted in descending order
  }
  else if(cabs_singles_h0_ == string("complete"))
  {
//    RefSCVector bcopy = b->copy();
    RefSCVector X = b->clone();
    X.assign(0.0);
  #if DEBUGG
    b.print(ExEnv::out0());
  #endif
    RefSymmSCMatrix Bsymm = B.kit()->symmmatrix(dimiA);
    Bsymm.assign_subblock(B, 0, ni*nA-1, 0, ni*nA-1);
    const double rcond = lapack_linsolv_symmnondef(Bsymm, X, b);
    if (rcond < 1e-8)
      ExEnv::out0() << indent << "[1]_S wfn eqs rcond = " << std::setprecision(12) << rcond << std::endl;
//    B.solve_lin(b);
  #if DEBUGG
    b.print(string("b").c_str());
  #endif
    E = -1.0 * (X.dot(b));
  }
  return E;
}
*/

RefSymmSCMatrix PT2R12::density()
{
  throw FeatureNotImplemented("PT2R12::density() not yet implemented");
}

void PT2R12::print(std::ostream & os) const
{
  os << indent << "PT2R12:" << endl;
  os << incindent;
  os << indent << "nfzc = " << nfzc_ << std::endl;
  os << indent << "omit_uocc = " << (omit_uocc_ ? "true" : "false") << std::endl;
  r12world()->print(os);
  Wavefunction::print(os);
  os << decindent;
}

RefSymmSCMatrix PT2R12::_rdm2_to_gg(RefSymmSCMatrix rdm)
{
  Ref<OrbitalSpace> orbs = rdm2_->orbs();
  Ref<OrbitalSpace> gspace = r12eval_->ggspace(AnySpinCase1);
  // if density is already in required spaces?
  if (*orbs == *gspace)
    return rdm;

  Ref<LocalSCMatrixKit> local_kit = new LocalSCMatrixKit;
  RefSymmSCMatrix result = local_kit->symmmatrix(r12eval_->dim_gg(AnySpinCase2));
  result.assign(0.0);
  // it's possible for gspace to be a superset of orbs
  std::vector<int> map1 = map(*orbs, *gspace);
  SpinMOPairIter UV_iter(gspace->rank(),gspace->rank(),AlphaBeta);
  SpinMOPairIter PQ_iter(gspace->rank(),gspace->rank(),AlphaBeta);
  const int nmo = orbs->rank();

  for(PQ_iter.start(); int(PQ_iter); PQ_iter.next()) {
    const int P = PQ_iter.i();
    const int Q = PQ_iter.j();
    const int PQ = PQ_iter.ij();
    const int pp = map1[P];
    const int qq = map1[Q];
    if (pp == -1 || qq == -1) continue;   // skip if these indices are not in the source rdm2
    double pfac_pq = 1.0;
    int pq = pp * nmo + qq;

    for(UV_iter.start(); int(UV_iter); UV_iter.next()) {
      const int U = UV_iter.i();
      const int V = UV_iter.j();
      const int UV = UV_iter.ij();
      const int uu = map1[U];
      const int vv = map1[V];
      if (uu == -1 || vv == -1) continue;   // skip if these indices are not in the source rdm2
      double pfac_uv = 1.0;
      int uv = uu * nmo + vv; break;

      const double rdm_PQ_UV = pfac_pq * pfac_uv * rdm.get_element(pq, uv);
      result.set_element(PQ, UV, rdm_PQ_UV);
    }
  }

  return result;
}

int PT2R12::nelectron()
{
  return r12world()->refwfn()->nelectron();
}

RefSymmSCMatrix PT2R12::X_transformed_by_C() {
  RefSCMatrix T = C();
  RefSymmSCMatrix X = r12eval_->X();
  RefSCDimension transformed_dim = T.coldim();
  RefSymmSCMatrix XT = T.kit()->symmmatrix(transformed_dim);
  XT.assign(0.0);
  XT.accumulate_transform(T, X, SCMatrix::TransposeTransform);
  return(XT);
}

double
PT2R12::energy_recomputed_from_densities() {
//  ExEnv::out0() << std::endl<< indent << "Entered energy_recomputed_from_densities" << std::endl;
  double twoparticle_energy;
  double oneparticle_energy;

  {
    const RefSymmSCMatrix H = compute_obints<&Integral::hcore>(rdm1_->orbs());
    RefSymmSCMatrix opdm = rdm1();
    // H and opdm might use different SCMatrixKits -> copy H into another matrix with the same kit as opdm
    RefSymmSCMatrix hh = opdm.clone();
    hh->convert(H);
    hh->element_op(new SCElementScaleDiagonal(0.5));

    Ref<SCElementScalarProduct> trace_op = new SCElementScalarProduct;
    hh->element_op(trace_op, opdm);
    oneparticle_energy = trace_op->result() * 2.0;
  }

  {
    const Ref<OrbitalSpace>& space = rdm2_->orbs();

    const RefSymmSCMatrix tpdm = rdm2();
    RefSCMatrix G = g(space, space);
    // G and opdm might use different SCMatrixKits -> copy G into another matrix with the same kit as opdm
    RefSymmSCMatrix gg = tpdm.clone();  gg.assign(0.0);
    gg->accumulate_subblock(G.pointer(), 0, G.nrow()-1, 0, G.ncol() - 1);  // his automatically scales off diagonal elements by 2
                                                                           // no need to scale the diagonal when taking scalar product
    G = 0;

    Ref<SCElementScalarProduct> trace_op = new SCElementScalarProduct;
    gg->element_op(trace_op, tpdm);
    twoparticle_energy = trace_op->result();
    twoparticle_energy /= 2.0; // 1/2 comes from the formula
  }

#define ENERGY_CONVENTIONAL_PRINT_CONTRIBUTIONS 1

#if ENERGY_CONVENTIONAL_PRINT_CONTRIBUTIONS
  ExEnv::out0() << "1-e energy: "
                << setprecision(12) << oneparticle_energy << endl;
  ExEnv::out0() << "2-e energy: "
                << setprecision(12) << twoparticle_energy << endl;
#endif
  double energy = oneparticle_energy + twoparticle_energy;
  energy += this->r12world()->refwfn()->basis()->molecule()->nuclear_repulsion_energy();
//  ExEnv::out0() << std::endl<< indent << "Exited energy_recomputed_from_densities" << std::endl << std::endl << std::endl;

  return(energy);
}

RefSCMatrix PT2R12::g(const Ref<OrbitalSpace>& space1,
                      const Ref<OrbitalSpace>& space2) {
  const Ref<SCMatrixKit> localkit = new LocalSCMatrixKit;
  Ref<SpinMOPairIter> PQ_iter = new SpinMOPairIter(space1->rank(),space2->rank(),AlphaBeta);
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

  std::vector<string> tforms;
  {
    const string tform_key = ParsedTwoBodyFourCenterIntKey::key(s1->id(),s2->id(),
                                                                     s1->id(),s2->id(),
                                                                     string("ERI"),
                                                                     string(TwoBodyIntLayout::b1b2_k1k2));
    tforms.push_back(tform_key);
  }

  r12eval_->compute_tbint_tensor<ManyBodyTensors::I_to_T,false,false>(G,TwoBodyOper::eri,s1,s1,s2,s2,
      false,tforms);

  if (registered_space1) oreg->remove(space1->id());
  if (registered_space2) oreg->remove(space2->id());

  return(G);
}

RefSCMatrix PT2R12::g(const Ref<OrbitalSpace>& bra1,
                      const Ref<OrbitalSpace>& bra2,
                      const Ref<OrbitalSpace>& ket1,
                      const Ref<OrbitalSpace>& ket2)
{
      const Ref<SCMatrixKit> localkit = new LocalSCMatrixKit;
      Ref<SpinMOPairIter> braiter12 = new SpinMOPairIter(bra1->rank(),bra2->rank(),AlphaBeta);
      Ref<SpinMOPairIter> ketiter12 = new SpinMOPairIter(ket1->rank(),ket2->rank(),AlphaBeta);
      RefSCMatrix G = localkit->matrix(new SCDimension(braiter12->nij()),
                                       new SCDimension(ketiter12->nij()));
      G.assign(0.0);

      // find equivalent spaces in the registry
      Ref<OrbitalSpaceRegistry> oreg = this->r12world()->world()->tfactory()->orbital_registry();
      if (!oreg->value_exists(bra1) || !oreg->value_exists(bra2) || !oreg->value_exists(ket1) || !oreg->value_exists(ket2))
        throw ProgrammingError("PT2R12::g() -- spaces must be registered",__FILE__,__LINE__);


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

      for(braiter12->start(); *braiter12; braiter12->next())
      {
        const int b1 = braiter12->i();
        const int b2 = braiter12->j();
        const int b12 = braiter12->ij();

        const double* blk_b1b2 = da4_b1k1_b2k2->retrieve_pair_block(b1, b2, 0);

        for(ketiter12->start(); *ketiter12; ketiter12->next())
        {
          const int k1 = ketiter12->i();
          const int k2 = ketiter12->j();
          const int k12 = ketiter12->ij();

          G.set_element( b12, k12, blk_b1b2[k12] );
        }

        da4_b1k1_b2k2->release_pair_block(b1, b2, 0);
      }
  return(G);
}

RefSCMatrix PT2R12::f() {
  Ref<OrbitalSpace> space = rdm1_->orbs();
  Ref<OrbitalSpaceRegistry> oreg = r12world()->world()->tfactory()->orbital_registry();
  if (!oreg->value_exists(space)) {
    oreg->add(make_keyspace_pair(space));
  }
  const string key = oreg->key(space);
  space = oreg->value(key);
  RefSCMatrix fmat = r12eval_->fock(space,space,AnySpinCase1);

  return(fmat);
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:


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

#ifdef __GNUG__
#pragma implementation
#endif

#include <chemistry/qc/mbptr12/pt2r12.h>
#include <chemistry/qc/mbptr12/print.h>
#include <math/scmat/local.h>
#include <chemistry/qc/mbptr12/compute_tbint_tensor.h>

using namespace std;
using namespace sc;

namespace {

  int triang_half_INDEX_ordered(int i, int j) {
    return(i*(i+1)/2+j);
  }

  int triang_half_INDEX(int i, int j) {
    return((i>j) ? triang_half_INDEX_ordered(i,j) : triang_half_INDEX_ordered(j,i));
  }

  int ordinary_INDEX(int i, int j, int coldim) {
    return(i*coldim+j);
  }

  /// tpdm_index and init_ioff: for indexing of density matrices.
  int tpdm_index(int i, int j, int k, int l,int coldim){
    //int ind_half1=triang_half_INDEX(i,j);
    //int ind_half2=triang_half_INDEX(k,l);
    int ind_half1=ordinary_INDEX(i,j,coldim);
    int ind_half2=ordinary_INDEX(k,l,coldim);
    return(triang_half_INDEX(ind_half1,ind_half2));
  }

  void vector_to_symmmatrix(RefSymmSCMatrix &matrix, const RefSCVector &vector) {
    int dim = matrix.dim().n();
    for(int i=0; i<dim; i++){
      for(int j=0; j<=i; j++) {
        matrix.set_element(i,j,vector.get_element(triang_half_INDEX(i,j)));
      }
    }
  }

  void symmmatrix_to_vector(RefSCVector &vector, const RefSymmSCMatrix &matrix) {
    int dim = matrix.dim().n();
    for(int i=0; i<dim; i++){
      for(int j=0; j<=i; j++) {
        vector.set_element(triang_half_INDEX(i,j),matrix.get_element(i,j));
      }
    }
  }

  void vector_to_matrix(RefSCMatrix &matrix, const RefSCVector &vector) {
    int dim1 = matrix.rowdim().n();
    int dim2 = matrix.coldim().n();
    for(int i=0; i<dim1; i++) {
      for(int j=0; j<dim2; j++) {
        matrix.set_element(i,j,vector.get_element(ordinary_INDEX(i,j,dim2)));
      }
    }
  }

  int lowerupper_index(int p, int q);

  void vector_to_matrix(RefSCMatrix &matrix,const RefSCVector &vector,const SpinCase2 &pairspin) {
    int dim1 = matrix.rowdim().n();
    int dim2 = matrix.coldim().n();
    if(pairspin==AlphaBeta) {
      for(int i=0; i<dim1; i++) {
        for(int j=0; j<dim2; j++) {
          matrix.set_element(i,j,vector.get_element(ordinary_INDEX(i,j,dim2)));
        }
      }
    }
    else {  // pairspin==AlphaAlpha || pairspin==BetaBeta
      matrix->assign(0.0);
      for(int i=0; i<dim1; i++) {
        for(int j=0; j<i; j++) {
          const double value = vector.get_element(lowerupper_index(i,j));
          matrix.set_element(i,j,value);
          matrix.set_element(j,i,-value);
        }
      }
    }
  }

  void matrix_to_vector(RefSCVector &vector, const RefSCMatrix &matrix) {
    int dim1 = matrix.rowdim().n();
    int dim2 = matrix.coldim().n();
    for(int i=0; i<dim1; i++) {
      for(int j=0; j<dim2; j++) {
        vector.set_element(ordinary_INDEX(i,j,dim2),matrix.get_element(i,j));
      }
    }
  }

  void matrix_to_vector(RefSCVector &vector, const RefSCMatrix &matrix,const SpinCase2 &pairspin) {
    int dim1 = matrix.rowdim().n();
    int dim2 = matrix.coldim().n();
    if(pairspin==AlphaBeta) {
      for(int i=0; i<dim1; i++) {
        for(int j=0; j<dim2; j++) {
          vector.set_element(ordinary_INDEX(i,j,dim2),matrix.get_element(i,j));
        }
      }
    }
    else {  // pairspin==AlphaAlpha || pairspin==BetaBeta
      for(int i=0; i<dim1; i++) {
        for(int j=0; j<i; j++) {
          vector.set_element(lowerupper_index(i,j),matrix.get_element(i,j));
        }
      }
    }
  }

  int lowertriang_index(int p,int q) {
    if(q>=p){
      throw ProgrammingError("lowertriang_index(p,q) -- q must be smaller than p.",__FILE__,__LINE__);
    }
    int index=p*(p+1)/2+q-p;
    return(index);
  }

  int lowerupper_index(int p, int q) {
    if(p>q) {
      return(lowertriang_index(p,q));
    }
    else if(q>p) {
      return(lowertriang_index(q,p));
    }
    else {
      throw ProgrammingError("lowerupper_index(p,q) -- p and q are not allowed to be equal.",__FILE__,__LINE__);
    }
  }

  double indexsizeorder_sign(int p,int q) {
    if(p>q) {
      return(1.0);
    }
    else if(q>p) {
      return(-1.0);
    }
    else {
      return(0.0);
    }
  }

}

////////////////////

extern Ref<R12RefWavefunction>
make_PsiSD_R12RefWavefunction(const Ref<WavefunctionWorld>& world,
                              const Ref<Wavefunction>& wfn,
                              bool spin_restricted = true,
                              unsigned int nfzc = 0,
                              unsigned int nfzv = 0,
                              Ref<OrbitalSpace> vir_space = 0);

static ClassDesc PT2R12_cd(typeid(PT2R12),"PT2R12",
                           1,"public Wavefunction",0,create<PT2R12>,create<PT2R12>);

PT2R12::PT2R12(const Ref<KeyVal> &keyval) : Wavefunction(keyval)
{
  reference_ = require_dynamic_cast<Wavefunction*>(
        keyval->describedclassvalue("reference").pointer(),
        "PT2R12::PT2R12\n"
        );
  rdm2_ = require_dynamic_cast<RDM<Two>*>(
        keyval->describedclassvalue("rdm2").pointer(),
        "PT2R12::PT2R12\n"
        );
  assert(reference_ == rdm2_->wfn());

  Ref<WavefunctionWorld> world = new WavefunctionWorld(keyval, this);
  //world->memory(memory);
  const bool spin_restricted = true;
#if 0
  RefSymmSCMatrix rdm1_a = reference_->alpha_ao_density();
  RefSymmSCMatrix rdm1_b = reference_->spin_polarized() ? reference_->beta_ao_density() : rdm1_a;
  Ref<R12RefWavefunction> ref = new ORDM_R12RefWavefunction(world, basis(),
                                                            integral(),
                                                            rdm1_a, rdm1_b,
                                                            spin_restricted,
                                                            nfzc_,
                                                            omit_virtuals_);
#endif
  Ref<R12RefWavefunction> ref = make_PsiSD_R12RefWavefunction(world,
                                                              reference_,
                                                              spin_restricted,
                                                              nfzc_,
                                                              0);
  r12world_ = new R12WavefunctionWorld(keyval, ref);
  r12eval_ = new R12IntEval(r12world_);

  tpdm_from_opdms_ = keyval->booleanvalue("tpdm_from_opdms",KeyValValueboolean(false));
  debug_ = keyval->intvalue("debug", KeyValValueint(0));
  r12eval_->debug(debug_);
}

PT2R12::PT2R12(StateIn &s) : Wavefunction(s) {
  reference_ << SavableState::restore_state(s);
  r12world_ << SavableState::restore_state(s);
  r12eval_ << SavableState::restore_state(s);
  s.get(tpdm_from_opdms_);
  s.get(debug_);
}

PT2R12::~PT2R12() {}

void PT2R12::save_data_state(StateOut &s) {
  Wavefunction::save_data_state(s);
  SavableState::save_state(reference_, s);
  SavableState::save_state(r12world_, s);
  SavableState::save_state(r12eval_, s);
  s.put(tpdm_from_opdms_);
  s.put(debug_);
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

RefSymmSCMatrix PT2R12::hcore_mo(SpinCase1 spin) {
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

RefSymmSCMatrix PT2R12::moints() {
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

RefSCMatrix PT2R12::moints(SpinCase2 pairspin) {
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

RefSCMatrix PT2R12::g(SpinCase2 pairspin) {
  const Ref<SCMatrixKit> localkit = new LocalSCMatrixKit;
  const SpinCase1 spin1 = case1(pairspin);
  const SpinCase1 spin2 = case2(pairspin);
  const Ref<OrbitalSpace> space1 = r12eval_->orbs(spin1);
  const Ref<OrbitalSpace> space2 = r12eval_->orbs(spin2);
  Ref<SpinMOPairIter> PQ_iter = new SpinMOPairIter(space1,space2,pairspin);
  const int nmo1 = space1->rank();
  const int nmo2 = space2->rank();
  int nmo;
  if(nmo1==nmo2) {
    nmo = nmo1;
  }
  else {
    throw ProgrammingError("PT2R12::g(SpinCase2) -- nmo dimension of the spaces does not fit.",__FILE__,__LINE__);
  }
  const int pairrank = PQ_iter->nij();
  const RefSCDimension pairdim = new SCDimension(pairrank);
  RefSCMatrix G = localkit->matrix(pairdim,pairdim);
  G.assign(0.0);

  const bool antisymmetrize = (pairspin==AlphaBeta) ? false : true;
  std::vector<std::string> tforms;
  {
    const std::string tform_key = ParsedTwoBodyFourCenterIntKey::key(space1->id(),space2->id(),
                                                                     space1->id(),space2->id(),
                                                                     std::string("ERI"),
                                                                     std::string(TwoBodyIntLayout::b1b2_k1k2));
    tforms.push_back(tform_key);
  }

  r12eval_->compute_tbint_tensor<ManyBodyTensors::I_to_T,false,false>(G,TwoBodyOper::eri,space1,space1,space2,space2,
      antisymmetrize,tforms);

  return(G);
}

RefSCMatrix PT2R12::f(SpinCase1 spin) {
  const Ref<OrbitalSpace> bra_space = r12eval_->orbs(spin);
  const Ref<OrbitalSpace> ket_space = r12eval_->orbs(spin);
  RefSCMatrix fmat = r12eval_->fock(bra_space,ket_space,spin);

  return(fmat);
}


RefSCMatrix PT2R12::C(SpinCase2 S) {
  //return(r12energy_->C(S));
  Ref<LocalSCMatrixKit> local_matrix_kit = new LocalSCMatrixKit();
  RefSCMatrix Cmat = local_matrix_kit->matrix(r12eval_->dim_GG(S),r12eval_->dim_gg(S));
  if(S==AlphaBeta) {
    SpinMOPairIter OW_iter(r12eval_->GGspace(Alpha), r12eval_->GGspace(Beta), S );
    SpinMOPairIter PQ_iter(r12eval_->ggspace(Alpha), r12eval_->GGspace(Beta), S );
    Ref<LinearR12::GeminalDescriptor> geminaldesc = r12world()->r12tech()->corrfactor()->geminaldescriptor();
    CuspConsistentGeminalCoefficient coeff_gen(S,geminaldesc);
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
    SpinMOPairIter OW_iter(r12eval_->GGspace(spin), r12eval_->GGspace(spin), S );
    SpinMOPairIter PQ_iter(r12eval_->ggspace(spin), r12eval_->GGspace(spin), S );
    Ref<LinearR12::GeminalDescriptor> geminaldesc = r12world()->r12tech()->corrfactor()->geminaldescriptor();
    CuspConsistentGeminalCoefficient coeff_gen(S,geminaldesc);
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

RefSCMatrix PT2R12::V_genref_projector2(SpinCase2 pairspin) {
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

#if 0
  V_intermed.print(prepend_spincase(pairspin, "V").c_str());
  tpdm.print(prepend_spincase(pairspin, "gamma").c_str());
  V_genref.print(prepend_spincase(pairspin, "V.gamma").c_str());
#endif

  return(V_genref);
}

RefSCMatrix PT2R12::V_transformed_by_C(SpinCase2 pairspin) {
  RefSCMatrix T = C(pairspin);
  RefSCMatrix V = r12eval_->V(pairspin);
  RefSCMatrix Vtransposed = V.t();
  RefSCMatrix VT = Vtransposed*T;
  return(VT);
}

RefSymmSCMatrix PT2R12::X_transformed_by_C(SpinCase2 pairspin) {
  RefSCMatrix T = C(pairspin);
  RefSymmSCMatrix X = r12eval_->X(pairspin);
  RefSCDimension transformed_dim = T.coldim();
  RefSymmSCMatrix XT = T.kit()->symmmatrix(transformed_dim);
  XT.assign(0.0);
  XT.accumulate_transform(T,X,SCMatrix::TransposeTransform);

  return(XT);
}

RefSymmSCMatrix PT2R12::B_transformed_by_C(SpinCase2 pairspin) {
  RefSymmSCMatrix B = r12eval_->B(pairspin);
  RefSCMatrix T = C(pairspin);
  RefSCDimension transformed_dim = T.coldim();
  RefSymmSCMatrix BT = T.kit()->symmmatrix(transformed_dim);
  BT.assign(0.0);
  BT.accumulate_transform(T,B,SCMatrix::TransposeTransform);

  return(BT);
}

#if 0
namespace {

RefSCVector dipolemoments(const Ref<DipoleData> &dipoledata) {
  Ref<SCMatrixKit> localkit = new LocalSCMatrixKit;
  Ref<GaussianBasisSet> basis = this->basis();
  Ref<OrbitalSpace> space = r12eval_->orbs(Alpha);
  const int nmo = reference_->nmo();

  // computing dipole integrals in MO basis
  RefSCMatrix MX, MY, MZ;
  RefSCMatrix MXX, MYY, MZZ, MXY, MXZ, MYZ;
  compute_multipole_ints(space,space,MX,MY,MZ,MXX,MYY,MZZ,MXY,MXZ,MYZ);

  // computing oneparticle density matrix in MO basis
  RefSymmSCMatrix opdm = onepdm_transformed();
  //opdm.print("oneparticle density matrix in MO basis");

  RefSCDimension M_dim = opdm.dim();
  RefSCMatrix MX_nb = localkit->matrix(M_dim,M_dim);
  RefSCMatrix MY_nb = localkit->matrix(M_dim,M_dim);
  RefSCMatrix MZ_nb = localkit->matrix(M_dim,M_dim);
  for(int i=0; i<M_dim.n(); i++) {
    for(int j=0; j<M_dim.n(); j++) {
      MX_nb.set_element(i,j,MX.get_element(i,j));
      MY_nb.set_element(i,j,MY.get_element(i,j));
      MZ_nb.set_element(i,j,MZ.get_element(i,j));
    }
  }

  RefSCMatrix opdm_x = MX_nb * opdm;
  RefSCMatrix opdm_y = MY_nb * opdm;
  RefSCMatrix opdm_z = MZ_nb * opdm;

  RefSCVector dipoles = localkit->vector(RefSCDimension(new SCDimension(3)));
  dipoles.assign(0.0);

  double dx=0.0, dy=0.0, dz=0.0;
  for(int p=0; p<nmo; p++) {
    dx += opdm_x.get_element(p,p);
    dy += opdm_y.get_element(p,p);
    dz += opdm_z.get_element(p,p);
  }

  int natom = reference_->molecule()->natom();
  double dx_nuc=0.0, dy_nuc=0.0, dz_nuc=0.0;
  for(int i=0; i<natom; i++) {
    dx_nuc += reference_->molecule()->r(i,0) * reference_->molecule()->Z(i);
    dy_nuc += reference_->molecule()->r(i,1) * reference_->molecule()->Z(i);
    dz_nuc += reference_->molecule()->r(i,2) * reference_->molecule()->Z(i);
  }

  dipoles[0] = dx + dx_nuc;
  dipoles[1] = dy + dy_nuc;
  dipoles[2] = dz + dz_nuc;

  return(dipoles);
}


double energy_conventional() {
  double twoparticle_energy[NSpinCases2];
  double oneparticle_energy[NSpinCases1];
  const int npurespincases2 = (r12eval_->spin_polarized()) ? 3 : 2;
  const int npurespincases1 = (r12eval_->spin_polarized()) ? 2 : 1;
  const int nmo = r12world()->ref()->orbs()->rank();
  const RefSCDimension nmodim = new SCDimension(nmo);
  const Ref<SCMatrixKit> localkit = new LocalSCMatrixKit;

  for(int spincase1=0; spincase1<npurespincases1; spincase1++) {
    const SpinCase1 spin = static_cast<SpinCase1>(spincase1);
    const RefSymmSCMatrix opdm = onepdm_refmo(spin);
    const RefSymmSCMatrix hcore = hcore_mo(spin);
    oneparticle_energy[spin] = (hcore*opdm)->trace();
  }

  for(int spincase2=0; spincase2<npurespincases2; spincase2++) {
    const SpinCase2 pairspin = static_cast<SpinCase2>(spincase2);
    const SpinCase1 spin1 = case1(pairspin);
    const SpinCase1 spin2 = case2(pairspin);

    const RefSymmSCMatrix tpdm = twopdm_refmo(pairspin);
#if 0
    if(pairspin==AlphaBeta) {
      ExEnv::out0() << "start test in PT2R12::energy_conventional() ..." << endl;
      const double tol = 1.0e-6;
      const RefSymmSCMatrix tpdm_vgl = twopdm_transformed_dirac();
      for(int p=0; p<nmo; p++) {
        for(int q=0; q<nmo; q++) {
          const int ind_pq = ordinary_INDEX(p,q,nmo);
          for(int r=0; r<nmo; r++) {
            for(int s=0; s<nmo; s++) {
              const int ind_rs = ordinary_INDEX(r,s,nmo);
              const double val1 = 2.0 * tpdm->get_element(ind_pq,ind_rs);
              const double val2 = tpdm_vgl->get_element(ind_pq,ind_rs);
              if(fabs(val1-val2)>tol) {
                ExEnv::out0() << p << "  " << q << "  " << r << "  " << s
                              << "  " << setprecision(12) << val1
                              << "  " << setprecision(12) << val2 << endl;
              }
            }
          }
        }
      }
      ExEnv::out0() << "Test in PT2R12::energy_conventional() finished." << endl;
    }
#endif
    const RefSCMatrix G = g(pairspin);

    twoparticle_energy[pairspin] = (G*tpdm).trace();
  }

  if(!r12eval_->spin_polarized()) {
    twoparticle_energy[BetaBeta] = twoparticle_energy[AlphaAlpha];
    oneparticle_energy[Beta] = oneparticle_energy[Alpha];
  }

//#define ENERGY_CONVENTIONAL_PRINT_CONTRIBUTIONS

  double energy_hcore = 0.0;
  for(int i=0; i<NSpinCases1; i++) {
    const SpinCase1 spin = static_cast<SpinCase1>(i);
#ifdef ENERGY_CONVENTIONAL_PRINT_CONTRIBUTIONS
    ExEnv::out0() << ((spin==Alpha) ? "Alpha" : "Beta") << " hcore energy: "
                  << setprecision(12) << oneparticle_energy[i] << endl;
#endif
    energy_hcore += oneparticle_energy[i];
  }
#ifdef ENERGY_CONVENTIONAL_PRINT_CONTRIBUTIONS
  ExEnv::out0() << "hcore energy: " << setprecision(12) << energy_hcore << endl;
#endif
  double energy_twoelec = 0.0;
  for(int i=0; i<NSpinCases2; i++) {
#ifdef ENERGY_CONVENTIONAL_PRINT_CONTRIBUTIONS
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
#ifdef ENERGY_CONVENTIONAL_PRINT_CONTRIBUTIONS
  ExEnv::out0() << "two-electron energy: " << setprecision(12) << energy_twoelec << endl;
#endif
  double energy = energy_hcore + energy_twoelec;
  energy += reference_->nuclear_repulsion_energy();

  return(energy);
}

RefSCMatrix phi_HF(SpinCase2 pairspin) {
  RefSCMatrix fmat[NSpinCases1];
  for(int i=0; i<NSpinCases1; i++) {
    SpinCase1 spin = static_cast<SpinCase1>(i);
    fmat[i] = f(spin);

    //fmat[i].print(prepend_spincase(spin,"Fock matrix").c_str());
  }
  int nmo = reference()->nmo();
  RefSCMatrix phi;
  phi = C(pairspin).kit()->matrix(r12eval_->dim_gg(pairspin),r12eval_->dim_gg(pairspin));
  phi.assign(0.0);

  if(pairspin==AlphaBeta) {
    SpinMOPairIter UV_iter(r12eval_->orbs(Alpha),r12eval_->orbs(Beta),pairspin);
    SpinMOPairIter PQ_iter(r12eval_->orbs(Alpha),r12eval_->orbs(Beta),pairspin);
    for(PQ_iter.start(); int(PQ_iter); PQ_iter.next()) {
      const int P = PQ_iter.i();
      const int Q = PQ_iter.j();
      const int PQ = PQ_iter.ij();
      for(UV_iter.start(); int(UV_iter); UV_iter.next()) {
        const int U = UV_iter.i();
        const int V = UV_iter.j();
        const int UV = UV_iter.ij();
        if(U==P) {
          phi->accumulate_element(PQ,UV,-fmat[Beta].get_element(Q,V));
        }
        if(Q==V) {
          phi->accumulate_element(PQ,UV,-fmat[Alpha].get_element(P,U));
        }
      }
    }
  }  // pairspin==AlphaBeta
  else {  // pairspin!=AlphaBeta
    SpinCase1 spin = (pairspin == AlphaAlpha) ? Alpha : Beta;
    SpinCase1 otherspin = (spin==Alpha) ? Beta : Alpha;
    SpinMOPairIter UV_iter(r12eval_->orbs(spin),r12eval_->orbs(spin),pairspin);
    SpinMOPairIter PQ_iter(r12eval_->orbs(spin),r12eval_->orbs(spin),pairspin);

    for(PQ_iter.start(); int(PQ_iter); PQ_iter.next()) {
      const int P = PQ_iter.i();
      const int Q = PQ_iter.j();
      const int PQ = PQ_iter.ij();
      for(UV_iter.start(); int(UV_iter); UV_iter.next()) {
        const int U = UV_iter.i();
        const int V = UV_iter.j();
        const int UV = UV_iter.ij();

        if(V==Q) {
          phi->accumulate_element(PQ,UV,-1.0*fmat[spin].get_element(P,U));
        }
        if(U==Q) {
          phi->accumulate_element(PQ,UV,fmat[spin].get_element(P,V));
        }
        if(V==P) {
          phi->accumulate_element(PQ,UV,fmat[spin].get_element(Q,U));
        }
        if(U==P) {
          phi->accumulate_element(PQ,UV,-1.0*fmat[spin].get_element(Q,V));
        }

      }
    }
  }  // pairspin!=AlphaBeta
  return(phi);
}

RefSCMatrix phi_twoelec(SpinCase2 pairspin) {
  double energ_PT2R12[NSpinCases2];
  RefSCMatrix fmat[NSpinCases1];
  RefSymmSCMatrix opdm[NSpinCases1];
  const RefSymmSCMatrix tpdm = twopdm_refmo(pairspin);

  for(int i=0; i<NSpinCases1; i++) {
    SpinCase1 spin = static_cast<SpinCase1>(i);
    fmat[i] = f(spin);
    opdm[i] = onepdm_refmo(spin);
  }
  const int nmo = opdm[0].dim().n();

  double J = 0.0;
  for(int i=0; i<NSpinCases1; i++) {
    for(int z=0; z<nmo; z++) {
      for(int y=0; y<nmo; y++) {
        J += fmat[i].get_element(z,y)*opdm[i].get_element(y,z);
      }
    }
  }

  RefSCMatrix phi;
  phi = C(pairspin).kit()->matrix(r12eval_->dim_gg(pairspin),r12eval_->dim_gg(pairspin));
  phi.assign(0.0);
  phi.accumulate(tpdm*J);

  // redefined the sign for phi to correspond directly to the Fock operator
  phi.scale(-1.0);

  if(debug_>=DefaultPrintThresholds::mostO4) {
    phi->print(prepend_spincase(pairspin,"phi").c_str());
  }

  return(phi);
}


//#define DEBUG_PHI
// In this function, the statements commented out with "#if 0" would be more
// efficient, but are not tested yet.
RefSCMatrix PT2R12::phi(SpinCase2 pairspin) {
  RefSCMatrix fmat[NSpinCases1];
  RefSymmSCMatrix opdm[NSpinCases1];
  RefSymmSCMatrix tpdm[NSpinCases2];

  for(int i=0; i<NSpinCases1; i++) {
    SpinCase1 spin = static_cast<SpinCase1>(i);
    fmat[i] = f(spin);
    opdm[i] = rdm1(spin);
  }
  const int nmo = opdm[0].dim().n();

  // R12 intermediates, transformed with geminal amplitude matrix
  for(int i=0; i<NSpinCases2; i++) {
    SpinCase2 S = static_cast<SpinCase2>(i);
    tpdm[i] = rdm2(S);
  }

  // intermediates
  RefSCMatrix I[NSpinCases1];
  RefSCMatrix K[NSpinCases1];
  double J = 0.0;
  RefSCMatrix L[NSpinCases1];
  RefSCMatrix M[NSpinCases1];
  for(int i=0; i<NSpinCases1; i++) {
#if 1
    I[i] = fmat[i].clone();
    I[i].assign(0.0);
    for(int q=0; q<nmo; q++) {
      for(int v=0; v<nmo; v++) {
        for(int z=0; z<nmo; z++) {
          for(int y=0; y<nmo; y++) {
            I[i].accumulate_element(q,v,opdm[i].get_element(q,z)*fmat[i].get_element(z,y)*opdm[i].get_element(y,v));
          }
        }
      }
    }
#endif
#if 0
    I[i] = opdm[i] * fmat[i] * opdm[i];
#endif

#if 1
    for(int z=0; z<nmo; z++) {
      for(int y=0; y<nmo; y++) {
        J += fmat[i].get_element(z,y)*opdm[i].get_element(y,z);
      }
    }
#endif
#if 0
    J += (fmat[i] * opdm[i]).trace();
#endif

#if 1
    K[i] = fmat[i].clone();
    K[i].assign(0.0);
    for(int z=0; z<nmo; z++) {
      for(int u=0; u<nmo; u++) {
        for(int y=0; y<nmo; y++) {
          K[i].accumulate_element(z,u,opdm[i].get_element(y,u)*fmat[i].get_element(z,y));
        }
      }
    }
#endif
#if 0
    K[i] = opdm[i] * fmat[i].t();
#endif

#if 1
    L[i] = fmat[i].clone();
    L[i].assign(0.0);
    for(int p=0; p<nmo; p++) {
      for(int y=0; y<nmo; y++) {
        for(int z=0; z<nmo; z++) {
          L[i].accumulate_element(p,y,opdm[i].get_element(p,z)*fmat[i].get_element(z,y));
        }
      }
    }
#endif
#if 0
    L[i] = opdm[i] *  fmat[i];
#endif

    M[i] = fmat[i].clone();
    M[i].assign(0.0);
    if(i==Alpha) {
      for(int q=0; q<nmo; q++) {
        for(int y=0; y<nmo; y++) {
          for(int v=0; v<nmo; v++) {
            for(int z=0; z<nmo; z++) {
              int qy_alphabeta = q*nmo + y;
              int vz_alphabeta = v*nmo + z;
              M[i].accumulate_element(q,v,fmat[Beta].get_element(z,y)*tpdm[AlphaBeta].get_element(qy_alphabeta,vz_alphabeta));
              if((q!=y) && (v!=z)) {
                int qy_alphaalpha = lowerupper_index(q,y);
                int vz_alphaalpha = lowerupper_index(v,z);
                int sign = indexsizeorder_sign(q,y)*indexsizeorder_sign(v,z);
                M[i].accumulate_element(q,v,sign*fmat[Alpha].get_element(z,y)*tpdm[AlphaAlpha].get_element(qy_alphaalpha,vz_alphaalpha));
              }
            }
          }
        }
      }

    }  // i==Alpha
    else {  // i!=Alpha
      for(int q=0; q<nmo; q++) {
        for(int v=0; v<nmo; v++) {
          for(int z=0; z<nmo; z++) {
            for(int y=0; y<nmo; y++) {
              int yq_alphabeta = y*nmo + q;
              int zv_alphabeta = z*nmo + v;
              M[i].accumulate_element(q,v,fmat[Alpha].get_element(z,y)*tpdm[AlphaBeta].get_element(yq_alphabeta,zv_alphabeta));
              if((q!=y) && (v!=z)) {
                int qy_betabeta = lowerupper_index(q,y);
                int vz_betabeta = lowerupper_index(v,z);
                int sign = indexsizeorder_sign(q,y)*indexsizeorder_sign(v,z);
                M[i].accumulate_element(q,v,sign*fmat[Beta].get_element(z,y)*tpdm[BetaBeta].get_element(qy_betabeta,vz_betabeta));
              }
            }
          }
        }
      }
    }  // i!=Alpha

    const SpinCase1 spin = static_cast<SpinCase1>(i);
    L[i].print(prepend_spincase(spin,"K").c_str());
    I[i].print(prepend_spincase(spin,"I").c_str());
    M[i].print(prepend_spincase(spin,"M").c_str());

  }

  // computing phi
  RefSCMatrix phi;
  phi = C(pairspin).kit()->matrix(r12eval_->dim_gg(pairspin),r12eval_->dim_gg(pairspin));
  phi.assign(0.0);
#ifdef DEBUG_PHI
  RefSCMatrix term1 = phi.clone();
  RefSCMatrix term2 = phi.clone();
  RefSCMatrix term3 = phi.clone();
  RefSCMatrix term4 = phi.clone();
  RefSCMatrix term5 = phi.clone();
  term1->assign(0.0);
  term2->assign(0.0);
  term3->assign(0.0);
  term4->assign(0.0);
  term5->assign(0.0);
#endif
  if(pairspin==AlphaBeta) {
    SpinMOPairIter UV_iter(r12eval_->orbs(Alpha),r12eval_->orbs(Beta),pairspin);
    SpinMOPairIter PQ_iter(r12eval_->orbs(Alpha),r12eval_->orbs(Beta),pairspin);
    for(PQ_iter.start(); int(PQ_iter); PQ_iter.next()) {
      int P = PQ_iter.i();
      int Q = PQ_iter.j();
      int PQ = PQ_iter.ij();
      for(UV_iter.start(); int(UV_iter); UV_iter.next()) {
        int U = UV_iter.i();
        int V = UV_iter.j();
        int UV = UV_iter.ij();

        phi.accumulate_element(PQ,UV,2.0*opdm[Alpha].get_element(P,U)*I[Beta].get_element(Q,V)
                                    +2.0*opdm[Beta].get_element(Q,V)*I[Alpha].get_element(P,U));

        phi.accumulate_element(PQ,UV,-2.0*opdm[Alpha].get_element(P,U)*opdm[Beta].get_element(Q,V)*J);

        for(int Z=0; Z<PQ_iter.ni(); Z++) {
          int ZV = Z*nmo + V;
          int UZ = U*nmo + Z;
          phi.accumulate_element(PQ,UV,-K[Alpha].get_element(Z,U)*tpdm[AlphaBeta].get_element(PQ,ZV)
                                       -K[Beta].get_element(Z,V)*tpdm[AlphaBeta].get_element(PQ,UZ));
        }

        for(int Y=0; Y<PQ_iter.ni(); Y++) {
          int YQ = Y*nmo + Q;
          int PY = P*nmo + Y;
          phi.accumulate_element(PQ,UV,-L[Alpha].get_element(P,Y)*tpdm[AlphaBeta].get_element(YQ,UV)
                                       -L[Beta].get_element(Q,Y)*tpdm[AlphaBeta].get_element(PY,UV));
        }

        phi.accumulate_element(PQ,UV,opdm[Alpha].get_element(P,U)*M[Beta].get_element(Q,V)
                                    +opdm[Beta].get_element(Q,V)*M[Alpha].get_element(P,U));

      }
    }
  }  // pairspin==AlphaBeta
  else {  // pairspin!=AlphaBeta
    SpinCase1 spin;
    SpinCase1 otherspin;
    if(pairspin==AlphaAlpha) {
      spin = Alpha;
      otherspin = Beta;
    }
    else {
      spin = Beta;
      otherspin = Alpha;
    }
    SpinMOPairIter UV_iter(r12eval_->orbs(spin),r12eval_->orbs(spin),pairspin);
    SpinMOPairIter PQ_iter(r12eval_->orbs(spin),r12eval_->orbs(spin),pairspin);
    for(PQ_iter.start(); int(PQ_iter); PQ_iter.next()) {
      int P = PQ_iter.i();
      int Q = PQ_iter.j();
      int PQ = PQ_iter.ij();
      double sign_PQ = indexsizeorder_sign(P,Q);
      int QP = lowerupper_index(Q,P);
      double sign_QP = indexsizeorder_sign(Q,P);

      for(UV_iter.start(); int(UV_iter); UV_iter.next()) {
        int U = UV_iter.i();
        int V = UV_iter.j();
        int UV = UV_iter.ij();

        if((U!=V) && (P!=Q)) {
          phi.accumulate_element(PQ,UV,2.0*opdm[spin].get_element(P,U)*I[spin].get_element(Q,V)
                                      -2.0*opdm[spin].get_element(Q,U)*I[spin].get_element(P,V)
                                      -2.0*opdm[spin].get_element(P,V)*I[spin].get_element(Q,U)
                                      +2.0*opdm[spin].get_element(Q,V)*I[spin].get_element(P,U));
#ifdef DEBUG_PHI
          term1.accumulate_element(PQ,UV,2.0*opdm[spin].get_element(P,U)*I[spin].get_element(Q,V)
                                      -2.0*opdm[spin].get_element(Q,U)*I[spin].get_element(P,V)
                                      -2.0*opdm[spin].get_element(P,V)*I[spin].get_element(Q,U)
                                      +2.0*opdm[spin].get_element(Q,V)*I[spin].get_element(P,U));
#endif

          phi.accumulate_element(PQ,UV,opdm[spin].get_element(P,V)*opdm[spin].get_element(Q,U)*J
                                      -opdm[spin].get_element(Q,V)*opdm[spin].get_element(P,U)*J
                                      -opdm[spin].get_element(P,U)*opdm[spin].get_element(Q,V)*J
                                      +opdm[spin].get_element(Q,U)*opdm[spin].get_element(P,V)*J);
#ifdef DEBUG_PHI
          term2.accumulate_element(PQ,UV,opdm[spin].get_element(P,V)*opdm[spin].get_element(Q,U)*J
                                        -opdm[spin].get_element(Q,V)*opdm[spin].get_element(P,U)*J
                                        -opdm[spin].get_element(P,U)*opdm[spin].get_element(Q,V)*J
                                        +opdm[spin].get_element(Q,U)*opdm[spin].get_element(P,V)*J);
#endif

          for(int Z=0; Z<UV_iter.ni(); Z++) {
            if(V!=Z) {
              int VZ = lowerupper_index(V,Z);
              double sign_VZ = indexsizeorder_sign(V,Z);
              phi.accumulate_element(PQ,UV,sign_VZ*K[spin].get_element(Z,U)*tpdm[pairspin].get_element(PQ,VZ));
#ifdef DEBUG_PHI
              term3.accumulate_element(PQ,UV,sign_VZ*K[spin].get_element(Z,U)*tpdm[pairspin].get_element(PQ,VZ));
#endif
            }
            if(U!=Z) {
              int UZ = lowerupper_index(U,Z);
              double sign_UZ = indexsizeorder_sign(U,Z);
              phi.accumulate_element(PQ,UV,-sign_UZ*K[spin].get_element(Z,V)*tpdm[pairspin].get_element(PQ,UZ));
#ifdef DEBUG_PHI
              term3.accumulate_element(PQ,UV,-sign_UZ*K[spin].get_element(Z,V)*tpdm[pairspin].get_element(PQ,UZ));
#endif
            }
          }

          for(int Y=0; Y<PQ_iter.ni(); Y++) {
            double sign_UV = indexsizeorder_sign(U,V);
            int VU = lowerupper_index(V,U);
            double sign_VU = indexsizeorder_sign(V,U);
            if(Q!=Y) {
              int QY = lowerupper_index(Q,Y);
              double sign_QY = indexsizeorder_sign(Q,Y);
              phi.accumulate_element(PQ,UV,sign_QY*L[spin].get_element(P,Y)*tpdm[pairspin].get_element(QY,UV));
#ifdef DEBUG_PHI
              term4.accumulate_element(PQ,UV,sign_QY*L[spin].get_element(P,Y)*tpdm[pairspin].get_element(QY,UV));
#endif
            }
            if(P!=Y) {
              int PY = lowerupper_index(P,Y);
              double sign_PY = indexsizeorder_sign(P,Y);
              phi.accumulate_element(PQ,UV,-sign_PY*L[spin].get_element(Q,Y)*tpdm[pairspin].get_element(PY,UV));
#ifdef DEBUG_PHI
              term4.accumulate_element(PQ,UV,-sign_PY*L[spin].get_element(Q,Y)*tpdm[pairspin].get_element(PY,UV));
#endif
            }
          }

          phi.accumulate_element(PQ,UV,opdm[spin].get_element(P,U)*M[spin].get_element(Q,V)
                                      -opdm[spin].get_element(Q,U)*M[spin].get_element(P,V)
                                      -opdm[spin].get_element(P,V)*M[spin].get_element(Q,U)
                                      +opdm[spin].get_element(Q,V)*M[spin].get_element(P,U));
#ifdef DEBUG_PHI
          term5.accumulate_element(PQ,UV,opdm[spin].get_element(P,U)*M[spin].get_element(Q,V)
                                      -opdm[spin].get_element(Q,U)*M[spin].get_element(P,V)
                                      -opdm[spin].get_element(P,V)*M[spin].get_element(Q,U)
                                      +opdm[spin].get_element(Q,V)*M[spin].get_element(P,U));
#endif

        }
      }
    }
#ifdef DEBUG_PHI
    term1->print(prepend_spincase(pairspin,"term1 of phi").c_str());
    term2->print(prepend_spincase(pairspin,"term2 of phi").c_str());
    term3->print(prepend_spincase(pairspin,"term3 of phi").c_str());
    term4->print(prepend_spincase(pairspin,"term4 of phi").c_str());
    term5->print(prepend_spincase(pairspin,"term5 of phi").c_str());
    RefSCMatrix term_sum = term1 + term2 + term3 + term4 + term5;
    term_sum->print(prepend_spincase(pairspin,"term_sum of phi").c_str());
#endif
  }  // pairspin!=AlphaBeta

  // redefined sign so that phi corresponds to the Fock matrix
  phi.scale(-1.0);

  if(debug_>=DefaultPrintThresholds::mostO4) {
    phi->print(prepend_spincase(pairspin,"phi").c_str());
  }

  return(phi);
}

} // end of anonymous namespace
#endif

RefSymmSCMatrix PT2R12::phi_cumulant(SpinCase2 spin12) {
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
    }
  }

  // compute phi:
  // \phi^{u v}_{p q} = & P(p q) P(u v) (\gamma^u_p \gamma^{q_3}_q f^{q_2}_{q_3} \gamma^v_{q_2}
  //                                     + \frac{1}{2} \gamma^v_{q_2} f^{q_2}_{q_3} \lambda^{u q_3}_{p q}
  //                                     + \frac{1}{2} \gamma^{q_2}_p f^{q_3}_{q_2} \lambda^{u v}_{q_3 q}
  //                                     - \gamma^u_p f^{q_2}_{q_3} \lambda^{q_3 v}_{q_2 q})
  Ref<LocalSCMatrixKit> local_kit = new LocalSCMatrixKit();
  RefSymmSCMatrix phi = local_kit->symmmatrix(r12eval_->dim_gg(spin12));
  phi.assign(0.0);
  const SpinCase1 spin1 = case1(spin12);
  const SpinCase1 spin2 = case2(spin12);
  Ref<OrbitalSpace> orbs1 = r12eval_->orbs(spin1);
  Ref<OrbitalSpace> orbs2 = r12eval_->orbs(spin2);
  Ref<OrbitalSpace> gspace1 = r12eval_->ggspace(spin1);
  Ref<OrbitalSpace> gspace2 = r12eval_->ggspace(spin2);
  MOIndexMap map1 = (*orbs1 << *gspace1);
  MOIndexMap map2 = (*orbs2 << *gspace2);
  SpinMOPairIter UV_iter(gspace1,gspace2,spin12);
  SpinMOPairIter PQ_iter(gspace1,gspace2,spin12);

  if (spin12 == AlphaBeta) {
    for(PQ_iter.start(); int(PQ_iter); PQ_iter.next()) {
      const int P = PQ_iter.i();
      const int Q = PQ_iter.j();
      const int PQ = PQ_iter.ij();
      const int pp = map1[P];
      const int qq = map2[Q];
      const int pq = pp * nmo + qq;

      for(UV_iter.start(); int(UV_iter); UV_iter.next()) {
        const int U = UV_iter.i();
        const int V = UV_iter.j();
        const int UV = UV_iter.ij();
        const int uu = map1[U];
        const int vv = map2[V];
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
      const int pp = map1[P];
      const int qq = map2[Q];
      assert(pp >= qq);
      const int pq = pp * (pp - 1) / 2 + qq;

      for(UV_iter.start(); int(UV_iter); UV_iter.next()) {
        const int U = UV_iter.i();
        const int V = UV_iter.j();
        const int UV = UV_iter.ij();
        const int uu = map1[U];
        const int vv = map2[V];
        assert(uu >= vv);
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
    assert(false);
  }

  if(debug_>=DefaultPrintThresholds::mostO4) {
    phi->print(prepend_spincase(spin12,"phi (new)").c_str());
  }

  return phi;
}

double PT2R12::energy_PT2R12_projector1(SpinCase2 pairspin) {
  const int nelectron = reference_->nelectron();
  SpinCase1 spin1 = case1(pairspin);
  SpinCase1 spin2 = case2(pairspin);
  Ref<OrbitalSpace> gg1space = r12eval_->ggspace(spin1);
  Ref<OrbitalSpace> gg2space = r12eval_->ggspace(spin2);
  SpinMOPairIter gg_iter(gg1space,gg2space,pairspin);

  RefSymmSCMatrix tpdm = rdm2_gg(pairspin);
  RefSCMatrix Phi = phi(pairspin);
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

double PT2R12::energy_PT2R12_projector2(SpinCase2 pairspin) {
  const int nelectron = reference_->nelectron();
  SpinCase1 spin1 = case1(pairspin);
  SpinCase1 spin2 = case2(pairspin);
  Ref<OrbitalSpace> gg1space = r12eval_->ggspace(spin1);
  Ref<OrbitalSpace> gg2space = r12eval_->ggspace(spin2);
  SpinMOPairIter gg_iter(gg1space,gg2space,pairspin);

  RefSymmSCMatrix TBT = B_transformed_by_C(pairspin);
  RefSymmSCMatrix tpdm = rdm2_gg(pairspin);
  RefSCMatrix HylleraasMatrix = TBT*tpdm;
  TBT=0;
  tpdm=0;

  RefSymmSCMatrix TXT = X_transformed_by_C(pairspin);
  RefSymmSCMatrix Phi = phi_cumulant(pairspin);
  RefSCMatrix TXT_t_Phi = TXT*Phi; TXT_t_Phi.scale(-1.0);
  HylleraasMatrix.accumulate(TXT_t_Phi);
  TXT=0;
  Phi=0;
  TXT_t_Phi=0;

  RefSCMatrix V_genref = V_genref_projector2(pairspin);
  RefSCMatrix T = C(pairspin);
  HylleraasMatrix.accumulate(2.0*V_genref.t()*T);
  T=0;
  V_genref=0;

  const double energy = this->compute_energy(HylleraasMatrix, pairspin);
  return(energy);
}

namespace {
  RefSymmSCMatrix convert_to_local_kit(const RefSymmSCMatrix& A) {
    RefSymmSCMatrix result;
    Ref<LocalSCMatrixKit> kit_cast_to_local; kit_cast_to_local << A.kit();
    if (kit_cast_to_local.null()) {
      Ref<LocalSCMatrixKit> local_kit = new LocalSCMatrixKit();
      RefSymmSCMatrix A_local = local_kit->symmmatrix(A.dim());
      A_local->convert(A);
      result = A_local;
    }
    else
      result = A;
    return result;
  }
}

RefSymmSCMatrix sc::PT2R12::rdm1(SpinCase1 spin)
{
  return r12world_->ref()->ordm_orbs_sb(spin);
}

RefSymmSCMatrix sc::PT2R12::rdm2(SpinCase2 spin)
{
  // since LocalSCMatrixKit is used everywhere, convert to Local kit
  return convert_to_local_kit(rdm2_->scmat(spin));
}

RefSymmSCMatrix sc::PT2R12::lambda2(SpinCase2 spin)
{
  // since LocalSCMatrixKit is used everywhere, convert to Local kit
  return convert_to_local_kit(rdm2_->cumulant()->scmat(spin));
}

RefSymmSCMatrix sc::PT2R12::rdm2_gg(SpinCase2 spin)
{
  Ref<LocalSCMatrixKit> local_kit = new LocalSCMatrixKit;
  RefSymmSCMatrix result = local_kit->symmmatrix(r12eval_->dim_gg(spin));
  RefSymmSCMatrix rdm = this->rdm2(spin);

  const SpinCase1 spin1 = case1(spin);
  const SpinCase1 spin2 = case2(spin);
  Ref<OrbitalSpace> orbs1 = r12eval_->orbs(spin1);
  Ref<OrbitalSpace> orbs2 = r12eval_->orbs(spin2);
  Ref<OrbitalSpace> gspace1 = r12eval_->ggspace(spin1);
  Ref<OrbitalSpace> gspace2 = r12eval_->ggspace(spin2);
  MOIndexMap map1 = (*orbs1 << *gspace1);
  MOIndexMap map2 = (*orbs2 << *gspace2);
  SpinMOPairIter UV_iter(gspace1,gspace2,spin);
  SpinMOPairIter PQ_iter(gspace1,gspace2,spin);
  const int nmo = orbs1->rank();

  for(PQ_iter.start(); int(PQ_iter); PQ_iter.next()) {
    const int P = PQ_iter.i();
    const int Q = PQ_iter.j();
    const int PQ = PQ_iter.ij();
    const int pp = map1[P];
    const int qq = map2[Q];
    const int pq = ( (spin == AlphaBeta) ? (pp * nmo + qq) : (pp * (pp-1)/2 + qq));

    for(UV_iter.start(); int(UV_iter); UV_iter.next()) {
      const int U = UV_iter.i();
      const int V = UV_iter.j();
      const int UV = UV_iter.ij();
      const int uu = map1[U];
      const int vv = map2[V];
      const int uv = ( (spin == AlphaBeta) ? (uu * nmo + vv) : (uu * (uu-1)/2 + vv) );

      // first term is easy
      const double rdm_PQ_UV = rdm.get_element(pq, uv);
      result.set_element(PQ, UV, rdm_PQ_UV);
    }
  }

#if 0
  rdm.print(prepend_spincase(spin, "2-rdm (full)").c_str());
  result.print(prepend_spincase(spin, "2-rdm (occ)").c_str());
#endif

  return result;
}

RefSymmSCMatrix sc::PT2R12::density()
{
  throw FeatureNotImplemented("PT2R12::density() not yet implemented");
}

void sc::PT2R12::print(std::ostream & os)
{
  os << indent << "PT2R12:" << endl;
  os << incindent;
  reference_->print(os);
  r12world()->print(os);
  os << decindent;
}

void sc::PT2R12::compute()
{
  double energy_correction_r12 = 0.0;
  double energy_pt2r12[NSpinCases2];
  const bool spin_polarized = r12world()->ref()->spin_polarized();
  for(int i=0; i<NSpinCases2; i++) {
    SpinCase2 pairspin = static_cast<SpinCase2>(i);
    double scale = 1.0;
    if (pairspin == BetaBeta && !spin_polarized) continue;
    if (pairspin == AlphaAlpha && !spin_polarized)
      scale = 2.0;
    switch (r12world()->r12tech()->ansatz()->projector()) {
      case LinearR12::Projector_1:
        energy_pt2r12[i] = scale * energy_PT2R12_projector1(pairspin);
        break;
      case LinearR12::Projector_2:
        energy_pt2r12[i] = scale * energy_PT2R12_projector2(pairspin);
        break;
      default:
        abort();
    }
    energy_correction_r12 +=  energy_pt2r12[i];
  }

  const double energy = reference_->energy() + energy_correction_r12;
  ExEnv::out0() << indent << "PT2R12 energy correction: " << setprecision(20) << energy_correction_r12 << endl;
  ExEnv::out0() << indent << "reference energy: " << setprecision(20) << reference_->energy() << endl;
  ExEnv::out0() << indent << "total energy: " << setprecision(20) << energy << endl;

  set_energy(energy);
}

int sc::PT2R12::nelectron()
{
  return reference_->nelectron();
}

int sc::PT2R12::spin_polarized()
{
  return reference_->spin_polarized();
}

double PT2R12::compute_energy(const RefSCMatrix &hmat,
                              SpinCase2 pairspin,
                              bool print_pair_energies) {
  SpinCase1 spin1 = case1(pairspin);
  SpinCase1 spin2 = case2(pairspin);
  Ref<OrbitalSpace> gg1space = r12eval_->ggspace(spin1);
  Ref<OrbitalSpace> gg2space = r12eval_->ggspace(spin2);
  SpinMOPairIter gg_iter(gg1space,gg2space,pairspin);
  double energy = 0.0;
  if (print_pair_energies) ExEnv::out0() << indent
    << prepend_spincase(pairspin, " [2]_R12 pair energies:") << endl;
  for(gg_iter.start(); int(gg_iter); gg_iter.next()) {
    int i = gg_iter.i();
    int j = gg_iter.j();
    int ij = gg_iter.ij();
    if (print_pair_energies)
      ExEnv::out0() << setw(6) << i << setw(6) << j
                    << setw(20) << setprecision(12) << hmat.get_element(ij,ij) << endl;
    if((i>=nfzc_) && (j>=nfzc_)) {
      energy+=hmat.get_element(ij,ij);
    }
  }
  return energy;
}


/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:

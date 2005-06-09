//
// r12int_eval.cc
//
// Copyright (C) 2004 Edward Valeev
//
// Author: Edward Valeev <edward.valeev@chemistry.gatech.edu>
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

#include <util/misc/formio.h>
#include <util/ref/ref.h>
#include <util/state/state_bin.h>
#include <math/scmat/local.h>
#include <chemistry/qc/mbptr12/vxb_eval_info.h>
#include <chemistry/qc/mbptr12/pairiter.h>
#include <chemistry/qc/mbptr12/r12int_eval.h>

using namespace std;
using namespace sc;

#define TEST_FOCK 0

inline int max(int a,int b) { return (a > b) ? a : b;}

/*-----------------
  R12IntEval
 -----------------*/
static ClassDesc R12IntEval_cd(
  typeid(R12IntEval),"R12IntEval",1,"virtual public SavableState",
  0, 0, 0);

R12IntEval::R12IntEval(const Ref<R12IntEvalInfo>& r12info, bool gbc, bool ebc,
                       LinearR12::ABSMethod abs_method,
                       LinearR12::StandardApproximation stdapprox) :
  r12info_(r12info), gbc_(gbc), ebc_(ebc), abs_method_(abs_method),
  stdapprox_(stdapprox), spinadapted_(false), include_mp1_(false), evaluated_(false),
  debug_(0)
{
    int nocc_act = r12info_->nocc_act();
    int nvir_act = r12info_->nvir_act();
    dim_ij_aa_ = new SCDimension((nocc_act*(nocc_act-1))/2);
    dim_ij_ab_ = new SCDimension(nocc_act*nocc_act);
    dim_ij_s_ = new SCDimension((nocc_act*(nocc_act+1))/2);
    dim_ij_t_ = dim_ij_aa_;
    dim_ab_aa_ = new SCDimension((nvir_act*(nvir_act-1))/2);
    dim_ab_ab_ = new SCDimension(nvir_act*nvir_act);

    Ref<LocalSCMatrixKit> local_matrix_kit = new LocalSCMatrixKit();
    Vaa_ = local_matrix_kit->matrix(dim_ij_aa_,dim_ij_aa_);
    Vab_ = local_matrix_kit->matrix(dim_ij_ab_,dim_ij_ab_);
    Xaa_ = local_matrix_kit->matrix(dim_ij_aa_,dim_ij_aa_);
    Xab_ = local_matrix_kit->matrix(dim_ij_ab_,dim_ij_ab_);
    Baa_ = local_matrix_kit->matrix(dim_ij_aa_,dim_ij_aa_);
    Bab_ = local_matrix_kit->matrix(dim_ij_ab_,dim_ij_ab_);
    if (ebc_ == false) {
      Aaa_ = local_matrix_kit->matrix(dim_ij_aa_,dim_ab_aa_);
      Aab_ = local_matrix_kit->matrix(dim_ij_ab_,dim_ab_ab_);
      Ac_aa_ = local_matrix_kit->matrix(dim_ij_aa_,dim_ab_aa_);
      Ac_ab_ = local_matrix_kit->matrix(dim_ij_ab_,dim_ab_ab_);
      T2aa_ = local_matrix_kit->matrix(dim_ij_aa_,dim_ab_aa_);
      T2ab_ = local_matrix_kit->matrix(dim_ij_ab_,dim_ab_ab_);
      Raa_ = local_matrix_kit->matrix(dim_ij_aa_,dim_ab_aa_);
      Rab_ = local_matrix_kit->matrix(dim_ij_ab_,dim_ab_ab_);
    }
    emp2pair_aa_ = local_matrix_kit->vector(dim_ij_aa_);
    emp2pair_ab_ = local_matrix_kit->vector(dim_ij_ab_);
    
    init_intermeds_();
    init_tforms_();
}

R12IntEval::R12IntEval(StateIn& si) : SavableState(si)
{
  int gbc; si.get(gbc); gbc_ = (bool) gbc;
  int ebc; si.get(ebc); ebc_ = (bool) ebc;
  int absmethod; si.get(absmethod); abs_method_ = (LinearR12::ABSMethod) absmethod;
  int stdapprox; si.get(stdapprox); stdapprox_ = (LinearR12::StandardApproximation) stdapprox;

  r12info_ << SavableState::restore_state(si);
  dim_ij_aa_ << SavableState::restore_state(si);
  dim_ij_ab_ << SavableState::restore_state(si);
  dim_ij_s_ << SavableState::restore_state(si);
  dim_ij_t_ << SavableState::restore_state(si);
  dim_ab_aa_ << SavableState::restore_state(si);
  dim_ab_ab_ << SavableState::restore_state(si);

  Ref<LocalSCMatrixKit> local_matrix_kit = new LocalSCMatrixKit();
  Vaa_ = local_matrix_kit->matrix(dim_ij_aa_,dim_ij_aa_);
  Vab_ = local_matrix_kit->matrix(dim_ij_ab_,dim_ij_ab_);
  Xaa_ = local_matrix_kit->matrix(dim_ij_aa_,dim_ij_aa_);
  Xab_ = local_matrix_kit->matrix(dim_ij_ab_,dim_ij_ab_);
  Baa_ = local_matrix_kit->matrix(dim_ij_aa_,dim_ij_aa_);
  Bab_ = local_matrix_kit->matrix(dim_ij_ab_,dim_ij_ab_);
  if (ebc_ == false) {
    Aaa_ = local_matrix_kit->matrix(dim_ij_aa_,dim_ab_aa_);
    Aab_ = local_matrix_kit->matrix(dim_ij_ab_,dim_ab_ab_);
    T2aa_ = local_matrix_kit->matrix(dim_ij_aa_,dim_ab_aa_);
    T2ab_ = local_matrix_kit->matrix(dim_ij_ab_,dim_ab_ab_);
    Raa_ = local_matrix_kit->matrix(dim_ij_aa_,dim_ab_aa_);
    Rab_ = local_matrix_kit->matrix(dim_ij_ab_,dim_ab_ab_);
  }
  emp2pair_aa_ = local_matrix_kit->vector(dim_ij_aa_);
  emp2pair_ab_ = local_matrix_kit->vector(dim_ij_ab_);

  Vaa_.restore(si);
  Vab_.restore(si);
  Xaa_.restore(si);
  Xab_.restore(si);
  Baa_.restore(si);
  Bab_.restore(si);
  if (ebc_ == false) {
    Aaa_.restore(si);
    Aab_.restore(si);
    T2aa_.restore(si);
    T2ab_.restore(si);
    Raa_.restore(si);
    Rab_.restore(si);
  }
  emp2pair_aa_.restore(si);
  emp2pair_ab_.restore(si);

  int num_tforms;
  si.get(num_tforms);
  for(int t=0; t<num_tforms; t++) {
    std::string tform_name;
    si.get(tform_name);
    Ref<TwoBodyMOIntsTransform> tform;
    tform << SavableState::restore_state(si);
    tform_map_[tform_name] = tform;
  }

  int spinadapted; si.get(spinadapted); spinadapted_ = (bool) spinadapted;
  int evaluated; si.get(evaluated); evaluated_ = (bool) evaluated;
  si.get(debug_);

  init_tforms_();
}

R12IntEval::~R12IntEval()
{
  r12info_ = 0;
  dim_ij_aa_ = 0;
  dim_ij_ab_ = 0;
  dim_ij_s_ = 0;
  dim_ij_t_ = 0;
  dim_ab_aa_ = 0;
  dim_ab_ab_ = 0;
}

void
R12IntEval::save_data_state(StateOut& so)
{
  so.put((int)gbc_);
  so.put((int)ebc_);
  so.put((int)abs_method_);
  so.put((int)stdapprox_);

  SavableState::save_state(r12info_.pointer(),so);
  SavableState::save_state(dim_ij_aa_.pointer(),so);
  SavableState::save_state(dim_ij_ab_.pointer(),so);
  SavableState::save_state(dim_ij_s_.pointer(),so);
  SavableState::save_state(dim_ij_t_.pointer(),so);
  SavableState::save_state(dim_ab_aa_.pointer(),so);
  SavableState::save_state(dim_ab_ab_.pointer(),so);

  Vaa_.save(so);
  Vab_.save(so);
  Xaa_.save(so);
  Xab_.save(so);
  Baa_.save(so);
  Bab_.save(so);
  if (ebc_ == false) {
    Aaa_.save(so);
    Aab_.save(so);
    T2aa_.save(so);
    T2ab_.save(so);
    Raa_.save(so);
    Rab_.save(so);
  }
  emp2pair_aa_.save(so);
  emp2pair_ab_.save(so);

  int num_tforms = tform_map_.size();
  so.put(num_tforms);
  TformMap::iterator first_tform = tform_map_.begin();
  TformMap::iterator last_tform = tform_map_.end();
  for(TformMap::iterator t=first_tform; t!=last_tform; t++) {
    so.put((*t).first);
    SavableState::save_state((*t).second.pointer(),so);
  }

  so.put((int)spinadapted_);
  so.put((int)evaluated_);
  so.put(debug_);
}

void
R12IntEval::obsolete()
{
  evaluated_ = false;

  // make all transforms obsolete
  TformMap::iterator first_tform = tform_map_.begin();
  TformMap::iterator last_tform = tform_map_.end();
  for(TformMap::iterator t=first_tform; t!=last_tform; t++) {
    (*t).second->obsolete();
  }

  init_intermeds_();
}

void R12IntEval::include_mp1(bool include_mp1) { include_mp1_ = include_mp1; };
void R12IntEval::set_debug(int debug) { if (debug >= 0) { debug_ = debug; r12info_->set_debug_level(debug_); }};
void R12IntEval::set_dynamic(bool dynamic) { r12info_->set_dynamic(dynamic); };
void R12IntEval::set_print_percent(double pp) { r12info_->set_print_percent(pp); };
void R12IntEval::set_memory(size_t nbytes) { r12info_->set_memory(nbytes); };

Ref<R12IntEvalInfo> R12IntEval::r12info() const { return r12info_; };
RefSCDimension R12IntEval::dim_oo_aa() const { return dim_ij_aa_; };
RefSCDimension R12IntEval::dim_oo_ab() const { return dim_ij_ab_; };
RefSCDimension R12IntEval::dim_oo_s() const { return dim_ij_s_; };
RefSCDimension R12IntEval::dim_oo_t() const { return dim_ij_t_; };
RefSCDimension R12IntEval::dim_vv_aa() const { return dim_ab_aa_; };
RefSCDimension R12IntEval::dim_vv_ab() const { return dim_ab_ab_; };

RefSCMatrix R12IntEval::V_aa()
{
  compute();
  return Vaa_;
}

RefSCMatrix R12IntEval::X_aa()
{
  compute();
  return Xaa_;
}

RefSymmSCMatrix R12IntEval::B_aa()
{
  compute();

  // Extract lower triangle of the matrix
  Ref<SCMatrixKit> kit = Baa_.kit();
  RefSymmSCMatrix Baa = kit->symmmatrix(Baa_.rowdim());
  int naa = Baa_.nrow();
  double* baa = new double[naa*naa];
  Baa_.convert(baa);
  const double* baa_ptr = baa;
  for(int i=0; i<naa; i++, baa_ptr += i)
    for(int j=i; j<naa; j++, baa_ptr++)
      Baa.set_element(i,j,*baa_ptr);
  delete[] baa;

  return Baa;
}

RefSCMatrix R12IntEval::A_aa()
{
  if (ebc_ == false)
    compute();
  return Aaa_;
}

RefSCMatrix R12IntEval::Ac_aa()
{
  if (ebc_ == false)
    compute();
  return Ac_aa_;
}

RefSCMatrix R12IntEval::T2_aa()
{
  if (ebc_ == false)
    compute();
  return T2aa_;
}

RefSCMatrix R12IntEval::V_ab()
{
  compute();
  return Vab_;
}

RefSCMatrix R12IntEval::X_ab()
{
  compute();
  return Xab_;
}

RefSymmSCMatrix R12IntEval::B_ab()
{
  compute();

  // Extract lower triangle of the matrix
  Ref<SCMatrixKit> kit = Bab_.kit();
  RefSymmSCMatrix Bab = kit->symmmatrix(Bab_.rowdim());
  int nab = Bab_.nrow();
  double* bab = new double[nab*nab];
  Bab_.convert(bab);
  const double* bab_ptr = bab;
  for(int i=0; i<nab; i++, bab_ptr += i)
    for(int j=i; j<nab; j++, bab_ptr++)
      Bab.set_element(i,j,*bab_ptr);
  delete[] bab;

  return Bab;
}

RefSCMatrix R12IntEval::A_ab()
{
  if (ebc_ == false)
    compute();
  return Aab_;
}

RefSCMatrix R12IntEval::Ac_ab()
{
  if (ebc_ == false)
    compute();
  return Ac_ab_;
}

RefSCMatrix R12IntEval::T2_ab()
{
  if (ebc_ == false)
    compute();
  return T2ab_;
}

RefSCVector R12IntEval::emp2_aa()
{
  compute();
  return emp2pair_aa_;
}

RefSCVector R12IntEval::emp2_ab()
{
  compute();
  return emp2pair_ab_;
}

RefDiagSCMatrix R12IntEval::evals() const { return r12info_->obs_space()->evals(); };

void
R12IntEval::checkpoint_() const
{
  int me = r12info_->msg()->me();
  Wavefunction* wfn = r12info_->wfn();

  if (me == 0 && wfn->if_to_checkpoint()) {
    StateOutBin stateout(wfn->checkpoint_file());
    SavableState::save_state(wfn,stateout);
    ExEnv::out0() << indent << "Checkpointed the wave function" << endl;
  }
}

void
R12IntEval::init_tforms_()
{
  Ref<MOIntsTransformFactory> tfactory = r12info_->tfactory();
  tfactory->set_ints_method((MOIntsTransformFactory::StoreMethod)r12info_->ints_method());

  const std::string ipjq_name = "(ip|jq)";
  Ref<TwoBodyMOIntsTransform> ipjq_tform = tform_map_[ipjq_name];
  if (ipjq_tform.null()) {
    tfactory->set_spaces(r12info_->act_occ_space(),r12info_->obs_space(),
                         r12info_->act_occ_space(),r12info_->obs_space());
    ipjq_tform = tfactory->twobody_transform_13(ipjq_name);
    tform_map_[ipjq_name] = ipjq_tform;
    tform_map_[ipjq_name]->set_num_te_types(3);
  }

  const std::string iajb_name = "(ia|jb)";
  Ref<TwoBodyMOIntsTransform> iajb_tform = tform_map_[iajb_name];
  if (iajb_tform.null()) {
    tfactory->set_spaces(r12info_->act_occ_space(),r12info_->act_vir_space(),
                         r12info_->act_occ_space(),r12info_->act_vir_space());
    iajb_tform = tfactory->twobody_transform_13(iajb_name);
    tform_map_[iajb_name] = iajb_tform;
    tform_map_[iajb_name]->set_num_te_types(3);
  }

  const std::string imja_name = "(im|ja)";
  Ref<TwoBodyMOIntsTransform> imja_tform = tform_map_[imja_name];
  if (imja_tform.null()) {
    tfactory->set_spaces(r12info_->act_occ_space(),r12info_->occ_space(),
                         r12info_->act_occ_space(),r12info_->act_vir_space());
    imja_tform = tfactory->twobody_transform_13(imja_name);
    tform_map_[imja_name] = imja_tform;
    tform_map_[imja_name]->set_num_te_types(4);
  }

  const std::string imjn_name = "(im|jn)";
  Ref<TwoBodyMOIntsTransform> imjn_tform = tform_map_[imjn_name];
  if (imjn_tform.null()) {
    tfactory->set_spaces(r12info_->act_occ_space(),r12info_->occ_space(),
                         r12info_->act_occ_space(),r12info_->occ_space());
    imjn_tform = tfactory->twobody_transform_13(imjn_name);
    tform_map_[imjn_name] = imjn_tform;
    tform_map_[imjn_name]->set_num_te_types(3);
  }

  const std::string imjy_name = "(im|jy)";
  Ref<TwoBodyMOIntsTransform> imjy_tform = tform_map_[imjy_name];
  if (imjy_tform.null()) {
    tfactory->set_spaces(r12info_->act_occ_space(),r12info_->occ_space(),
                         r12info_->act_occ_space(),r12info_->ribs_space());
    imjy_tform = tfactory->twobody_transform_13(imjy_name);
    tform_map_[imjy_name] = imjy_tform;
    tform_map_[imjy_name]->set_num_te_types(4);
  }

  iajb_tform = tform_map_[iajb_name];
  imjn_tform = tform_map_[imjn_name];
  ipjq_tform = tform_map_[ipjq_name];
}

Ref<TwoBodyMOIntsTransform>
R12IntEval::get_tform_(const std::string& tform_name)
{
  TformMap::const_iterator tform_iter = tform_map_.find(tform_name);
  TwoBodyMOIntsTransform* tform = (*tform_iter).second.pointer();
  if (tform == NULL) {
    std::string errmsg = "R12IntEval::get_tform_() -- transform " + tform_name + " is not known";
    throw std::runtime_error(errmsg.c_str());
  }
  tform->compute();

  return tform;
}

void
R12IntEval::init_intermeds_()
{
  if (r12info_->msg()->me() == 0) {
    Vaa_->unit();
    Vab_->unit();
    Baa_->unit();
    Bab_->unit();
  }
  else {
    Vaa_.assign(0.0);
    Vab_.assign(0.0);
    Baa_.assign(0.0);
    Bab_.assign(0.0);
  }
  if (ebc_ == false) {
    Aaa_.assign(0.0);
    Aab_.assign(0.0);
    Ac_aa_.assign(0.0);
    Ac_ab_.assign(0.0);
    T2aa_.assign(0.0);
    T2ab_.assign(0.0);
    Raa_.assign(0.0);
    Rab_.assign(0.0);
  }

  Xaa_.assign(0.0);
  Xab_.assign(0.0);
  //r2_contrib_to_X_orig_();
  r2_contrib_to_X_new_();

  emp2pair_aa_.assign(0.0);
  emp2pair_ab_.assign(0.0);
}

/// Compute <space1 space1|r_{12}^2|space1 space2>
RefSCMatrix
R12IntEval::compute_r2_(const Ref<MOIndexSpace>& space1, const Ref<MOIndexSpace>& space2)
{
  /*-----------------------------------------------------
    Compute overlap, dipole, quadrupole moment integrals
   -----------------------------------------------------*/
  RefSCMatrix S_11, MX_11, MY_11, MZ_11, MXX_11, MYY_11, MZZ_11;
  r12info_->compute_multipole_ints(space1, space1, MX_11, MY_11, MZ_11, MXX_11, MYY_11, MZZ_11);
  r12info_->compute_overlap_ints(space1, space1, S_11);

  RefSCMatrix S_12, MX_12, MY_12, MZ_12, MXX_12, MYY_12, MZZ_12;
  if (space1 == space2) {
    S_12 = S_11;
    MX_12 = MX_11;
    MY_12 = MY_11;
    MZ_12 = MZ_11;
    MXX_12 = MXX_11;
    MYY_12 = MYY_11;
    MZZ_12 = MZZ_11;
  }
  else {
    r12info_->compute_multipole_ints(space1, space2, MX_12, MY_12, MZ_12, MXX_12, MYY_12, MZZ_12);
    r12info_->compute_overlap_ints(space1, space2, S_12);
  }
  if (debug_)
    ExEnv::out0() << indent << "Computed overlap and multipole moment integrals" << endl;

  const int nproc = r12info_->msg()->n();
  const int me = r12info_->msg()->me();

  const int n1 = space1->rank();
  const int n2 = space2->rank();
  const int n12 = n1*n2;
  const int n1112 = n1*n1*n12;
  double* r2_array = new double[n1112];
  memset(r2_array,0,n1112*sizeof(double));

  int ij = 0;
  double* ijkl_ptr = r2_array;
  for(int i=0; i<n1; i++)
    for(int j=0; j<n1; j++, ij++) {

    int ij_proc = ij%nproc;
    if (ij_proc != me) {
      ijkl_ptr += n12;
      continue;
    }

    int kl=0;
    for(int k=0; k<n1; k++)
      for(int l=0; l<n2; l++, kl++, ijkl_ptr++) {

        double r1r1_ik = -1.0*(MXX_11->get_element(i,k) + MYY_11->get_element(i,k) + MZZ_11->get_element(i,k));
        double r1r1_jl = -1.0*(MXX_12->get_element(j,l) + MYY_12->get_element(j,l) + MZZ_12->get_element(j,l));
        double r1r2_ijkl = MX_11->get_element(i,k)*MX_12->get_element(j,l) +
          MY_11->get_element(i,k)*MY_12->get_element(j,l) +
          MZ_11->get_element(i,k)*MZ_12->get_element(j,l);
        double S_ik = S_11.get_element(i,k);
        double S_jl = S_12.get_element(j,l);
        
        double R2_ijkl = r1r1_ik * S_jl + r1r1_jl * S_ik - 2.0*r1r2_ijkl;
        *ijkl_ptr = R2_ijkl;
      }
    }

  r12info_->msg()->sum(r2_array,n1112);

  MOPairIterFactory pair_factory;
  RefSCDimension dim_ij = pair_factory.scdim_ab(space1,space1);
  RefSCDimension dim_kl = pair_factory.scdim_ab(space1,space2);

  Ref<LocalSCMatrixKit> local_matrix_kit = new LocalSCMatrixKit();
  RefSCMatrix R2 = local_matrix_kit->matrix(dim_ij, dim_kl);
  R2.assign(r2_array);
  delete[] r2_array;

  return R2;
}


void
R12IntEval::r2_contrib_to_X_orig_()
{
  /*---------------------------------------------------------------
    Compute dipole and quadrupole moment integrals in act MO basis
   ---------------------------------------------------------------*/
  RefSCMatrix MX, MY, MZ, MXX, MYY, MZZ;
  r12info_->compute_multipole_ints(r12info_->act_occ_space(),r12info_->act_occ_space(),MX,MY,MZ,MXX,MYY,MZZ);
  if (debug_)
    ExEnv::out0() << indent << "Computed multipole moment integrals" << endl;

  const int nproc = r12info_->msg()->n();
  const int me = r12info_->msg()->me();

  SpatialMOPairIter_eq ij_iter(r12info_->act_occ_space());
  SpatialMOPairIter_eq kl_iter(r12info_->act_occ_space());

  for(kl_iter.start();int(kl_iter);kl_iter.next()) {

    const int kl = kl_iter.ij();
    int kl_proc = kl%nproc;
    if (kl_proc != me)
      continue;
    const int k = kl_iter.i();
    const int l = kl_iter.j();
    const int kl_aa = kl_iter.ij_aa();
    const int kl_ab = kl_iter.ij_ab();
    const int lk_ab = kl_iter.ij_ba();

    for(ij_iter.start();int(ij_iter);ij_iter.next()) {

      const int i = ij_iter.i();
      const int j = ij_iter.j();
      const int ij_aa = ij_iter.ij_aa();
      const int ij_ab = ij_iter.ij_ab();
      const int ji_ab = ij_iter.ij_ba();

      /*----------------------------------
        Compute (r12)^2 contribution to X
       ----------------------------------*/
      double r1r1_ik = -1.0*(MXX->get_element(i,k) + MYY->get_element(i,k) + MZZ->get_element(i,k));
      double r1r1_il = -1.0*(MXX->get_element(i,l) + MYY->get_element(i,l) + MZZ->get_element(i,l));
      double r1r1_jk = -1.0*(MXX->get_element(j,k) + MYY->get_element(j,k) + MZZ->get_element(j,k));
      double r1r1_jl = -1.0*(MXX->get_element(j,l) + MYY->get_element(j,l) + MZZ->get_element(j,l));
      double r1r2_ijkl = MX->get_element(i,k)*MX->get_element(j,l) +
        MY->get_element(i,k)*MY->get_element(j,l) +
        MZ->get_element(i,k)*MZ->get_element(j,l);
      double r1r2_ijlk = MX->get_element(i,l)*MX->get_element(j,k) +
        MY->get_element(i,l)*MY->get_element(j,k) +
        MZ->get_element(i,l)*MZ->get_element(j,k);
      double delta_ik = (i==k ? 1.0 : 0.0);
      double delta_il = (i==l ? 1.0 : 0.0);
      double delta_jk = (j==k ? 1.0 : 0.0);
      double delta_jl = (j==l ? 1.0 : 0.0);
      
      double Xab_ijkl = r1r1_ik * delta_jl + r1r1_jl * delta_ik - 2.0*r1r2_ijkl;
      Xab_.accumulate_element(ij_ab,kl_ab,Xab_ijkl);

      if (ij_ab != ji_ab) {
        double Xab_jikl = r1r1_jk * delta_il + r1r1_il * delta_jk - 2.0*r1r2_ijlk;
        Xab_.accumulate_element(ji_ab,kl_ab,Xab_jikl);
      }

      if (kl_ab != lk_ab) {
        double Xab_ijlk = r1r1_il * delta_jk + r1r1_jk * delta_il - 2.0*r1r2_ijlk;
        Xab_.accumulate_element(ij_ab,lk_ab,Xab_ijlk);
      }

      if (ij_ab != ji_ab && kl_ab != lk_ab) {
        double Xab_jilk = r1r1_ik * delta_jl + r1r1_jl * delta_ik - 2.0*r1r2_ijkl;
        Xab_.accumulate_element(ji_ab,lk_ab,Xab_jilk);
      }

      if (ij_aa != -1 && kl_aa != -1) {
        double Xaa_ijkl = r1r1_ik * delta_jl + r1r1_jl * delta_ik - 2.0*r1r2_ijkl -
          r1r1_jk * delta_il - r1r1_il * delta_jk + 2.0*r1r2_ijlk;
        Xaa_.accumulate_element(ij_aa,kl_aa,Xaa_ijkl);
      }

    }
  }          
}

void
R12IntEval::r2_contrib_to_X_new_()
{
  unsigned int me = r12info_->msg()->me();

  // compute r_{12}^2 operator in act.occ.pair/act.occ.pair basis
  RefSCMatrix R2 = compute_r2_(r12info_->act_occ_space(),r12info_->act_occ_space());

  if (me != 0)
    return;
  Xab_.accumulate(R2);

  SpatialMOPairIter_eq ij_iter(r12info_->act_occ_space());
  SpatialMOPairIter_eq kl_iter(r12info_->act_occ_space());

  for(kl_iter.start();int(kl_iter);kl_iter.next()) {

    const int kl_aa = kl_iter.ij_aa();
    if (kl_aa == -1)
      continue;
    const int kl_ab = kl_iter.ij_ab();

    for(ij_iter.start();int(ij_iter);ij_iter.next()) {

      const int ij_aa = ij_iter.ij_aa();
      const int ij_ab = ij_iter.ij_ab();
      const int ji_ab = ij_iter.ij_ba();

      if (ij_aa != -1) {
        double Xaa_ijkl = R2.get_element(ij_ab,kl_ab) - R2.get_element(ji_ab,kl_ab);
        Xaa_.accumulate_element(ij_aa,kl_aa,Xaa_ijkl);
      }

    }
  }
}

void
R12IntEval::form_focc_space_()
{
  // compute the Fock matrix between the complement and all occupieds and
  // create the new Fock-weighted space
  if (focc_space_.null()) {
    Ref<MOIndexSpace> occ_space = r12info_->occ_space();
    Ref<MOIndexSpace> ribs_space = r12info_->ribs_space();
    
    RefSCMatrix F_ri_o = fock_(occ_space,ribs_space,occ_space);
    if (debug_ > 1)
      F_ri_o.print("Fock matrix (RI-BS/occ.)");
    focc_space_ = new MOIndexSpace("Fock-weighted occupied MOs sorted by energy",
                                   occ_space, ribs_space->coefs()*F_ri_o, ribs_space->basis());
  }
}

void
R12IntEval::form_factocc_space_()
{
  // compute the Fock matrix between the complement and active occupieds and
  // create the new Fock-weighted space
  if (factocc_space_.null()) {
    Ref<MOIndexSpace> occ_space = r12info_->occ_space();
    Ref<MOIndexSpace> act_occ_space = r12info_->act_occ_space();
    Ref<MOIndexSpace> ribs_space = r12info_->ribs_space();
    
    RefSCMatrix F_ri_ao = fock_(occ_space,ribs_space,act_occ_space);
    if (debug_ > 1)
      F_ri_ao.print("Fock matrix (RI-BS/act.occ.)");
    factocc_space_ = new MOIndexSpace("Fock-weighted active occupied MOs sorted by energy",
                                      act_occ_space, ribs_space->coefs()*F_ri_ao, ribs_space->basis());
  }
}

void
R12IntEval::form_canonvir_space_()
{
  // Create a complement space to all occupieds
  // Fock operator is diagonal in this space
  if (canonvir_space_.null()) {

    if (r12info_->basis_vir()->equiv(r12info_->basis())) {
      canonvir_space_ = r12info_->vir_space();
      return;
    }

    const Ref<MOIndexSpace> mo_space = r12info_->mo_space();
    Ref<MOIndexSpace> vir_space = r12info_->vir_space_symblk();
    RefSCMatrix F_vir = fock_(r12info_->occ_space(),vir_space,vir_space);

    int nrow = vir_space->rank();
    double* F_full = new double[nrow*nrow];
    double* F_lowtri = new double [nrow*(nrow+1)/2];
    F_vir->convert(F_full);
    int ij = 0;
    for(int row=0; row<nrow; row++) {
      int rc = row*nrow;
      for(int col=0; col<=row; col++, rc++, ij++) {
        F_lowtri[ij] = F_full[rc];
      }
    }
    RefSymmSCMatrix F_vir_lt(F_vir.rowdim(),F_vir->kit());
    F_vir_lt->assign(F_lowtri);
    F_vir = 0;
    delete[] F_full;
    delete[] F_lowtri;

    Ref<MOIndexSpace> canonvir_space_symblk = new MOIndexSpace("Virt. MOs symmetry-blocked",
                                                               vir_space, vir_space->coefs()*F_vir_lt.eigvecs(),
                                                               vir_space->basis());
    RefDiagSCMatrix F_vir_evals = F_vir_lt.eigvals();
    canonvir_space_ = new MOIndexSpace("Virt. MOs sorted by energy",
                                       canonvir_space_symblk->coefs(), canonvir_space_symblk->basis(),
                                       F_vir_evals, 0, 0,
                                       MOIndexSpace::energy);
  }
}

const int
R12IntEval::tasks_with_ints_(const Ref<R12IntsAcc> ints_acc, vector<int>& map_to_twi)
{
  int nproc = r12info_->msg()->n();
  
  // Compute the number of tasks that have full access to the integrals
  // and split the work among them
  int nproc_with_ints = 0;
  for(int proc=0;proc<nproc;proc++)
    if (ints_acc->has_access(proc)) nproc_with_ints++;
 
  map_to_twi.resize(nproc);
  int count = 0;
  for(int proc=0;proc<nproc;proc++)
    if (ints_acc->has_access(proc)) {
      map_to_twi[proc] = count;
      count++;
    }
      else
        map_to_twi[proc] = -1;

  ExEnv::out0() << indent << "Computing intermediates on " << nproc_with_ints
    << " processors" << endl;
    
  return nproc_with_ints;
}


void
R12IntEval::compute()
{
  if (evaluated_)
    return;
  
  if (r12info_->basis_vir()->equiv(r12info_->basis())) {
    obs_contrib_to_VXB_gebc_vbseqobs_();
    if (debug_ > 1) {
      Baa_.print("Alpha-alpha B(OBS) contribution");
      Bab_.print("Alpha-beta B(OBS) contribution");
    }
    if (r12info_->basis() != r12info_->basis_ri())
      abs1_contrib_to_VXB_gebc_();
    if (debug_ > 1) {
      Baa_.print("Alpha-alpha B(OBS+ABS) contribution");
      Bab_.print("Alpha-beta B(OBS+ABS) contribution");
    }
  }
  else {
    contrib_to_VXB_gebc_vbsneqobs_();
    compute_dualEmp2_();
    if (include_mp1_)
      compute_dualEmp1_();
  }

  
#if TEST_FOCK
  if (!evaluated_) {
    RefSCMatrix F = fock_(r12info_->occ_space(),r12info_->obs_space(),r12info_->obs_space());
    F.print("Fock matrix in OBS");
    r12info_->obs_space()->evals().print("OBS eigenvalues");

    r12info_->ribs_space()->coefs().print("Orthonormal RI-BS");
    RefSCMatrix S_ri;
    r12info_->compute_overlap_ints(r12info_->ribs_space(),r12info_->ribs_space(),S_ri);
    S_ri.print("Overlap in RI-BS");
    RefSCMatrix F_ri = fock_(r12info_->occ_space(),r12info_->ribs_space(),r12info_->ribs_space());
    F_ri.print("Fock matrix in RI-BS");
    RefSymmSCMatrix F_ri_symm = F_ri.kit()->symmmatrix(F_ri.rowdim());
    int nrow = F_ri.rowdim().n();
    for(int r=0; r<nrow; r++)
      for(int c=0; c<nrow; c++)
        F_ri_symm.set_element(r,c,F_ri.get_element(r,c));
    F_ri_symm.eigvals().print("Eigenvalues of the Fock matrix (RI-BS)");

    RefSCMatrix F_obs_ri = fock_(r12info_->occ_space(),r12info_->obs_space(),r12info_->ribs_space());
    F_obs_ri.print("Fock matrix in OBS/RI-BS");
  }
#endif

  if (!ebc_) {
    // These functions assume that virtuals are expanded in the same basis
    // as the occupied orbitals
    if (!r12info_->basis_vir()->equiv(r12info_->basis()))
      throw std::runtime_error("R12IntEval::compute() -- ebc=false is only supported when basis_vir == basis");

    compute_A_simple_();
    compute_A_via_commutator_();
    //Aaa_.assign(Ac_aa_);
    //Aab_.assign(Ac_ab_);
    compute_T2_();
    AT2_contrib_to_V_();
    compute_R_();
    AR_contrib_to_B_();
  }
  
  if (!gbc_) {
    // These functions assume that virtuals are expanded in the same basis
    // as the occupied orbitals
    if (!r12info_->basis_vir()->equiv(r12info_->basis()))
      throw std::runtime_error("R12IntEval::compute() -- gbc=false is only supported when basis_vir == basis");

    compute_B_gbc_1_();
    if (debug_ > 1) {
      Baa_.print("Alpha-alpha B(OBS+ABS+GBC1) contribution");
      Bab_.print("Alpha-beta B(OBS+ABS+GBC1) contribution");
    }
    compute_B_gbc_2_();
    if (debug_ > 1) {
      Baa_.print("Alpha-alpha B(OBS+ABS+GBC1+GBC2) contribution");
      Bab_.print("Alpha-beta B(OBS+ABS+GBC1+GBC2) contribution");
    }
  }

  // Distribute the final intermediates to every node
  globally_sum_intermeds_(true);

  evaluated_ = true;
}

void
R12IntEval::globally_sum_scmatrix_(RefSCMatrix& A, bool to_all_tasks, bool to_average)
{
  Ref<MessageGrp> msg = r12info_->msg();
  unsigned int ntasks = msg->n();
  // If there's only one task then there's nothing to do
  if (ntasks == 1)
    return;

  const int nelem = A.ncol() * A.nrow();
  double *A_array = new double[nelem];
  A.convert(A_array);
  if (to_all_tasks)
    msg->sum(A_array,nelem,0,-1);
  else
    msg->sum(A_array,nelem,0,0);
  A.assign(A_array);
  if (to_average)
    A.scale(1.0/(double)ntasks);
  if (!to_all_tasks && msg->me() != 0)
    A.assign(0.0);

  delete[] A_array;
}

void
R12IntEval::globally_sum_scvector_(RefSCVector& A, bool to_all_tasks, bool to_average)
{
  Ref<MessageGrp> msg = r12info_->msg();
  unsigned int ntasks = msg->n();
  // If there's only one task then there's nothing to do
  if (ntasks == 1)
    return;

  const int nelem = A.dim().n();
  double *A_array = new double[nelem];
  A.convert(A_array);
  if (to_all_tasks)
    msg->sum(A_array,nelem,0,-1);
  else
    msg->sum(A_array,nelem,0,0);
  A.assign(A_array);
  if (to_average)
    A.scale(1.0/(double)ntasks);
  if (!to_all_tasks && msg->me() != 0)
    A.assign(0.0);

  delete[] A_array;
}

void
R12IntEval::globally_sum_intermeds_(bool to_all_tasks)
{
  globally_sum_scmatrix_(Vaa_,to_all_tasks);
  globally_sum_scmatrix_(Vab_,to_all_tasks);

  globally_sum_scmatrix_(Xaa_,to_all_tasks);
  globally_sum_scmatrix_(Xab_,to_all_tasks);

  globally_sum_scmatrix_(Baa_,to_all_tasks);
  globally_sum_scmatrix_(Bab_,to_all_tasks);

  if (ebc_ == false) {
    globally_sum_scmatrix_(Aaa_,to_all_tasks);
    globally_sum_scmatrix_(Aab_,to_all_tasks);

    globally_sum_scmatrix_(Ac_aa_,to_all_tasks);
    globally_sum_scmatrix_(Ac_ab_,to_all_tasks);
    
    globally_sum_scmatrix_(T2aa_,to_all_tasks);
    globally_sum_scmatrix_(T2ab_,to_all_tasks);
    
    globally_sum_scmatrix_(Raa_,to_all_tasks);
    globally_sum_scmatrix_(Rab_,to_all_tasks);
  }

  globally_sum_scvector_(emp2pair_aa_,to_all_tasks);
  globally_sum_scvector_(emp2pair_ab_,to_all_tasks);

  if (debug_) {
    ExEnv::out0() << indent << "Collected contributions to the intermediates from all tasks";
    if (to_all_tasks)
      ExEnv::out0() << " and distributed to every task" << endl;
    else
      ExEnv::out0() << " on task 0" << endl;
  }
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:

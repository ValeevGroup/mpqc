//
// ccr12_info_fill_in.cc
//
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki.toru@gmail.com>
// Maintainer: EFV & TS
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

#include <string>
#include <cassert>
#include <vector>
#include <algorithm>
#include <chemistry/qc/wfn/orbitalspace.h>
#include <chemistry/qc/ccr12/ccr12_info.h>
#include <chemistry/qc/ccr12/mtensor.h>

using namespace std;
using namespace sc;

void CCR12_Info::compute_corr_space() {

  const std::string id("pp(corr)");
  const std::string label("Active reference spin-orbitals in correlated order");

  // Prepare OrbitalSpace that corresponds to the CorrelatedMOOrder assumed by SMITH
  Ref<OrbitalSpace> occs_a = r12world()->refwfn()->occ_sb(Alpha);
  Ref<OrbitalSpace> occs_b = r12world()->refwfn()->occ_sb(Beta);
  Ref<OrbitalSpace> orbs_a = r12world()->refwfn()->orbs_sb(Alpha);
  Ref<OrbitalSpace> orbs_b = r12world()->refwfn()->orbs_sb(Beta);
  const size_t norbs = orbs_a->rank();
  RefDiagSCMatrix occnums_a = orbs_a->evals().clone();
  occnums_a.assign(0.0);
  RefDiagSCMatrix occnums_b = orbs_b->evals().clone();
  occnums_b.assign(0.0);
  const unsigned int nirreps = occnums_a.dim()->blocks()->nblock();
  const std::vector<unsigned int>& occpi_a = occs_a->block_sizes();
  const std::vector<unsigned int>& occpi_b = occs_b->block_sizes();
  for (unsigned int irrep = 0; irrep < nirreps; ++irrep) {
    const unsigned int nocc_a = occpi_a[irrep];
    const unsigned int occ_a_begin = occnums_a.dim()->blocks()->start(irrep);
    const unsigned int occ_a_end = occ_a_begin + nocc_a;
    for (unsigned int o = occ_a_begin; o < occ_a_end; ++o)
      occnums_a( o) = 1.0;
    const unsigned int nocc_b = occpi_b[irrep];
    const unsigned int occ_b_begin = occnums_b.dim()->blocks()->start(irrep);
    const unsigned int occ_b_end = occ_b_begin + nocc_b;
    for (unsigned int o = occ_b_begin; o < occ_b_end; ++o)
      occnums_b( o) = 1.0;
  }

  Ref<OrbitalSpace>
      full_space =
          new OrderedSpinOrbitalSpace<CorrelatedSpinMOOrder> (
                                                              std::string(
                                                                          "pp(corr)"),
                                                              std::string(
                                                                          "Active reference spin-orbitals in correlated order"),
                                                              orbs_a->basis(),
                                                              orbs_a->integral(),
                                                              orbs_a->coefs(),
                                                              orbs_b->coefs(),
                                                              orbs_a->evals(),
                                                              orbs_b->evals(),
                                                              occnums_a,
                                                              occnums_b,
                                                              orbs_a->orbsym(),
                                                              orbs_b->orbsym(),
                                                              CorrelatedSpinMOOrder(
                                                                                    nirreps));
  corr_space_ = full_space;

  // minus frozen core
  {
    typedef MolecularOrbitalMask<double, RefDiagSCMatrix> FCMask;
    FCMask fcmask(2 * nfzc_, corr_space_->evals());
    Ref<OrbitalSpace>
        full_space_minus_fc =
            new MaskedOrbitalSpace(
                                   std::string("pp(corr)"),
                                   std::string(
                                               "Active reference spin-orbitals in correlated order"),
                                   corr_space_, fcmask.mask());
    corr_space_ = full_space_minus_fc;
  }
  // minus frozen virtuals
  if (nfzv_ != 0) {
    typedef MolecularOrbitalMask<double, RefDiagSCMatrix, std::greater<double> >
        FVMask;
    FVMask fvmask(2 * nfzv_, corr_space_->evals());
    Ref<OrbitalSpace>
        full_space_minus_fv =
            new MaskedOrbitalSpace(
                                   std::string("pp(corr)"),
                                   std::string(
                                               "Active reference spin-orbitals in correlated order"),
                                   corr_space_, fvmask.mask());
    corr_space_ = full_space_minus_fv;
  }

  // append CABS alpha and beta spaces
  {
    const bool merge_blocks = false;
    Ref<OrbitalSpace> plus_cabs_a = new OrbitalSpaceUnion(id,label,
                                                          *(corr_space_),*(r12world()->cabs_space(Alpha)),
                                                          merge_blocks);
    Ref<OrbitalSpace> plus_cabs_ab = new OrbitalSpaceUnion(id,label,
                                                          *(plus_cabs_a),*(r12world()->cabs_space(Beta)),
                                                          merge_blocks);
    corr_space_ = plus_cabs_ab;
  }

  // register this space so that I can use it with "smart" runtimes
  const Ref<OrbitalSpaceRegistry> idxreg = r12world()->world()->moints_runtime()->factory()->orbital_registry();
  idxreg->add(make_keyspace_pair(corr_space_));

  //
  // also compute CABS-free correlated spaces for both spins
  //
  {
    std::string id = ParsedOrbitalSpaceKey::key("p(corr)", restricted_ ? AnySpinCase1 : Alpha);
    std::string name("Alpha MOs in correlated order");
    aobs_space_ =
          new OrderedOrbitalSpace<CorrelatedMOOrder> (id, name,
              orbs_a->basis(),
              orbs_a->integral(),
              orbs_a->coefs(),
              orbs_a->evals(),
              occnums_a,
              orbs_a->orbsym(),
              CorrelatedMOOrder(nirreps));
    // minus frozen core
    {
      typedef MolecularOrbitalMask<double, RefDiagSCMatrix> FCMask;
      FCMask fcmask(nfzc_, aobs_space_->evals());
      aobs_space_ =
            new MaskedOrbitalSpace(id,name,
                                   aobs_space_, fcmask.mask());
    }
    // minus frozen virtuals
    {
      typedef MolecularOrbitalMask<double, RefDiagSCMatrix, std::greater<double> > FVMask;
      FVMask fvmask(nfzv_, aobs_space_->evals());
      aobs_space_ =
            new MaskedOrbitalSpace(id, name,
                                   aobs_space_, fvmask.mask());
    }
    idxreg->add(make_keyspace_pair(aobs_space_));
  }

  if (!restricted_) {
    std::string id = ParsedOrbitalSpaceKey::key("p(corr)",Beta);
    std::string name("Beta MOs in correlated order");
    bobs_space_ =
          new OrderedOrbitalSpace<CorrelatedMOOrder> (id, name,
              orbs_b->basis(),
              orbs_b->integral(),
              orbs_b->coefs(),
              orbs_b->evals(),
              occnums_b,
              orbs_b->orbsym(),
              CorrelatedMOOrder(nirreps));
    // minus frozen core
    {
      typedef MolecularOrbitalMask<double, RefDiagSCMatrix> FCMask;
      FCMask fcmask(nfzc_, bobs_space_->evals());
      bobs_space_ =
            new MaskedOrbitalSpace(id,name,
                                   bobs_space_, fcmask.mask());
    }
    // minus frozen virtuals
    {
      typedef MolecularOrbitalMask<double, RefDiagSCMatrix, std::greater<double> > FVMask;
      FVMask fvmask(nfzv_, bobs_space_->evals());
      bobs_space_ =
            new MaskedOrbitalSpace(id, name,
                                   bobs_space_, fvmask.mask());
    }
    idxreg->add(make_keyspace_pair(bobs_space_));
  }
  else // spin-restricted
    bobs_space_ = aobs_space_;

}

void
CCR12_Info::compute_source_integrals_rhf_r12() {
  MPQC_ASSERT(need_w1());

  Ref<OrbitalSpace> occs = r12world()->refwfn()->occ_sb();
  Ref<OrbitalSpace> orbs = r12world()->refwfn()->orbs_sb();
  Ref<OrbitalSpace> cabs = r12world()->refwfn()->orbs_sb();

  abort();
}

void
CCR12_Info::compute_source_integrals_rhf() {

  // this function should only be called for rhf
  if (r12world()->refwfn()->spin_polarized())
    throw ProgrammingError("CCR12_Info::compute_source_integrals_rhf -- reference wfn is spin-polarized",
                           __FILE__,__LINE__);

  // and CCSD
  if (theory_ != std::string("CCSD") && theory_ != std::string("CCSDT") && theory_ != std::string("CCSDTQ")
   && theory_ != std::string("CCSD-R12") && theory_ != std::string("CCSD(R12)"))
    throw FeatureNotImplemented("Can only import integrals for CCSD, CCSD-R12, CCSDT, and CCSDTQ methods",
                                    __FILE__,__LINE__);

  // compute Fock matrices in pitzer-order spaces because their occ-virtual blocks may be nonzero
  Ref<OrbitalSpace> aobs_space = r12world()->refwfn()->orbs_sb(Alpha);
  Ref<FockBuildRuntime> fb_rtime = r12world()->world()->fockbuild_runtime();
  const std::string fkey = ParsedOneBodyIntKey::key(aobs_space->id(),aobs_space->id(),std::string("F"));
  RefSCMatrix F = fb_rtime->get(fkey);
  F_[Alpha] = F;
  F_[Beta] = F_[Alpha];

  const std::string skey = aobs_space_->id();
  const std::string
      tkey = ParsedTwoBodyFourCenterIntKey::key(skey, skey, skey, skey,
                                                std::string("ERI"), std::string(""),
                                                TwoBodyIntLayout::b1b2_k1k2);

  Ref<TwoBodyMOIntsTransform> pppp_tform =
      r12world()->world()->moints_runtime()->runtime_4c()->get(tkey);
  pppp_tform->compute();
  pppp_acc_[AlphaBeta] = pppp_tform->ints_distarray4();
  pppp_acc_[AlphaBeta]->activate();

  // more integrals are needed for the (R12) calculations
  if (this->need_w1() || this->need_w2()) {
    // get the CABS space
    Ref<OrbitalSpace> cabs_space = r12world()->cabs_space(Alpha);

    r12int_eval_ = new R12IntEval(r12world_);
    r12int_eval_->compute();

    Vgg_[AlphaBeta] = r12int_eval_->V(AlphaBeta, aobs_space_, bobs_space_);

    {
      const std::string fkey_pA = ParsedOneBodyIntKey::key(aobs_space->id(),cabs_space->id(),std::string("F"));
      RefSCMatrix F = fb_rtime->get(fkey_pA);
      FpA_[Alpha] = F;
    }
    {
      const std::string fkey_Ap = ParsedOneBodyIntKey::key(cabs_space->id(),aobs_space->id(),std::string("F"));
      RefSCMatrix F = fb_rtime->get(fkey_Ap);
      FAp_[Alpha] = F;
    }

    {
      const std::string
          tkey = ParsedTwoBodyFourCenterIntKey::key(skey, skey, skey, cabs_space->id(),
                                                    std::string("ERI"), std::string(""),
                                                    TwoBodyIntLayout::b1b2_k1k2);

      Ref<TwoBodyMOIntsTransform> pppA_tform =
          r12world()->world()->moints_runtime()->runtime_4c()->get(tkey);
      pppA_tform->compute();
      pppA_acc_[AlphaBeta] = pppA_tform->ints_distarray4();
      pppA_acc_[AlphaBeta]->activate();
    }

    Ref<R12Technology::CorrelationFactor> corrfactor = r12world()->r12tech()->corrfactor();
    // only 1 correlation factor can be handled
    MPQC_ASSERT(corrfactor->nfunctions() == 1);
    Ref<TwoBodyIntDescr> tbdescr = corrfactor->tbintdescr(r12world()->integral(),0);
    const std::string descr_key = r12world()->world()->moints_runtime4()->descr_key(tbdescr);
    {
      Ref<OrbitalSpace> aocc_space = r12world()->refwfn()->occ_act_sb(Alpha);
      Ref<OrbitalSpace> avir_space = r12world()->refwfn()->uocc_act_sb(Alpha);
      const std::string
        tkey = ParsedTwoBodyFourCenterIntKey::key(aocc_space->id(), aocc_space->id(),
                                                  avir_space->id(), cabs_space->id(),
                                                  descr_key, TwoBodyIntLayout::b1b2_k1k2);
      Ref<TwoBodyMOIntsTransform> iiaA_tform =
        r12world()->world()->moints_runtime()->runtime_4c()->get(tkey);
      iiaA_tform->compute();
      iiaA_acc_[AlphaBeta] = iiaA_tform->ints_distarray4();
      iiaA_acc_[AlphaBeta]->activate();
    }
  } // (R12) contributions

}

void
CCR12_Info::compute_source_integrals_uhf() {

  MPQC_ASSERT(false);

  // this function should only be called for rohf or uhf
  if (!r12world()->refwfn()->spin_polarized())
    throw ProgrammingError("CCR12_Info::compute_source_integrals_uhf -- reference wfn is not spin-polarized",
                           __FILE__,__LINE__);

  // and CCSD
  if (theory_ != std::string("CCSD") && theory_ != std::string("CCSDT") && theory_ != std::string("CCSDTQ")
   && theory_ != std::string("CCSD-R12") && theory_ != std::string("CCSD(R12)"))
    throw FeatureNotImplemented("Can only import integrals for CCSD, CCSD-R12, CCSDT, and CCSDTQ methods",
                                    __FILE__,__LINE__);

  // Prepare OrbitalSpace that corresponds to the CorrelatedMOOrder assumed by SMITH

  Ref<OrbitalSpace> occs_a = r12world()->refwfn()->occ_sb(Alpha);
  Ref<OrbitalSpace> occs_b = r12world()->refwfn()->occ_sb(Beta);
  Ref<OrbitalSpace> orbs_a = r12world()->refwfn()->orbs_sb(Alpha);
  Ref<OrbitalSpace> orbs_b = r12world()->refwfn()->orbs_sb(Beta);
  const size_t norbs = orbs_a->rank();
  RefDiagSCMatrix occnums_a = orbs_a->evals().clone();
  occnums_a.assign(0.0);
  RefDiagSCMatrix occnums_b = orbs_b->evals().clone();
  occnums_b.assign(0.0);
  const unsigned int nirreps = occnums_a.dim()->blocks()->nblock();
  const std::vector<unsigned int>& occpi_a = occs_a->block_sizes();
  const std::vector<unsigned int>& occpi_b = occs_b->block_sizes();
  for (unsigned int irrep = 0; irrep < nirreps; ++irrep) {
    const unsigned int nocc_a = occpi_a[irrep];
    const unsigned int occ_a_begin = occnums_a.dim()->blocks()->start(irrep);
    const unsigned int occ_a_end = occ_a_begin + nocc_a;
    for (unsigned int o = occ_a_begin; o < occ_a_end; ++o)
      occnums_a( o) = 1.0;
    const unsigned int nocc_b = occpi_b[irrep];
    const unsigned int occ_b_begin = occnums_b.dim()->blocks()->start(irrep);
    const unsigned int occ_b_end = occ_b_begin + nocc_b;
    for (unsigned int o = occ_b_begin; o < occ_b_end; ++o)
      occnums_b( o) = 1.0;
  }

  Ref<OrbitalSpace>
      full_space =
          new OrderedSpinOrbitalSpace<CorrelatedSpinMOOrder> (
                                                              std::string(
                                                                          "pp(corr)"),
                                                              std::string(
                                                                          "Active reference spin-orbitals in correlated order"),
                                                              orbs_a->basis(),
                                                              orbs_a->integral(),
                                                              orbs_a->coefs(),
                                                              orbs_b->coefs(),
                                                              orbs_a->evals(),
                                                              orbs_b->evals(),
                                                              occnums_a,
                                                              occnums_b,
                                                              orbs_a->orbsym(),
                                                              orbs_b->orbsym(),
                                                              CorrelatedSpinMOOrder(
                                                                                    nirreps));
  corr_space_ = full_space;

  // minus frozen core
  {
    typedef MolecularOrbitalMask<double, RefDiagSCMatrix> FCMask;
    FCMask fcmask(2 * nfzc_, corr_space_->evals());
    Ref<OrbitalSpace>
        full_space_minus_fc =
            new MaskedOrbitalSpace(
                                   std::string("pp(corr)"),
                                   std::string(
                                               "Active reference spin-orbitals in correlated order"),
                                   corr_space_, fcmask.mask());
    corr_space_ = full_space_minus_fc;
  }
  // minus frozen virtuals
  {
    typedef MolecularOrbitalMask<double, RefDiagSCMatrix, std::greater<double> >
        FVMask;
    FVMask fvmask(2 * nfzv_, corr_space_->evals());
    Ref<OrbitalSpace>
        full_space_minus_fv =
            new MaskedOrbitalSpace(
                                   std::string("pp(corr)"),
                                   std::string(
                                               "Active reference spin-orbitals in correlated order"),
                                   corr_space_, fvmask.mask());
    corr_space_ = full_space_minus_fv;
  }
  const Ref<OrbitalSpaceRegistry> idxreg = r12world()->world()->moints_runtime()->factory()->orbital_registry();
  idxreg->add(make_keyspace_pair(corr_space_));

  const std::string pspace_a_id = ParsedOrbitalSpaceKey::key(std::string("p(corr)"),Alpha);
  const std::string pspace_a_name("Alpha active reference orbitals in correlated order");
  Ref<OrbitalSpace> pspace_a =
      new OrderedOrbitalSpace<CorrelatedMOOrder> (pspace_a_id, pspace_a_name,
                                                  orbs_a->basis(),
                                                  orbs_a->integral(),
                                                  orbs_a->coefs(),
                                                  orbs_a->evals(), occnums_a,
                                                  orbs_a->orbsym(),
                                                  CorrelatedMOOrder(nirreps));
  // minus frozen core
  {
    typedef MolecularOrbitalMask<double,RefDiagSCMatrix> FCMask;
    FCMask fcmask(nfzc_,pspace_a->evals());
    pspace_a = new MaskedOrbitalSpace(pspace_a_id, pspace_a_name,
                                      pspace_a,
                                      fcmask.mask());
  }
  // minus frozen virtuals
  {
    typedef MolecularOrbitalMask<double,RefDiagSCMatrix,std::greater<double> > FVMask;
    FVMask fvmask(nfzv_,pspace_a->evals());
    pspace_a = new MaskedOrbitalSpace(pspace_a_id, pspace_a_name,
                                      pspace_a,
                                      fvmask.mask());
  }
  // get rid of the blocking in case the Fock matrix is not diagonal
  pspace_a = new NonblockedOrbitalSpace(pspace_a_id,pspace_a_name,pspace_a);
  idxreg->add(make_keyspace_pair(pspace_a));

  const std::string pspace_b_id = ParsedOrbitalSpaceKey::key(std::string("p(corr)"),Beta);
  const std::string pspace_b_name("Beta active reference orbitals in correlated order");
  Ref<OrbitalSpace> pspace_b =
      new OrderedOrbitalSpace<CorrelatedMOOrder> (pspace_b_id, pspace_b_name,
                                                  orbs_b->basis(),
                                                  orbs_b->integral(),
                                                  orbs_b->coefs(),
                                                  orbs_b->evals(), occnums_b,
                                                  orbs_b->orbsym(),
                                                  CorrelatedMOOrder(nirreps));
  // minus frozen core
  {
    typedef MolecularOrbitalMask<double,RefDiagSCMatrix> FCMask;
    FCMask fcmask(nfzc_,pspace_b->evals());
    pspace_b = new MaskedOrbitalSpace(pspace_b_id, pspace_b_name,
                                      pspace_b,
                                      fcmask.mask());
  }
  // minus frozen virtuals
  {
    typedef MolecularOrbitalMask<double,RefDiagSCMatrix,std::greater<double> > FVMask;
    FVMask fvmask(nfzv_,pspace_b->evals());
    pspace_b = new MaskedOrbitalSpace(pspace_b_id, pspace_b_name,
                                      pspace_b,
                                      fvmask.mask());
  }
  // get rid of the blocking in case the Fock matrix is not diagonal
  pspace_b = new NonblockedOrbitalSpace(pspace_b_id,pspace_b_name,pspace_b);
  idxreg->add(make_keyspace_pair(pspace_b));

#if 0
  corr_space_->print_detail();
  pspace_a->print_detail();
  pspace_b->print_detail();
#endif

  Ref<FockBuildRuntime> fb_rtime = r12world()->world()->fockbuild_runtime();
  const std::string fkey_a = ParsedOneBodyIntKey::key(pspace_a->id(),pspace_a->id(),std::string("F"),Alpha);
  F_[Alpha] = fb_rtime->get(fkey_a);
  const std::string fkey_b = ParsedOneBodyIntKey::key(pspace_b->id(),pspace_b->id(),std::string("F"),Beta);
  F_[Beta] = fb_rtime->get(fkey_b);
  const std::string hkey_a = ParsedOneBodyIntKey::key(pspace_a->id(),pspace_a->id(),std::string("H"));
  RefSCMatrix H_a = fb_rtime->get(hkey_a);  F_[Alpha].accumulate(H_a);
  const std::string hkey_b = ParsedOneBodyIntKey::key(pspace_b->id(),pspace_b->id(),std::string("H"));
  RefSCMatrix H_b = fb_rtime->get(hkey_b);  F_[Beta].accumulate(H_b);

#if 0
  F_[Alpha].print("Alpha Fock matrix");
  F_[Beta].print("Beta Fock matrix");
#endif

  const std::string skey = corr_space_->id();
  const std::string
      tkey = ParsedTwoBodyFourCenterIntKey::key(skey, skey, skey, skey,
                                                std::string("ERI"), std::string(""),
                                                TwoBodyIntLayout::b1b2_k1k2);

  Ref<TwoBodyMOIntsTransform> pppp_tform =
      r12world()->world()->moints_runtime()->runtime_4c()->get(tkey);
  pppp_tform->compute();
  pppp_acc_[AlphaBeta] = pppp_tform->ints_distarray4();
  pppp_acc_[AlphaBeta]->activate();

}

#define LOCAL_DEBUG 0

void CCR12_Info::fill_in_f1(){

  Ref<OrbitalSpace> aobs_space = r12world()->refwfn()->orbs_sb(Alpha);
  Ref<OrbitalSpace> bobs_space = r12world()->refwfn()->orbs_sb(Beta);

  // compute map from indices in full spin-orbital space to indices in the respective spin spaces
  vector<long> amap;
  {
    vector<int> smith_to_obs = sc::map(*aobs_space, *corr_space_, false);
    amap.resize(smith_to_obs.size());
    std::copy(smith_to_obs.begin(), smith_to_obs.end(), amap.begin());
  }
  vector<long> bmap;
  {
    vector<int> smith_to_obs = map(*bobs_space, *corr_space_, false);
    bmap.resize(smith_to_obs.size());
    std::copy(smith_to_obs.begin(), smith_to_obs.end(), bmap.begin());
  }

  MTensor<2>::tile_ranges pp(2, MTensor<2>::tile_range(0, this->nab()));
  MTensor<2> f1(this,d_f1.pointer(),pp);
  // assuming RHF!
  const int norbs = aobs_space->rank();

  MTensor<2>::element_ranges pp_erange(2, MTensor<2>::element_range(0, norbs) );
  f1.convert(F_[Alpha], amap, amap, &pp_erange);
  if (this->ref()->spin_polarized())
    f1.convert(F_[Beta], bmap, bmap, &pp_erange);

  if (need_w1() || need_w2()) {

    // get the CABS space
    Ref<OrbitalSpace> acabs_space = r12world()->cabs_space(Alpha);
    Ref<OrbitalSpace> bcabs_space = r12world()->cabs_space(Beta);

    vector<long> acabsmap;
    {
      vector<int> smith_to_cabs( map(*acabs_space, *corr_space_) );
      acabsmap.resize(smith_to_cabs.size());
      std::copy(smith_to_cabs.begin(), smith_to_cabs.end(), acabsmap.begin());
    }
    vector<long> bcabsmap;
    {
      vector<int> smith_to_cabs( map(*bcabs_space, *corr_space_) );
      bcabsmap.resize(smith_to_cabs.size());
      std::copy(smith_to_cabs.begin(), smith_to_cabs.end(), bcabsmap.begin());
    }

    const int ncabs = acabs_space->rank();

    {
      MTensor<2>::element_ranges pA_erange(2);
      pA_erange[0] = MTensor<2>::element_range(0, norbs);
      pA_erange[1] = MTensor<2>::element_range(0, ncabs);
      f1.convert(FpA_[Alpha], amap, acabsmap, &pA_erange);
    }
    {
      MTensor<2>::element_ranges Ap_erange(2);
      Ap_erange[0] = MTensor<2>::element_range(0, ncabs);
      Ap_erange[1] = MTensor<2>::element_range(0, norbs);
      f1.convert(FAp_[Alpha], acabsmap, amap, &Ap_erange);
    }

    if (ref()->spin_polarized()) {
    {
      MTensor<2>::element_ranges pA_erange(2);
      pA_erange[0] = MTensor<2>::element_range(0, norbs);
      pA_erange[1] = MTensor<2>::element_range(0, ncabs);
      f1.convert(FpA_[Beta], bmap, bcabsmap, &pA_erange);
    }
    {
      MTensor<2>::element_ranges Ap_erange(2);
      Ap_erange[0] = MTensor<2>::element_range(0, ncabs);
      Ap_erange[1] = MTensor<2>::element_range(0, norbs);
      f1.convert(FAp_[Beta], bcabsmap, bmap, &Ap_erange);
    }
    }

  }
}

void CCR12_Info::fill_in_v2() {

  MPQC_ASSERT(restricted_);

  // compute map from indices in full spin-orbital space to indices in the respective spin spaces
  vector<long> amap;
  {
    vector<int> smith_to_obs = sc::map(*aobs_space_, *corr_space_, false);
    amap.resize(smith_to_obs.size());
    std::copy(smith_to_obs.begin(), smith_to_obs.end(), amap.begin());
  }
  vector<long> bmap;
  {
    vector<int> smith_to_obs = map(*bobs_space_, *corr_space_, false);
    bmap.resize(smith_to_obs.size());
    std::copy(smith_to_obs.begin(), smith_to_obs.end(), bmap.begin());
  }

  const unsigned int norbs = aobs_space_->rank();

  MTensor<4>::tile_ranges pppp(4, MTensor<4>::tile_range(0, this->nab()));
  MTensor<4> v2(this,d_v2.pointer(),pppp);
  {
    MTensor<4>::element_ranges pppp_erange(4, MTensor<4>::element_range(0, norbs) );
    v2.convert(pppp_acc_[AlphaBeta], 0u, amap, amap, amap, amap, &pppp_erange);
  }

  if (need_w1()) { // (R12)

    // compute map from indices in full spin-orbital space to indices in the respective spin spaces
    vector<long> aAmap;
    {
      vector<int> smith_to_obs = sc::map(*(r12world()->cabs_space(Alpha)), *corr_space_, false);
      aAmap.resize(smith_to_obs.size());
      std::copy(smith_to_obs.begin(), smith_to_obs.end(), aAmap.begin());
    }
    vector<long> bAmap;
    {
      vector<int> smith_to_obs = map(*(r12world()->cabs_space(Beta)), *corr_space_, false);
      bAmap.resize(smith_to_obs.size());
      std::copy(smith_to_obs.begin(), smith_to_obs.end(), bAmap.begin());
    }
    const unsigned int ncabs = r12world()->cabs_space(Alpha)->rank();
    // pppA
    {
      MTensor<4>::element_ranges pppA_erange(3, MTensor<4>::element_range(0, norbs) );
      pppA_erange.push_back( MTensor<4>::element_range(0, ncabs) );
      v2.convert(pppA_acc_[AlphaBeta], 0u, amap, amap, amap, aAmap, &pppA_erange);
    }
    // pApp
    {
      MTensor<4>::element_ranges pApp_erange(4);
      pApp_erange[0] = ( MTensor<4>::element_range(0, norbs) );
      pApp_erange[1] = ( MTensor<4>::element_range(0, ncabs) );
      pApp_erange[2] = ( MTensor<4>::element_range(0, norbs) );
      pApp_erange[3] = ( MTensor<4>::element_range(0, norbs) );
      const bool transpose_23_01 = true;
      v2.convert(pppA_acc_[AlphaBeta], 0u, amap, aAmap, amap, amap, &pApp_erange, transpose_23_01);
    }

  }

}


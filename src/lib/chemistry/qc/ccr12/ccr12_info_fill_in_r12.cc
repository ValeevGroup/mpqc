//
// ccr12_info_fill_in_r12.cc
//
// Copyright (C) 2008 Toru Shiozaki
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


void CCR12_Info::fill_in_iiii() {

  // only RHF is handled right now
  MPQC_ASSERT(restricted_);

  // fill in X, B and P intermediates, as well as geminal T amplitudes
  MPQC_ASSERT(need_w1());

  // intermediates are computed in GGspace(Alpha) space
  // map corr_space_ to it
  // compute map from indices in full spin-orbital space to indices in the respective spin spaces
  vector<long> amap;
  {
    vector<int> intmap = sc::map(*(r12eval()->GGspace(Alpha)), *corr_space_, false);
    amap.resize(intmap.size());
    std::copy(intmap.begin(), intmap.end(), amap.begin());
  }

  // assuming RHF!
  const int nG = r12eval()->GGspace(Alpha)->rank();
  MTensor<4>::tile_ranges iiii(4, MTensor<4>::tile_range(0, this->noab()));
  MTensor<4>::element_ranges iiii_erange(4, MTensor<4>::element_range(0, nG) );
  {
    MTensor<4> X(this,d_xs2.pointer(),iiii);
    X.convert(r12int_eval_->X(AlphaBeta), nG, nG,
              amap, amap, amap, amap, &iiii_erange);
  }

  {
    MTensor<4> B(this,d_bs2.pointer(),iiii);
    B.convert(r12int_eval_->B(AlphaBeta), nG, nG,
              amap, amap, amap, amap, &iiii_erange);
  }

  if (need_gt2()) {
    // (2)R12 methods do not need gt2 tensors (but need all others).
    // In fullopt calculations, we still use it for guess functions.
    MTensor<4> GT2(this,d_gt2.pointer(),iiii);

    // columns of GT2 are computed in product ggspace(Alpha) space
    // map corr_space_ to it
    vector<long> amap_o;
    {
      vector<int> intmap = sc::map(*(r12eval()->occ_act(Alpha)), *corr_space_, false);
      amap_o.resize(intmap.size());
      std::copy(intmap.begin(), intmap.end(), amap_o.begin());
    }
    const int no = r12eval()->occ_act(Alpha)->rank();

    // compute fixed geminal coefficients
    RefSCMatrix gt2 = X_.kit()->matrix(r12int_eval_->dim_oo(AlphaBeta),
                                       r12int_eval_->dim_oo(AlphaBeta));
    firstorder_cusp_coefficients(AlphaBeta,
                                 gt2,
                                 r12int_eval_->occ_act(Alpha),
                                 r12int_eval_->occ_act(Beta),
                                 r12int_eval_->r12world()->r12tech()->corrfactor());
    MTensor<4>::element_ranges iiii_erange(4, MTensor<4>::element_range(0, no) );
    GT2.convert(gt2, no, no,
                amap_o, amap_o, amap_o, amap_o, &iiii_erange);

#if 0
    ExEnv::out0() << "TEST: mp2-r12 energy = " << mp2r12_energy->energy() << endl;
    {
      double C_dot_V_tr = 0.0;

#if 1
      {
        SpinCase2 spin = AlphaAlpha;
        RefSCMatrix C_dot_V = mp2r12_energy->C(spin) * r12int_eval_->V(spin);
        C_dot_V.print("C.V(AA)");
        const int noo = r12int_eval_->dim_oo(spin).n();
        for(int ij=0; ij<noo; ++ij) {
          C_dot_V_tr += C_dot_V(ij,ij);
        }
        C_dot_V_tr *= 2.0;
      }
#endif
#if 1
      {
        SpinCase2 spin = AlphaBeta;
        RefSCMatrix C_dot_V = mp2r12_energy->C(spin) * r12int_eval_->V(spin);
        C_dot_V.print("C.V(AB)");
        const int noo = r12int_eval_->dim_oo(spin).n();
        for(int ij=0; ij<noo; ++ij) {
          C_dot_V_tr += C_dot_V(ij,ij);
        }
      }
#endif

      ExEnv::out0() << "TEST: tr(C.V) = " << C_dot_V_tr << endl;
    }
#endif

  } // GT2

}


void CCR12_Info::fill_in_vr_and_vd() {

  // only RHF is handled so far
  MPQC_ASSERT(restricted_);

  // for V intermediate (denoted as d_vr2 in smith; its transpose is d_vd2)
  MPQC_ASSERT(need_w1());

  // intermediates are computed using i indices from GGspace(Alpha) space
  // map corr_space_ to it
  vector<long> aimap;
  {
    vector<int> intmap = sc::map(*(r12eval()->GGspace(Alpha)), *corr_space_, false);
    aimap.resize(intmap.size());
    std::copy(intmap.begin(), intmap.end(), aimap.begin());
  }
  // g indices are in aobs_space_, map to it also
  vector<long> agmap;
  {
    vector<int> intmap = sc::map(*aobs_space_, *corr_space_, false);
    agmap.resize(intmap.size());
    std::copy(intmap.begin(), intmap.end(), agmap.begin());
  }


  // assuming RHF!
  const int nG = r12eval()->GGspace(Alpha)->rank();
  const int norbs_act = aobs_space_->rank();

  // d_vr2 = ggii
  {
    MTensor<4>::tile_ranges ggii(4);
    ggii[0] = MTensor<4>::tile_range(0, this->nab());
    ggii[1] = MTensor<4>::tile_range(0, this->nab());
    ggii[2] = MTensor<4>::tile_range(0, this->noab());
    ggii[3] = MTensor<4>::tile_range(0, this->noab());
    MTensor<4> V(this,d_vr2.pointer(),ggii);
    MTensor<4>::element_ranges ggii_erange(4);
    ggii_erange[0] = MTensor<4>::element_range(0, norbs_act);
    ggii_erange[1] = MTensor<4>::element_range(0, norbs_act);
    ggii_erange[2] = MTensor<4>::element_range(0, nG);
    ggii_erange[3] = MTensor<4>::element_range(0, nG);
    V.convert<RefSCMatrix>(Vgg_[AlphaBeta].t(), norbs_act, nG,
              agmap, agmap, aimap, aimap, &ggii_erange);
  }

  // d_vd2 = iigg
  {
    MTensor<4>::tile_ranges iigg(4);
    iigg[0] = MTensor<4>::tile_range(0, this->noab());
    iigg[1] = MTensor<4>::tile_range(0, this->noab());
    iigg[2] = MTensor<4>::tile_range(0, this->nab());
    iigg[3] = MTensor<4>::tile_range(0, this->nab());
    MTensor<4> V(this,d_vd2.pointer(),iigg);
    MTensor<4>::element_ranges iigg_erange(4);
    iigg_erange[0] = MTensor<4>::element_range(0, nG);
    iigg_erange[1] = MTensor<4>::element_range(0, nG);
    iigg_erange[2] = MTensor<4>::element_range(0, norbs_act);
    iigg_erange[3] = MTensor<4>::element_range(0, norbs_act);
    V.convert<RefSCMatrix>(Vgg_[AlphaBeta], nG, norbs_act,
              aimap, aimap, agmap, agmap, &iigg_erange);
  }

}

void CCR12_Info::fill_in_vd2_gen(bool need_cabs, bool need_xx) {

  // only RHF is handled so far
  MPQC_ASSERT(restricted_);

  // intermediates are computed using i indices from GGspace(Alpha) space
  // map corr_space_ to it
  vector<long> aimap;
  {
    vector<int> intmap = sc::map(*(r12eval()->GGspace(Alpha)), *corr_space_, false);
    aimap.resize(intmap.size());
    std::copy(intmap.begin(), intmap.end(), aimap.begin());
  }
  // g indices are in aobs_space_, map to it also
  vector<long> agmap;
  {
    vector<int> intmap = sc::map(*aobs_space_, *corr_space_, false);
    agmap.resize(intmap.size());
    std::copy(intmap.begin(), intmap.end(), agmap.begin());
  }

  // assuming RHF!
  const int nG = r12eval()->GGspace(Alpha)->rank();
  const int norbs_act = aobs_space_->rank();

  // d_vd2 = iigg
  {
    MTensor<4>::tile_ranges iigg(4);
    const int nobs_b = this->noab() + this->nvab();
    iigg[0] = MTensor<4>::tile_range(0, nobs_b);
    iigg[1] = MTensor<4>::tile_range(0, nobs_b);
    iigg[2] = MTensor<4>::tile_range(0, this->nab());
    iigg[3] = MTensor<4>::tile_range(0, this->nab());
    MTensor<4> V(this,d_vd2_gen.pointer(),iigg);
    MTensor<4>::element_ranges iigg_erange(4);
    iigg_erange[0] = MTensor<4>::element_range(0, nG);
    iigg_erange[1] = MTensor<4>::element_range(0, nG);
    iigg_erange[2] = MTensor<4>::element_range(0, norbs_act);
    iigg_erange[3] = MTensor<4>::element_range(0, norbs_act);
    V.convert<RefSCMatrix>(Vgg_[AlphaBeta], nG, norbs_act,
              aimap, aimap, agmap, agmap, &iigg_erange);
  }

}

void CCR12_Info::fill_in_fr_and_fd() {


  // only RHF is handled so far
  MPQC_ASSERT(restricted_);
  MPQC_ASSERT(need_w1());

  // intermediates are computed using i indices from occ_act_sb(Alpha) space
  // map corr_space_ to it
  vector<long> aimap;
  {
    vector<int> intmap = sc::map(*(r12eval()->occ_act(Alpha)), *corr_space_, false);
    aimap.resize(intmap.size());
    std::copy(intmap.begin(), intmap.end(), aimap.begin());
  }
  // a indices are in vir_act_sb(Alpha), map to it also
  vector<long> aamap;
  {
    vector<int> intmap = sc::map(*(r12eval()->vir_act(Alpha)), *corr_space_, false);
    aamap.resize(intmap.size());
    std::copy(intmap.begin(), intmap.end(), aamap.begin());
  }
  // A indices are in cabs_space(Alpha), map to it also
  vector<long> aAmap;
  {
    vector<int> intmap = sc::map(*(r12world()->cabs_space(Alpha)), *corr_space_, false);
    aAmap.resize(intmap.size());
    std::copy(intmap.begin(), intmap.end(), aAmap.begin());
  }

  // assuming RHF!
  const int nocc_act = r12eval()->occ_act(Alpha)->rank();
  const int nvir_act = r12eval()->vir_act(Alpha)->rank();
  const int ncabs = r12world()->cabs_space(Alpha)->rank();

  // d_fr2 = aAii
  {
    MTensor<4>::tile_ranges aAii(4);
    aAii[0] = MTensor<4>::tile_range(0, this->nab());
    aAii[1] = MTensor<4>::tile_range(0, this->nab());
    aAii[2] = MTensor<4>::tile_range(0, this->noab());
    aAii[3] = MTensor<4>::tile_range(0, this->noab());
    MTensor<4> F(this,d_fr2.pointer(),aAii);
    MTensor<4>::element_ranges aAii_erange(4);
    aAii_erange[0] = MTensor<4>::element_range(0, nvir_act);
    aAii_erange[1] = MTensor<4>::element_range(0, ncabs);
    aAii_erange[2] = MTensor<4>::element_range(0, nocc_act);
    aAii_erange[3] = MTensor<4>::element_range(0, nocc_act);
    const bool transpose_23_01 = true;
    Ref<TwoBodyIntDescr> tbintdescr = r12world()->r12tech()->corrfactor()->tbintdescr(r12world()->integral(), 0);
    MPQC_ASSERT(r12world()->r12tech()->corrfactor()->nfunctions() == 1);
    F.convert(iiaA_acc_[AlphaBeta], tbintdescr->intset( r12world()->r12tech()->corrfactor()->tbint_type_f12() ),
              aamap, aAmap, aimap, aimap, &aAii_erange, transpose_23_01);
  }

  // d_fd2 = iiaA
  {
    MTensor<4>::tile_ranges iiaA(4);
    iiaA[0] = MTensor<4>::tile_range(0, this->noab());
    iiaA[1] = MTensor<4>::tile_range(0, this->noab());
    iiaA[2] = MTensor<4>::tile_range(0, this->nab());
    iiaA[3] = MTensor<4>::tile_range(0, this->nab());
    MTensor<4> F(this,d_fd2.pointer(),iiaA);
    MTensor<4>::element_ranges iiaA_erange(4);
    iiaA_erange[0] = MTensor<4>::element_range(0, nocc_act);
    iiaA_erange[1] = MTensor<4>::element_range(0, nocc_act);
    iiaA_erange[2] = MTensor<4>::element_range(0, nvir_act);
    iiaA_erange[3] = MTensor<4>::element_range(0, ncabs);
    Ref<TwoBodyIntDescr> tbintdescr = r12world()->r12tech()->corrfactor()->tbintdescr(r12world()->integral(), 0);
    MPQC_ASSERT(r12world()->r12tech()->corrfactor()->nfunctions() == 1);
    F.convert(iiaA_acc_[AlphaBeta], tbintdescr->intset( r12world()->r12tech()->corrfactor()->tbint_type_f12() ),
              aimap, aimap, aamap, aAmap, &iiaA_erange);
  }

  //ExEnv::out0() << "Norms of fd2 and fr2 = " << RMS(*d_fd2) << " " << RMS(*d_fr2) << endl;

}



//
// sr_r12intermediates_util.h
//
// Copyright (C) 2013 Edward Valeev
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

// this should not be included directly, will be included by sr_r12intermediates.h

#ifndef _mpqc_src_lib_chemistry_qc_mbptr12_srr12intermediatesutil_h
#define _mpqc_src_lib_chemistry_qc_mbptr12_srr12intermediatesutil_h

#include <algorithm>
#include <numeric>
#include <TiledArray/bitset.h>

namespace sc {

  template <typename T>
  std::shared_ptr<typename SingleReference_R12Intermediates<T>::TArray22>
  SingleReference_R12Intermediates<T>::ab_O_cd(const std::string& key,
                                               const int te_type) {

    typedef typename TArray22::value_type value_type;

    Ref<TwoBodyMOIntsTransform> tform = r12world_->world()->moints_runtime4()->get(key);
    tform->compute();
    Ref<DistArray4> darray4 = tform->ints_distarray4();
    darray4->activate();

    // this only works if all nodes can read integrals
    std::vector<int> tasks_with_access;
    const int ntasks_with_access = darray4->tasks_with_access(tasks_with_access);
    MPQC_ASSERT(ntasks_with_access == world_.nproc());
    MPQC_ASSERT(ntasks_with_access == r12world_->world()->msg()->n());

    const size_t n1 = darray4->ni();
    const size_t n2 = darray4->nj();
    const size_t n3 = darray4->nx();
    const size_t n4 = darray4->ny();

    // make tiled ranges for occupied indices
    // N.B. using little trickery using partial sum to create tilesize=1 range
    std::vector<size_t> i1_hashmarks(n1+1, 1);
    i1_hashmarks[0] = 0;
    std::partial_sum(i1_hashmarks.begin(), i1_hashmarks.end(), i1_hashmarks.begin());
    std::vector<size_t> i2_hashmarks(n2+1, 1);
    i2_hashmarks[0] = 0;
    std::partial_sum(i2_hashmarks.begin(), i2_hashmarks.end(), i2_hashmarks.begin());

    std::vector<TA::TiledRange1> hashmarks;
    hashmarks.push_back(TiledArray::TiledRange1(i1_hashmarks.begin(), i1_hashmarks.end()));
    hashmarks.push_back(TiledArray::TiledRange1(i2_hashmarks.begin(), i2_hashmarks.end()));

    TiledArray::TiledRange i1i2_trange(hashmarks.begin(), hashmarks.end());

    // make an empty tensor now
    std::shared_ptr<TArray22> i1i2_g_p1p2(new TArray22(world_, i1i2_trange));

    // make ordinary ranges needed to make p1p2 tiles
    std::vector<size_t> p1p2_start(2, 0);
    std::vector<size_t> p1p2_finish(2);  p1p2_finish[0] = n3; p1p2_finish[1] = n4;
    TA::Range p1p2_range(p1p2_start, p1p2_finish);

    { // load all local tiles

      std::vector<double> p1p2_block(n3 * n4);

      std::vector<size_t> i1i2(2, 0);
      for(i1i2[0] = 0; i1i2[0]<n1; ++i1i2[0]) {
        for(i1i2[1] = 0; i1i2[1]<n2; ++i1i2[1]) {
          if (i1i2_g_p1p2->is_local(i1i2)) {

            darray4->retrieve_pair_block(i1i2[0], i1i2[1], te_type, &p1p2_block[0]);

            madness::Future < value_type >
              tile((value_type(i1i2_g_p1p2->trange().make_tile_range(i1i2),
                               (typename value_type::value_type(p1p2_range,
                                                                 &p1p2_block[0])
                                         )
                                        )
                  ));

            // Insert the tile into the array
            i1i2_g_p1p2->set(i1i2, tile);

            darray4->release_pair_block(i1i2[0], i1i2[1], te_type);

          }
        }
      }

    }

    return i1i2_g_p1p2;
  }

  template <typename T>
  typename SingleReference_R12Intermediates<T>::TArray22&
  SingleReference_R12Intermediates<T>::_(const std::string& key) {
    auto v = tarray22_registry_.find(key);
    if (tarray22_registry_.find(key) == tarray22_registry_.end()) {
      ParsedTwoBodyFourCenterIntKey pkey(key);

      /// map operator to the index within the operator set
      std::string operset_label;
      const unsigned int oper_idx = 0;
      if (pkey.oper() == "g") {
        operset_label = "ERI";
      }
      else if (pkey.oper() == "gr") {
        operset_label = "R12_m1_G12[0]";
      }
      else if (pkey.oper() == "r") {
        operset_label = "R12_0_G12[0]";
      }
      else if (pkey.oper() == "r2") {
        operset_label = "R12_0_G12[0,0]";
      }
      else if (pkey.oper() == "rTr") {
        operset_label = "G12_T1_G12[0,0]";
      }
      else
        throw ProgrammingError("SingleReference_R12Intermediates<T>::_ : invalid operator key",__FILE__,__LINE__);

      // canonicalize indices, may compute spaces if needed
      auto bra1 = to_space(pkey.bra1());
      auto bra2 = to_space(pkey.bra2());
      auto ket1 = to_space(pkey.ket1());
      auto ket2 = to_space(pkey.ket2());

      std::string tform_key = ParsedTwoBodyFourCenterIntKey::key(bra1, bra2, ket1, ket2,
                                                                 operset_label, "");

      tarray22_registry_[key] = ab_O_cd(tform_key, oper_idx);
      v = tarray22_registry_.find(key);
    }

    return *(v->second);
  }

  template <typename T>
  typename SingleReference_R12Intermediates<T>::TArray4d&
  SingleReference_R12Intermediates<T>::ijxy(const std::string& key) {

    auto v = tarray4_registry_.find(key);
    if (tarray4_registry_.find(key) == tarray4_registry_.end()) {

    ParsedTwoBodyFourCenterIntKey pkey(key);

    /// map operator to the index within the operator set
    std::string operset_label;
    bool rdm2 = false;
    bool t2 = false;
    bool l2 = false;
    const unsigned int oper_idx = 0;
    if (pkey.oper() == "g") {
      operset_label = "ERI";
    }
    else if (pkey.oper() == "gr") {
      operset_label = "R12_m1_G12[0]";
    }
    else if (pkey.oper() == "r") {
      operset_label = "R12_0_G12[0]";
    }
    else if (pkey.oper() == "r2") {
      operset_label = "R12_0_G12[0,0]";
    }
    else if (pkey.oper() == "rTr") {
      operset_label = "G12_T1_G12[0,0]";
    }
    else if (pkey.oper() == "gamma") {
      rdm2 = true;
    }
    else if (pkey.oper() == "T2") {
      t2 = true;
    }
    else if (pkey.oper() == "L2") {
      l2 = true;
    }
    else
      throw ProgrammingError("SingleReference_R12Intermediates<T>::ijxy : invalid operator key",__FILE__,__LINE__);

    // canonicalize indices, may compute spaces if needed
    auto bra1 = to_space(pkey.bra1());
    auto bra2 = to_space(pkey.bra2());
    auto ket1 = to_space(pkey.ket1());
    auto ket2 = to_space(pkey.ket2());
    auto oreg = r12world_->world()->tfactory()->orbital_registry();

    Ref<DistArray4> darray4;
    if (rdm2) { // rdm2
      if (rdm2_.null())
        throw ProgrammingError("SingleReference_R12Intermediates<T>::ijxy: asked for rdm2, but it had not been given");
      // if requested spaces don't match exactly, make a new DistArray4
      auto bra1_space = oreg->value(bra1);
      auto bra2_space = oreg->value(bra2);
      auto ket1_space = oreg->value(ket1);
      auto ket2_space = oreg->value(ket2);
      const auto& rdm_space = rdm2_->orbs();

      darray4 = rdm2_->da4();
      if (not (*bra1_space == *rdm_space &&
          *bra2_space == *rdm_space &&
          *ket1_space == *rdm_space &&
          *ket2_space == *rdm_space )) {
        DistArray4Dimensions dims(1, bra1_space->rank(), bra2_space->rank(), ket1_space->rank(), ket2_space->rank());
        Ref<DistArray4> darray4_mapped = darray4->clone(dims);
        std::vector<int> bra1_map = sc::map(*rdm_space, *bra1_space, true);
        std::vector<int> bra2_map = sc::map(*rdm_space, *bra2_space, true);
        std::vector<int> ket1_map = sc::map(*rdm_space, *ket1_space, true);
        std::vector<int> ket2_map = sc::map(*rdm_space, *ket2_space, true);

        std::vector<double> buf(rdm_space->rank() * rdm_space->rank());
        std::vector<double> map_buf(ket1_space->rank() * ket2_space->rank());

        darray4->activate();
        darray4_mapped->activate();

        for(int b1=0; b1<bra1_space->rank(); ++b1) {
          const int bb1 = bra1_map[b1];
          for(int b2=0; b2<bra2_space->rank(); ++b2) {
            const int bb2 = bra2_map[b2];

            darray4->retrieve_pair_block(bb1, bb2, 0, &buf[0]);

            const int nk1 = ket1_space->rank();
            const int nk2 = ket2_space->rank();
            for(int k1=0, k12=0; k1<nk1; ++k1) {
              const int kk1 = ket1_map[k1];
              const int kk1_offset = kk1 * rdm_space->rank();
              for(int k2=0; k2<nk2; ++k2, ++k12) {
                const int kk2 = ket2_map[k2];
                map_buf[k12] = buf[kk1_offset + kk2];
              }
            }

            darray4_mapped->store_pair_block(b1, b2, 0, &map_buf[0]);
            darray4->release_pair_block(bb1, bb2, 0);
          }
        }
        if (darray4->data_persistent()) darray4->deactivate();
        if (darray4_mapped->data_persistent()) darray4_mapped->deactivate();

        darray4 = darray4_mapped;
      }
    }
    else if (t2) {
      if (t2_[AlphaBeta].null())
        throw ProgrammingError("SingleReference_R12Intermediates<T>::ijxy: asked for T2, but it had not been given");

//      auto oreg = r12world_->world()->tfactory()->orbital_registry();
//      MPQC_ASSERT(oreg->value(bra1) == r12world_->refwfn()->occ_act(Alpha));
//      MPQC_ASSERT(oreg->value(bra2) == r12world_->refwfn()->occ_act(Beta));
//      MPQC_ASSERT(oreg->value(ket1) == r12world_->refwfn()->uocc_act(Alpha));
//      MPQC_ASSERT(oreg->value(ket2) == r12world_->refwfn()->uocc_act(Beta));
      darray4 = t2_[AlphaBeta];
    }
    else if (l2) {
      if (l2_[AlphaBeta].null())
        throw ProgrammingError("SingleReference_R12Intermediates<T>::ijxy: asked for L2, but it had not been given");

      darray4 = l2_[AlphaBeta];
    }
    else { // not rdm2 nor t2 nor l2
      std::string tform_key = ParsedTwoBodyFourCenterIntKey::key(bra1, bra2, ket1, ket2,
                                                                 operset_label, "");

      Ref<TwoBodyMOIntsTransform> tform = r12world_->world()->moints_runtime4()->get(tform_key);
      tform->compute();
      darray4 = tform->ints_distarray4();
    }

    darray4->activate();

    const size_t n1 = darray4->ni();
    const size_t n2 = darray4->nj();
    const size_t n3 = darray4->nx();
    const size_t n4 = darray4->ny();

    // this only works if all nodes can read integrals
    std::vector<int> tasks_with_access;
    const int ntasks_with_access = darray4->tasks_with_access(tasks_with_access);
    MPQC_ASSERT(ntasks_with_access == world_.nproc());
    MPQC_ASSERT(ntasks_with_access == r12world_->world()->msg()->n());

    // make tiled ranges
    std::vector<size_t> i_hashmarks = space_hashmarks(bra1);
    std::vector<size_t> j_hashmarks = space_hashmarks(bra2);
    std::vector<size_t> x_hashmarks = space_hashmarks(ket1);
    std::vector<size_t> y_hashmarks = space_hashmarks(ket2);

    std::vector<TA::TiledRange1> hashmarks;
    hashmarks.push_back(TiledArray::TiledRange1(i_hashmarks.begin(), i_hashmarks.end()));
    hashmarks.push_back(TiledArray::TiledRange1(j_hashmarks.begin(), j_hashmarks.end()));
    hashmarks.push_back(TiledArray::TiledRange1(x_hashmarks.begin(), x_hashmarks.end()));
    hashmarks.push_back(TiledArray::TiledRange1(y_hashmarks.begin(), y_hashmarks.end()));
    TiledArray::TiledRange ijxy_trange(hashmarks.begin(), hashmarks.end());

    std::shared_ptr<TArray4d> result(new TArray4d(world_, ijxy_trange) );
    // construct local tiles
    for(auto t=ijxy_trange.tiles().begin();
        t!=ijxy_trange.tiles().end();
        ++t)
      if (result->is_local(*t))
      {
        std::array<std::size_t, 4> index; std::copy(t->begin(), t->end(), index.begin());
        madness::Future < typename TArray4d::value_type >
          tile((DA4_Tile<T>(result.get(), index, darray4, oper_idx)
              ));

        // Insert the tile into the array
        result->set(*t, tile);
      }

      tarray4_registry_[key] = result;
      v = tarray4_registry_.find(key);
    }

    return *(v->second);
  }

  template <typename T>
  typename SingleReference_R12Intermediates<T>::TArray22d&
  SingleReference_R12Intermediates<T>::ij_xy(const std::string& key) {

    auto v = tarray22d_registry_.find(key);
    if (tarray22d_registry_.find(key) == tarray22d_registry_.end()) {

    ParsedTwoBodyFourCenterIntKey pkey(key);

    /// map operator to the index within the operator set
    std::string operset_label;
    const unsigned int oper_idx = 0;
    if (pkey.oper() == "g") {
      operset_label = "ERI";
    }
    else if (pkey.oper() == "gr") {
      operset_label = "R12_m1_G12[0]";
    }
    else if (pkey.oper() == "r") {
      operset_label = "R12_0_G12[0]";
    }
    else if (pkey.oper() == "r2") {
      operset_label = "R12_0_G12[0,0]";
    }
    else if (pkey.oper() == "rTr") {
      operset_label = "G12_T1_G12[0,0]";
    }
    else
      throw ProgrammingError("SingleReference_R12Intermediates<T>::ij_xy : invalid operator key",__FILE__,__LINE__);

    // canonicalize indices, may compute spaces if needed
    auto bra1 = to_space(pkey.bra1());
    auto bra2 = to_space(pkey.bra2());
    auto ket1 = to_space(pkey.ket1());
    auto ket2 = to_space(pkey.ket2());

    std::string tform_key = ParsedTwoBodyFourCenterIntKey::key(bra1, bra2, ket1, ket2,
                                                               operset_label, "");

    Ref<TwoBodyMOIntsTransform> tform = r12world_->world()->moints_runtime4()->get(tform_key);
    tform->compute();
    Ref<DistArray4> darray4 = tform->ints_distarray4();
    darray4->activate();

    const size_t n1 = darray4->ni();
    const size_t n2 = darray4->nj();
    const size_t n3 = darray4->nx();
    const size_t n4 = darray4->ny();

    // this only works if all nodes can read integrals
    std::vector<int> tasks_with_access;
    const int ntasks_with_access = darray4->tasks_with_access(tasks_with_access);
    MPQC_ASSERT(ntasks_with_access == world_.nproc());
    MPQC_ASSERT(ntasks_with_access == r12world_->world()->msg()->n());

    // make tiled ranges
    // N.B. using little trickery using partial sum to create tilesize=1 range
    std::vector<size_t> i_hashmarks(n1+1, 1);
    i_hashmarks[0] = 0;
    std::partial_sum(i_hashmarks.begin(), i_hashmarks.end(), i_hashmarks.begin());
    std::vector<size_t> j_hashmarks(n2+1, 1);
    j_hashmarks[0] = 0;
    std::partial_sum(j_hashmarks.begin(), j_hashmarks.end(), j_hashmarks.begin());

    std::vector<TA::TiledRange1> hashmarks;
    hashmarks.push_back(TiledArray::TiledRange1(i_hashmarks.begin(), i_hashmarks.end()));
    hashmarks.push_back(TiledArray::TiledRange1(j_hashmarks.begin(), j_hashmarks.end()));
    TiledArray::TiledRange ij_trange(hashmarks.begin(), hashmarks.end());

    std::shared_ptr<TArray22d> result(new TArray22d(world_, ij_trange) );
    // construct local tiles
    for(auto t=ij_trange.tiles().begin();
        t!=ij_trange.tiles().end();
        ++t)
      if (result->is_local(*t))
      {
        std::array<std::size_t, 2> index; std::copy(t->begin(), t->end(), index.begin());
        madness::Future < typename TArray22d::value_type >
          tile((DA4_Tile34<T>(result.get(), index, darray4, oper_idx)
              ));

        // Insert the tile into the array
        result->set(*t, tile);
      }

      tarray22d_registry_[key] = result;
      v = tarray22d_registry_.find(key);
    }

    return *(v->second);
  }

  template <typename T>
  typename SingleReference_R12Intermediates<T>::TArray2&
  SingleReference_R12Intermediates<T>::xy(const std::string& key) {

    auto v = tarray2_registry_.find(key);
    if (v == tarray2_registry_.end()) {

    ParsedOneBodyIntKey pkey(key);

    // canonicalize indices, may compute spaces if needed
    auto bra_id = to_space(pkey.bra());
    auto ket_id = to_space(pkey.ket());

    auto oreg = r12world_->world()->tfactory()->orbital_registry();
    auto freg = r12world_->world()->fockbuild_runtime();
    auto bra = oreg->value(bra_id);
    auto ket = oreg->value(ket_id);
    auto nbra = bra->rank();
    auto nket = ket->rank();

    // grab the matrix
    RefSCMatrix operator_matrix;
    if (pkey.oper() == "J")
      operator_matrix = freg->get(ParsedOneBodyIntKey::key(bra_id, ket_id, "J"));
    else if (pkey.oper() == "K")
      operator_matrix = freg->get(ParsedOneBodyIntKey::key(bra_id, ket_id, "K"));
    else if (pkey.oper() == "F")
      operator_matrix = freg->get(ParsedOneBodyIntKey::key(bra_id, ket_id, "F"));
    else if (pkey.oper() == "h")
      operator_matrix = freg->get(ParsedOneBodyIntKey::key(bra_id, ket_id, "H"));
    else if (pkey.oper() == "hJ")
      operator_matrix = freg->get(ParsedOneBodyIntKey::key(bra_id, ket_id, "H"))
                      + freg->get(ParsedOneBodyIntKey::key(bra_id, ket_id, "J"));
    else if (pkey.oper().find("mu_") == 0 || pkey.oper().find("q_") == 0 ||
             pkey.oper().find("dphi_") == 0 || pkey.oper().find("ddphi_") == 0) {
      operator_matrix = freg->get(ParsedOneBodyIntKey::key(bra_id, ket_id, pkey.oper()));
    }
    else if (pkey.oper() == "gamma")
      operator_matrix = this->rdm1(bra, ket);
    else if (pkey.oper() == "I") {
      operator_matrix = bra->coefs()->kit()->matrix(bra->coefs().coldim(), ket->coefs().coldim());
      if (bra == ket) {
        operator_matrix.assign(0.0);
        for(int i=0; i<nbra; ++i)
          operator_matrix.set_element(i,i,1.0);
      }
      else {
        operator_matrix = sc::overlap(*bra,*ket,operator_matrix->kit());
      }
    }
    else if (pkey.oper() == "T1") {
      if (ket_id == "a") { // if given t1_ explicitly (CC), make sure its size matches
        if (t1_) {
          if (t1_.ncol() != oreg->value(ket_id)->rank())
            throw ProgrammingError("SingleReference_R12Intermediates::xy() -- T1.ncol() != nvir_act",
                                   __FILE__, __LINE__);
          operator_matrix = t1_;
        }
        else if (t1_cabs_) { // if t1_cabs given, this means extract the occ_act x vir_act block
          if (t1_cabs_.ncol() == oreg->value("a'")->rank()) // doesn't make sense if T1 only include CABS
            throw ProgrammingError("SingleReference_R12Intermediates::xy() -- asked for <i|T1|a> but T1_cabs does not include conv. virtuals",
                                   __FILE__, __LINE__);
          Ref<OrbitalSpace> vir_act = oreg->value("a");
          Ref<OrbitalSpace> allvir = oreg->value("A'"); // should exist since t1_cabs_ spans more than just CABS
          MPQC_ASSERT(t1_cabs_.ncol() == allvir->rank()); // this should be the only other possibility
          RefSCMatrix t1_subblock = t1_cabs_.kit()->matrix(t1_cabs_.rowdim(),
                                                           vir_act->coefs()->coldim());
          std::vector<int> vir_to_allvir = sc::map(*allvir, *vir_act, false);
          const int nocc_act = t1_subblock.nrow();
          const int nvir_act = vir_act->rank();
          for(int i=0; i<nocc_act; ++i) {
            for(int a=0; a<nvir_act; ++a) {
              const int Ap = vir_to_allvir[a];
              t1_subblock.set_element(i, a, t1_cabs_.get_element(i, Ap));
            }
          }
          operator_matrix = t1_subblock;
        }
      }
      else {
        MPQC_ASSERT(ket_id == "a'" || ket_id == "A'");
        MPQC_ASSERT(t1_cabs_);
        MPQC_ASSERT(oreg->value(ket_id)->rank() == t1_cabs_.ncol());
        operator_matrix = t1_cabs_;
      }
    }
    else if (pkey.oper() == "L1") {
      MPQC_ASSERT(l1_);
      MPQC_ASSERT(oreg->value(bra_id)->rank() == l1_.nrow());
      MPQC_ASSERT(oreg->value(ket_id)->rank() == l1_.ncol());
      operator_matrix = l1_;
    }
    else {
      std::ostringstream oss;
      oss << "SingleReference_R12Intermediates<T>::xy -- operator " << pkey.oper() << " is not recognized";
      throw ProgrammingError(oss.str().c_str(),
                              __FILE__, __LINE__);
    }

    // make tiled ranges
    std::vector<size_t> bra_hashmarks = space_hashmarks(bra_id);
    std::vector<size_t> ket_hashmarks = space_hashmarks(ket_id);

    std::vector<TA::TiledRange1> hashmarks;
    hashmarks.push_back(TiledArray::TiledRange1(bra_hashmarks.begin(), bra_hashmarks.end()));
    hashmarks.push_back(TiledArray::TiledRange1(ket_hashmarks.begin(), ket_hashmarks.end()));
    TiledArray::TiledRange braket_trange(hashmarks.begin(), hashmarks.end());

    std::shared_ptr<TArray2> result(new TArray2(world_, braket_trange) );
    // construct local tiles
    for(auto t=result->get_pmap()->begin();
        t!=result->get_pmap()->end();
        ++t)
      {

        auto tile_range = result->trange().make_tile_range(*t);
        std::vector<double> ptr_data(tile_range.volume());

        for(size_t r=tile_range.lobound()[0], rc=0;
            r != tile_range.upbound()[0];
            ++r) {
          for(size_t c=tile_range.lobound()[1];
              c != tile_range.upbound()[1];
              ++c, ++rc) {
            ptr_data[rc] = operator_matrix->get_element(r,c);
          }
        }

        madness::Future < typename TArray2::value_type >
          tile(typename TArray2::value_type(tile_range, &ptr_data[0]));

        // Insert the tile into the array
        result->set(*t, tile);
      }

    tarray2_registry_[key] = result;
    v = tarray2_registry_.find(key);
    }

    return *(v->second);
  }

  template <typename T>
  TA::expressions::TsrExpr<const typename SingleReference_R12Intermediates<T>::TArray4d, true>
  SingleReference_R12Intermediates<T>::_4(const std::string& key) {
    const TArray4d& tarray4 = ijxy(key);
    ParsedTwoBodyFourCenterIntKey pkey(key);

    std::array<std::string, 4> ind = {{pkey.bra1(), pkey.bra2(), pkey.ket1(), pkey.ket2()}};
    for(auto v=ind.begin(); v!=ind.end(); ++v) {
      if (ParsedTransformedOrbitalSpaceKey::valid_key(*v) ) {
        ParsedTransformedOrbitalSpaceKey ptkey(*v);
        *v = ptkey.label();
      }
    }

    const std::string annotation = ind[0] + "," + ind[1] + "," + ind[2] + "," + ind[3];
    return tarray4(annotation);
  }

  template <typename T>
  TA::expressions::TsrExpr<const typename SingleReference_R12Intermediates<T>::TArray2, true>
  SingleReference_R12Intermediates<T>::_2(const std::string& key) {
    const TArray2& tarray2 = xy(key);
    ParsedOneBodyIntKey pkey(key);

    std::array<std::string, 2> ind = {{pkey.bra(), pkey.ket()}};
    for(auto v=ind.begin(); v!=ind.end(); ++v) {
      if (ParsedTransformedOrbitalSpaceKey::valid_key(*v) ) {
        ParsedTransformedOrbitalSpaceKey ptkey(*v);
        *v = ptkey.label();
      }
    }

    const std::string annotation = ind[0] + "," + ind[1];
    return tarray2(annotation);
  }

  template <typename T>
  expressions::TGeminalGenerator<T> SingleReference_R12Intermediates<T>::tg_s0_gen(0);
  template <typename T>
  expressions::TGeminalGenerator<T> SingleReference_R12Intermediates<T>::tg_s1_gen(1);

  template <typename T>
  TA::expressions::TsrExpr<const typename SingleReference_R12Intermediates<T>::TArray4Tg, true>
  SingleReference_R12Intermediates<T>::_Tg(const std::string& key) {

    ParsedTwoBodyFourCenterIntKey pkey(key);

    auto v = tarray4tg_registry_.find(key);
    if (v == tarray4tg_registry_.end()) {

    if (pkey.oper() != "Tg") {
      std::ostringstream oss;
      oss << "SingleReference_R12Intermediates<T>::_Tg : " << pkey.oper() << " is an invalid operator key";
      throw ProgrammingError(oss.str().c_str(), __FILE__, __LINE__);
    }

    // canonicalize indices, may compute spaces if needed
    auto bra1 = to_space(pkey.bra1());
    auto bra2 = to_space(pkey.bra2());
    auto ket1 = to_space(pkey.ket1());
    auto ket2 = to_space(pkey.ket2());

    if (not (bra1 == bra2 && bra1 == ket1 && bra1 == ket2)) {
      throw ProgrammingError("SingleReference_R12Intermediates<T>::_Tg : all spaces must be the same", __FILE__, __LINE__);
    }

    // make tiled ranges
    std::vector<size_t> i_hashmarks = space_hashmarks(bra1);
    std::vector<size_t> j_hashmarks = space_hashmarks(bra2);
    std::vector<size_t> x_hashmarks = space_hashmarks(ket1);
    std::vector<size_t> y_hashmarks = space_hashmarks(ket2);

    std::vector<TA::TiledRange1> hashmarks;
    hashmarks.push_back(TiledArray::TiledRange1(i_hashmarks.begin(), i_hashmarks.end()));
    hashmarks.push_back(TiledArray::TiledRange1(j_hashmarks.begin(), j_hashmarks.end()));
    hashmarks.push_back(TiledArray::TiledRange1(x_hashmarks.begin(), x_hashmarks.end()));
    hashmarks.push_back(TiledArray::TiledRange1(y_hashmarks.begin(), y_hashmarks.end()));
    TiledArray::TiledRange ijxy_trange(hashmarks.begin(), hashmarks.end());

    std::shared_ptr<TArray4Tg> result(new TArray4Tg(world_, ijxy_trange) );
    // construct local tiles
    for(auto t=ijxy_trange.tiles().begin();
        t!=ijxy_trange.tiles().end();
        ++t)
      if (result->is_local(*t))
      {
        std::array<std::size_t, 4> index; std::copy(t->begin(), t->end(), index.begin());
        madness::Future < typename TArray4Tg::value_type >
          tile((typename TArray4Tg::value_type(result.get(), index, &tg_s0_gen)
              ));

        // Insert the tile into the array
        result->set(*t, tile);
      }

    tarray4tg_registry_[key] = result;
    v = tarray4tg_registry_.find(key);
    }
    const TArray4Tg& tarray4 = *(v->second);

    std::array<std::string, 4> ind = {{pkey.bra1(), pkey.bra2(), pkey.ket1(), pkey.ket2()}};
    for(auto v=ind.begin(); v!=ind.end(); ++v) {
      if (ParsedTransformedOrbitalSpaceKey::valid_key(*v) ) {
        ParsedTransformedOrbitalSpaceKey ptkey(*v);
        *v = ptkey.label();
      }
    }

    const std::string annotation = ind[0] + "," + ind[1] + "," + ind[2] + "," + ind[3];
    return tarray4(annotation);
  }

  namespace expressions {

    template <typename ArgType, bool Transpose = false>
  struct trace_tensor2_op : public std::binary_function<ArgType, ArgType, TA::Tensor<typename ArgType::numeric_type> > {
      typedef typename ArgType::numeric_type numeric_type;
      typedef TA::Tensor<typename ArgType::numeric_type> result_type;
      typedef result_type value_type;

      static bool is_dense(bool first_dense, bool second_dense) {
        return true;
      }
      static bool is_zero(const bool first_zero, const bool second_zero) {
        return first_zero || second_zero;
      }

      template <typename Left, typename Right>
      static void shape(::TiledArray::detail::Bitset<>& result, const Left& left, const Right& right) {
        result = ::TiledArray::detail::Bitset<>(0);
      }

      // applies the scaling factor
      void scale(numeric_type s) {
        scale_ = s;
      }

      result_type
      operator()(const ArgType& arg1, const ArgType& arg2) const {
        TA_ASSERT(arg1.size() == 1);
        TA_ASSERT(arg2.size() == 1);
        TA_ASSERT(arg1.data()->size() == arg2.data()->size());
        TA_ASSERT(arg1.range() == arg2.range()); // structure of the arguments must be the same

        numeric_type dot_product;
        if (not Transpose) // straight dot
          dot_product = TA::math::dot(arg1.data()->size(), arg1.data()->data(), arg2.data()->data());
        else { // transposed dot done as series of row*column products
          const integer nrow1 = arg1.data()->range().size()[0];
          const integer ncol1 = arg1.data()->range().size()[1];
          TA_ASSERT(ncol1 == arg2.data()->range().size()[0]); // ncol1 == nrow2
          const integer unit_stride = 1;
          dot_product = 0;

          const numeric_type* ptr1 = arg1.data()->data();
          const numeric_type* ptr2 = arg2.data()->data();
          for(integer row1=0; row1<nrow1; ++row1, ptr1+=ncol1, ++ptr2) {
            dot_product += TA::math::dot(ncol1, ptr1, unit_stride, ptr2, nrow1);
          }
        }
        //std::cout << dot_product << std::endl;

        result_type result(arg1.range(), &dot_product);
        return result * scale_;
      }

      result_type operator()(const TiledArray::ZeroTensor&, const ArgType&) {
        TA_ASSERT(false); // Should not be used.
        return result_type();
      }

      result_type operator()(const ArgType&, const TiledArray::ZeroTensor&) {
        TA_ASSERT(false); // Should not be used.
        return result_type();
      }

    private:
      numeric_type scale_;

  };

  template <typename ArgType, bool Transpose = false>
  struct diag_tensor2_op : public std::unary_function<ArgType, TA::Tensor<typename ArgType::numeric_type> > {
      typedef typename ArgType::numeric_type numeric_type;
      typedef TA::Tensor<typename ArgType::numeric_type> result_type;
      typedef result_type value_type;

      static bool is_dense(bool arg_dense) {
        return true;
      }
      static bool is_zero(const bool arg_zero) {
        return arg_zero;
      }

      template <typename Arg>
      static void shape(::TiledArray::detail::Bitset<>& result, const Arg& arg) {
        result = ::TiledArray::detail::Bitset<>(0);
      }

      // applies the scaling factor
      void scale(numeric_type s) {
        scale_ = s;
      }

      result_type
      operator()(const ArgType& arg) const {
        TA_ASSERT(arg.size() == 1);

        numeric_type result_value = 0;
        const auto& ij = arg.range().start();
        const auto& pq_size = arg.data()->range().size();
        if (not Transpose) { // extract i,j
          result_value = arg.data()->data()[ij[0] * pq_size[1] + ij[1]];
        }
        else { // extract j,i
          result_value = arg.data()->data()[ij[1] * pq_size[1] + ij[0]];
        }

        result_type result(arg.range(), &result_value);
        return result * scale_;
      }

      result_type operator()(const TiledArray::ZeroTensor&) {
        TA_ASSERT(false); // Should not be used.
        return result_type();
      }

    private:
      numeric_type scale_;
  };

  }; // end of namespace sc::expressions

  template <typename T>
  std::string
  SingleReference_R12Intermediates<T>::to_space(const std::string& index) {

    bool tformed_space = ParsedTransformedOrbitalSpaceKey::valid_key(index);
    if (tformed_space) {
      ParsedTransformedOrbitalSpaceKey pkey(index);

      // make sure this space exists
      const std::string space_label = to_space(pkey.label());
      std::string space_key = ParsedTransformedOrbitalSpaceKey::key(space_label, pkey.spin(),
                                                                    pkey.support_label(), pkey.support_spin(),
                                                                    pkey.transform_operator());
      auto oreg = r12world_->world()->tfactory()->orbital_registry();
      auto freg = r12world_->world()->fockbuild_runtime();

      if (not oreg->key_exists(space_key)) { // compute if needed
        std::string support_key = ParsedOrbitalSpaceKey::key(to_space(pkey.support_label()), pkey.support_spin());
        std::string target_key = ParsedOrbitalSpaceKey::key(space_label, pkey.spin());

        Ref<OrbitalSpace> target_space = oreg->value(target_key);
        Ref<OrbitalSpace> support_space = oreg->value(support_key);
        RefSCMatrix support_coefs = support_space->coefs();

        RefSCMatrix operator_matrix;
        switch(pkey.transform_operator()) {
          case OneBodyOper::J:
          operator_matrix = freg->get(ParsedOneBodyIntKey::key(support_key, target_key, "J"));
          break;
          case OneBodyOper::K:
          operator_matrix = freg->get(ParsedOneBodyIntKey::key(support_key, target_key, "K"));
          break;
          case OneBodyOper::F:
          operator_matrix = freg->get(ParsedOneBodyIntKey::key(support_key, target_key, "F"));
          break;
          case OneBodyOper::hJ:
          operator_matrix = freg->get(ParsedOneBodyIntKey::key(support_key, target_key, "H"))
                          + freg->get(ParsedOneBodyIntKey::key(support_key, target_key, "J"));
          break;
          case OneBodyOper::gamma:
          operator_matrix = this->rdm1(support_space, target_space);
          break;
          default:
            throw ProgrammingError("SingleReference_R12Intermediates<T>::to_space -- transformed spaces only with Fock-like operators are supported now",
                                   __FILE__, __LINE__);
        }

        Ref<OrbitalSpace> tformed_space = new OrbitalSpace(space_key, space_key,
                                                           target_space, support_coefs * operator_matrix,
                                                           support_space->basis());
        operator_matrix = 0;
        oreg->add(make_keyspace_pair(tformed_space));
      }
      return space_key;
    }

    return to_space_(index);
  }

  template <typename T>
  RefSCMatrix SingleReference_R12Intermediates<T>::rdm1(const Ref<OrbitalSpace>& row_space,
                                                        const Ref<OrbitalSpace>& col_space) const {
    // get the reference 1-RDM
    RefSymmSCMatrix gamma_total = r12world_->refwfn()->ordm_occ_sb(Alpha) + r12world_->refwfn()->ordm_occ_sb(Beta);
    Ref<OrbitalSpace> occ = r12world_->refwfn()->occ_sb(Alpha);
    // map to the target space and support space
    std::vector<int> col_to_occ;
    std::vector<int> row_to_occ;
    try {
      const bool require_same_basis = true;
      col_to_occ = map(*occ, *col_space, require_same_basis);
      row_to_occ = map(*occ, *row_space, require_same_basis);
    }
    catch (sc::CannotConstructMap&) {
      throw ProgrammingError("requested RDM1 in spaces x and y, but x and y are not supported by the same basis",
                             __FILE__, __LINE__);
    }
    RefSCMatrix result = gamma_total.kit()->matrix(row_space->dim(), col_space->dim());
    result.assign(0.0);
    const int nr = result.nrow();
    const int nc = result.ncol();
    for(int r=0; r<nr; ++r) {
      const int rr = row_to_occ[r];
      if (rr >= 0)
        for(int c=0; c<nc; ++c) {
          const int cc = col_to_occ[c];
          if (cc >= 0)
            result.set_element(r, c, gamma_total.get_element(rr, cc));
        }
    }

    return result;
  }

  template<typename T>
  std::string SingleReference_R12Intermediates<T>::to_space_(std::string index) {

    index.erase(std::remove_if(index.begin(), index.end(),
                               [](char c){
                                 return c >= '0' && c<='9';
                               }
                              ), index.end());

    if (index == "i" || index == "j" || index == "k" || index == "l")
      return "i";
    if (index == "m" || index == "n")
      return "m";
    if (index == "i'" || index == "j'" || index == "k'" || index == "l'")
      return "i'";
    if (index == "a" || index == "b" || index == "c" || index == "d")
      return "a";
    if (index == "e" || index == "f")
      return "e";
    if (index == "a" || index == "b" || index == "c" || index == "d")
      return "a";
    if (index == "p" || index == "q" || index == "r" || index == "s")
      return "p";
    if (index == "p'" || index == "q'" || index == "r'" || index == "s'")
      return "p'";
    if (index == "a'" || index == "b'" || index == "c'" || index == "d'")
      return "a'";
    if (index == "A'" || index == "B'" || index == "C'" || index == "D'")
      return "A'";

    std::ostringstream oss;
    oss << "SingleReference_R12Intermediates::to_space(): index label " << index << " not in the current dictionary:";
    throw ProgrammingError(oss.str().c_str(),
                           __FILE__, __LINE__);
  }

  template <typename T>
  std::vector<size_t>
  SingleReference_R12Intermediates<T>::space_hashmarks(std::string space_label) const {

    const bool tformed_space = ParsedTransformedOrbitalSpaceKey::valid_key(space_label);
    if (tformed_space) {
      ParsedTransformedOrbitalSpaceKey ptkey(space_label);
      space_label = ptkey.label();
    }

    auto oreg = r12world_->world()->tfactory()->orbital_registry();
    size_t rank = oreg->value(space_label)->rank();
    std::vector<size_t> hashmarks(2,0); hashmarks[1] = rank;  // all spaces have 1 tile, except
    if (space_label == "i" || space_label == "m") {           // occupied spaces, whose size is determined heuristically for now
                                                              // TODO create a central registry of space tilings?
      const int nproc = r12world_->world()->msg()->n();
      const size_t desired_tilesize = std::max(4 ,
                                               std::min( (int)floor((double)rank / sqrt((double)nproc)), 10)
      );
      const size_t ntiles = (rank+desired_tilesize-1)/desired_tilesize;
      const size_t tilesize = (rank + ntiles - 1)/ ntiles;
      hashmarks.resize(ntiles+1);
      hashmarks[0] = 0;
      for(size_t i=1; i<ntiles; ++i)
        hashmarks[i] = hashmarks[i-1] + tilesize;
      hashmarks.back() = rank;
    }

    return hashmarks;
  }

#if 0
  namespace {
    std::vector<std::string>
    split(const std::string& key,
          char separator) {

      std::stringstream ss(key);
      std::string item;
      std::vector<std::string> tokens;
      while (std::getline(ss, item, separator)) {
        tokens.push_back(item);
      }

      return tokens;
    }

  }

  template <typename T> template <size_t NDIM>
  TA::TiledRange
  SingleReference_R12Intermediates<T>::make_trange(const std::array<std::string, NDIM>& spaces) const {

    std::vector<TA::TiledRange1> hashmarks;
    for(auto id=spaces.begin(); id!=spaces.end(); ++id) {
      std::vector<size_t> hashmark = space_hashmarks(*id);
      hashmarks.push_back(TiledArray::TiledRange1(hashmark.begin(), hashmark.end()));
    }
    return TA::TiledRange(hashmarks.begin(), hashmarks.end());
  }

  namespace detail {

    template <typename T, typename Index>
    void
    sieve_tensor(TA::Array<T, 2>* result,
                 Index t,
                 madness::Future<typename TA::Array<T, 2>::value_type > source_tile,
                 TA::Array<T, 2>* source,
                 Index t_source,
                 std::vector<int>* row_map,
                 std::vector<int>* col_map) {

      typedef TA::Array<T, 2> TArray2;

      auto result_trange = result->trange().make_tile_range(t);
      auto source_trange = source->trange().make_tile_range(t_source);

      std::vector<T> ptr_result(result_trange.volume(), T(0.0));

      for(size_t r=source_trange.start()[0], rc=0;
          r != source_trange.finish()[0];
          ++r) {
        const int rr = row_map->operator[](r);
        if (rr != -1) {
          for(size_t c=source_trange.start()[1];
              c != source_trange.finish()[1];
              ++c, ++rc) {

            const int cc = col_map->operator[](c);

            if (cc != -1) {

              MPQC_ASSERT(false);

            }
          }
        }
      }

      madness::Future < typename TArray2::value_type >
        tile(typename TArray2::value_type(result_trange, &ptr_result[0]));

      // Insert the tile into the array
      result->set(*t, tile);

    }
  }

  template <typename T>
  typename SingleReference_R12Intermediates<T>::TArray2&
  SingleReference_R12Intermediates<T>::sieve(const TArray2& input,
                                             const std::string& output_annotation) {

    // find the input key
    std::string input_key;
    for(auto i = tarray2_registry_.begin();
        i!=tarray2_registry_.end(); ++i) {
      if (i->second.get() == &input) {
        input_key = i->first;
        break;
      }
    }
    MPQC_ASSERT(input_key.empty() == false);

    std::vector<std::string> output_space_keys = split(output_annotation, ',');
    MPQC_ASSERT(output_space_keys.size() == 2);

    ParsedOneBodyIntKey input_pkey(input_key);
    ParsedOneBodyIntKey output_pkey(ParsedOneBodyIntKey::key(output_space_keys[0],
                                                             output_space_keys[1],
                                                             input_pkey.oper()));

    auto v = tarray2_registry_.find(output_pkey.key());
    if (v == tarray2_registry_.end()) {

      // canonicalize indices, may compute spaces if needed
      auto ibra_id = to_space(input_pkey.bra());
      auto iket_id = to_space(input_pkey.ket());
      auto obra_id = to_space(output_pkey.bra());
      auto oket_id = to_space(output_pkey.ket());

      auto oreg = r12world_->world()->tfactory()->orbital_registry();
      auto freg = r12world_->world()->fockbuild_runtime();
      auto ibra = oreg->value(ibra_id);
      auto iket = oreg->value(iket_id);
      auto obra = oreg->value(obra_id);
      auto oket = oreg->value(oket_id);

      auto nbra = obra->rank();
      auto nket = oket->rank();

      std::vector<int> bra_map = map(*obra, *ibra);
      std::vector<int> ket_map = map(*oket, *iket);

      // construct the output array
      std::array<std::string, 2> ospaces = {{obra_id, oket_id}};
      TA::TiledRange output_trange = make_trange(ospaces);
      std::shared_ptr<TArray2> result(new TArray2(world_, output_trange) );

      // construct tasks to fill in local tiles
      for(auto t=result->get_pmap()->begin();
          t!=result->get_pmap()->end();
          ++t)
        {

          // map t to the index of the tile that contains the data needed for t
          auto t_input = *t;
          MPQC_ASSERT(false);

          world_.taskq.add(& detail::sieve_tensor,
                           &result, *t, input.find(t_input), &bra_map, &ket_map);
        }

      tarray2_registry_[output_pkey.key()] = result;
      v = tarray2_registry_.find(output_pkey.key());

    }

    return *(v->second);
  }
#endif

  template <typename T>
  std::ostream& operator<<(std::ostream& os, const DA4_Tile<T>& t) {
    os << typename DA4_Tile<T>::eval_type(t);
    return os;
  }

  template <typename T>
  std::ostream& operator<<(std::ostream& os, const DA4_Tile34<T>& t) {
    os << typename DA4_Tile34<T>::eval_type(t);
    return os;
  }


}; // end of namespace sc

#endif // end of header guard


// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:

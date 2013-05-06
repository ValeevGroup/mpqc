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
    assert(ntasks_with_access == world_.nproc());
    assert(ntasks_with_access == r12world_->world()->msg()->n());

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

    if (darray4->data_persistent()) darray4->deactivate();

    return i1i2_g_p1p2;
  }

  template <typename T>
  typename SingleReference_R12Intermediates<T>::TArray22&
  SingleReference_R12Intermediates<T>::_(const std::string& key) {
    auto v = tarray22_registry_.find(key);
    if (tarray22_registry_.find(key) == tarray22_registry_.end()) {
      ParsedTwoBodyFourCenterIntKey pkey(key);

      /// map operator to the index within the operator set
      Ref<TwoBodyIntDescr> tbint_descr = r12world_->r12tech()->corrfactor()->tbintdescr(r12world_->integral(),0);
      std::string operset_label = "G12'[0]";
      unsigned int oper_idx;
      if (pkey.oper() == "g")
        oper_idx = tbint_descr->intset(TwoBodyOper::eri);
      else if (pkey.oper() == "gr")
        oper_idx = tbint_descr->intset(TwoBodyOper::r12_m1_g12);
      else if (pkey.oper() == "r")
        oper_idx = tbint_descr->intset(TwoBodyOper::r12_0_g12);
      else if (pkey.oper() == "r2") {
        oper_idx = tbint_descr->intset(TwoBodyOper::r12_0_g12);
        operset_label = "G12'[0,0]";
      }
      else if (pkey.oper() == "rTr") {
        oper_idx = tbint_descr->intset(TwoBodyOper::g12t1g12);
        operset_label = "G12'[0,0]";
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
    Ref<TwoBodyIntDescr> tbint_descr = r12world_->r12tech()->corrfactor()->tbintdescr(r12world_->integral(),0);
    std::string operset_label = "G12'[0]";
    unsigned int oper_idx;
    if (pkey.oper() == "g")
      oper_idx = tbint_descr->intset(TwoBodyOper::eri);
    else if (pkey.oper() == "gr")
      oper_idx = tbint_descr->intset(TwoBodyOper::r12_m1_g12);
    else if (pkey.oper() == "r")
      oper_idx = tbint_descr->intset(TwoBodyOper::r12_0_g12);
    else if (pkey.oper() == "r2") {
      oper_idx = tbint_descr->intset(TwoBodyOper::r12_0_g12);
      operset_label = "G12'[0,0]";
    }
    else if (pkey.oper() == "rTr") {
      oper_idx = tbint_descr->intset(TwoBodyOper::g12t1g12);
      operset_label = "G12'[0,0]";
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
    assert(ntasks_with_access == world_.nproc());
    assert(ntasks_with_access == r12world_->world()->msg()->n());

    // make tiled ranges
    // N.B. using little trickery using partial sum to create tilesize=1 range
    std::vector<size_t> i_hashmarks(n1+1, 1);
    i_hashmarks[0] = 0;
    std::partial_sum(i_hashmarks.begin(), i_hashmarks.end(), i_hashmarks.begin());
    std::vector<size_t> j_hashmarks(n2+1, 1);
    j_hashmarks[0] = 0;
    std::partial_sum(j_hashmarks.begin(), j_hashmarks.end(), j_hashmarks.begin());
    std::vector<size_t> x_hashmarks(2, 0); x_hashmarks[1] = n3;
    std::vector<size_t> y_hashmarks(2, 0); y_hashmarks[1] = n4;

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
    Ref<TwoBodyIntDescr> tbint_descr = r12world_->r12tech()->corrfactor()->tbintdescr(r12world_->integral(),0);
    std::string operset_label = "G12'[0]";
    unsigned int oper_idx;
    if (pkey.oper() == "g")
      oper_idx = tbint_descr->intset(TwoBodyOper::eri);
    else if (pkey.oper() == "gr")
      oper_idx = tbint_descr->intset(TwoBodyOper::r12_m1_g12);
    else if (pkey.oper() == "r")
      oper_idx = tbint_descr->intset(TwoBodyOper::r12_0_g12);
    else if (pkey.oper() == "r2") {
      oper_idx = tbint_descr->intset(TwoBodyOper::r12_0_g12);
      operset_label = "G12'[0,0]";
    }
    else if (pkey.oper() == "rTr") {
      oper_idx = tbint_descr->intset(TwoBodyOper::g12t1g12);
      operset_label = "G12'[0,0]";
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

    Ref<TwoBodyMOIntsTransform> tform = r12world_->world()->moints_runtime4()->get(tform_key);
    tform->compute();
    Ref<DistArray4> darray4 = tform->ints_distarray4();

    const size_t n1 = darray4->ni();
    const size_t n2 = darray4->nj();
    const size_t n3 = darray4->nx();
    const size_t n4 = darray4->ny();

    // this only works if all nodes can read integrals
    std::vector<int> tasks_with_access;
    const int ntasks_with_access = darray4->tasks_with_access(tasks_with_access);
    assert(ntasks_with_access == world_.nproc());
    assert(ntasks_with_access == r12world_->world()->msg()->n());

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
  TA::expressions::TensorExpression<TA::Tensor<T> >
  SingleReference_R12Intermediates<T>::_4(const std::string& key) {
    TArray4d& tarray4 = ijxy(key);
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

      result_type operator()(const TiledArray::expressions::ZeroTensor<value_type>&, const ArgType&) {
        TA_ASSERT(false); // Should not be used.
        return result_type();
      }

      result_type operator()(const ArgType&, const TiledArray::expressions::ZeroTensor<value_type>&) {
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

      result_type operator()(const TiledArray::expressions::ZeroTensor<value_type>&) {
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
          operator_matrix = freg->get(ParsedOneBodyIntKey::key(target_key, support_key, "J"));
          break;
          case OneBodyOper::K:
          operator_matrix = freg->get(ParsedOneBodyIntKey::key(target_key, support_key, "K"));
          break;
          case OneBodyOper::F:
          operator_matrix = freg->get(ParsedOneBodyIntKey::key(target_key, support_key, "F"));
          break;
          case OneBodyOper::hJ:
          operator_matrix = freg->get(ParsedOneBodyIntKey::key(target_key, support_key, "H"))
                          + freg->get(ParsedOneBodyIntKey::key(target_key, support_key, "J"));
          break;
          default:
            throw ProgrammingError("SingleReference_R12Intermediates<T>::to_space -- transformed spaces only with Fock-like operators are supported now",
                                   __FILE__, __LINE__);
        }

        Ref<OrbitalSpace> tformed_space = new OrbitalSpace(space_key, space_key,
                                                           target_space, support_coefs * operator_matrix.t(),
                                                           support_space->basis());
        oreg->add(make_keyspace_pair(tformed_space));
      }
      return space_key;
    }

    if (index == "i" || index == "j" || index == "k" || index == "l")
      return "i";
    if (index == "m" || index == "n")
      return "m";
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
    oss << "SingleReference_R12Intermediates::to_space(): index label " << index << " not in the dictionary";
    throw ProgrammingError(oss.str().c_str(),
                           __FILE__, __LINE__);
  }

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

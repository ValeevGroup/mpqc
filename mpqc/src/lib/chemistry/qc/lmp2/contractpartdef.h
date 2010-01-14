
/*
 * Copyright 2009 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 *
 * This file is a part of the MPQC LMP2 library.
 *
 * The MPQC LMP2 library is free software: you can redistribute it
 * and/or modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/>.
 *
 */

#ifndef _chemistry_qc_lmp2_contractpartdef_h
#define _chemistry_qc_lmp2_contractpartdef_h

namespace sc {

namespace sma2 {

class SumOperation {
  public:
    inline void operator()(double &res, double val) { res += val; }
};

class DivOperation {
  public:
    inline void operator()(double &res, double val) { res /= val; }
};

template <int N, int N2>
void get_all_indices(
    const ContractPart<N> &array,
    const ContractPart<N2> &other,
    std::vector<int> &fixed, std::vector<int> &fixed_values,
    std::vector<int> &external, std::vector<int> &external_on_other,
    bool same_external_index_order_required)
{
  for (int i=0; i<array.array().nindex(); i++) {
      bool found = false, isfixed = false;
      // see if the index is fixed
      if (array.index(i).has_value()) {
          isfixed = true;
          if (i != fixed.size()) {
              throw std::runtime_error("fixed indices must be first");
            }
          fixed.push_back(i);
          fixed_values.push_back(array.index(i).value());
        }
      // see if the index is external
      // note: it is possible for an index to be both fixed and external
      for (int j=0; j<other.array().nindex(); j++) {
          if (array.index(i).symbolically_equivalent(other.index(j))) {
              external.push_back(i);
              external_on_other.push_back(j);
              found = true;
              break;
            }
        }
      if (isfixed && !found) {
          if (!array.array().index(i).all_size_one()) {
              throw std::runtime_error("nonexternal fixed indices require blocksize == 1");
            }
        }
      if (!found && !isfixed) {
          throw std::runtime_error("index not found");
        }
    }

  if (same_external_index_order_required) {
      for (int i=1; i<external_on_other.size(); i++) {
          if (external_on_other[i] < external_on_other[i-1]) {
              throw
                  std::runtime_error("external indices must be in same order");
            }
        }
    }
}

template <int N>
template <int N2, class Op>
void ContractPart<N>::do_binary_op_with_fixed_indices(
    double f, const ContractPart<N2> &o, bool initarray, Op &op) const
{
  double alpha = f * o.factor() / factor();
  if (initarray) {
      array_.allocate_blocks();
      array_.zero();
    }

  // compute the index arrays for C
  sma2::Array<N> &C = array_;
  std::vector<int> C_fixed, C_fixed_values;
  std::vector<int> C_external, C_external_on_A;
  get_all_indices(*this, o, C_fixed, C_fixed_values, C_external,
                  C_external_on_A, false);

  // compute the index arrays for A
  const sma2::Array<N2> &A = o.array();
  std::vector<int> A_fixed, A_fixed_values;
  std::vector<int> A_external, A_external_on_C;
  get_all_indices(o, *this, A_fixed, A_fixed_values, A_external,
                  A_external_on_C, false);

  // Fill in the blocks in A that are fixed.  The remaining
  // indices will be filled from C's external indices.
  sma2::BlockInfo<N2> A_bi;
  for (int i=0; i<A_fixed.size(); i++) {
//        std::cout << "assigning A_bi.block(" << A_fixed[i] << ") to "
//                  << A_fixed_values[i]
//                  << std::endl;
      A_bi.block(A_fixed[i]) = A_fixed_values[i];
    }

  // These index lists are needed to set up A's BlockInfo
  sma2::IndexList A_external_il(A_external);
  sma2::IndexList A_external_on_C_il(A_external_on_C);

  // find the range for the blocks in C that have the correct
  // fixed indices
  sma2::BlockInfo<N> first_C_bi;
  sma2::BlockInfo<N> last_C_bi;
  for (int i=0; i<N; i++) {
      first_C_bi.block(i) = 0;
      last_C_bi.block(i) = C.index(i).nindex();
    }
  for (int i=0; i<C_fixed.size(); i++) {
      first_C_bi.block(C_fixed[i]) = C_fixed_values[i];
      last_C_bi.block(C_fixed[i]) = C_fixed_values[i];
    }
  typename sma2::Array<N>::blockmap_t::const_iterator
      C_iter_begin = C.blockmap().lower_bound(first_C_bi);
  typename sma2::Array<N>::blockmap_t::const_iterator
      C_iter_fence = C.blockmap().upper_bound(last_C_bi);

  bool same_index_order = A_external_on_C_il.is_identity();

//    std::cout << "in sum routine" << std::endl;

//    std::cout << "A_external_on_C_il:" << std::endl;
//    for (int i=0; i<A_external_on_C_il.n(); i++) {
//        std::cout << " " << A_external_on_C_il.i(i);
//      }
//    std::cout << std::endl;

//    std::cout << "A_external_il:" << std::endl;
//    for (int i=0; i<A_external_il.n(); i++) {
//        std::cout << " " << A_external_il.i(i);
//      }
//    std::cout << std::endl;

  // iterate through relevant blocks in C
  for (typename sma2::Array<N>::blockmap_t::const_iterator C_iter = C_iter_begin;
       C_iter != C_iter_fence;
       C_iter++) {
      const sma2::BlockInfo<N> &C_bi = C_iter->first;
      // find the contributing block in A
      A_bi.assign_blocks(A_external_il, C_bi, A_external_on_C_il);
//        std::cout << "C block " << C_bi
//                  << " A block " << A_bi
//                  << std::endl;
      typename sma2::Array<N2>::blockmap_t::const_iterator
          A_iter = A.blockmap().find(A_bi);
      if (A_iter == A.blockmap().end()) continue;
      double *C_data = C_iter->second;
      double *A_data = A_iter->second;
      int sz = C.block_size(C_bi);
//        std::cout << " sz = " << sz << std::endl;
      if (same_index_order) {
          for (int i=0; i<sz; i++) {
              op(C_data[i], alpha * A_data[i]);
            }
        }
      else {
          BlockIter<N> cbiter(C.indices(),C_bi);
          int coff = 0;
          for (cbiter.start(); cbiter.ready(); cbiter++,coff++) {
              int aoff = cbiter.subset_offset(A_external_on_C_il);
              op(C_data[coff], alpha * A_data[aoff]);
            }
        }
    }

  if (!skip_bounds_update_) C.compute_bounds();
}

template <int N>
template <class Op>
void ContractPart<N>::do_binary_op(double f, const ContractPart<N> &o,
                            bool initarray, Op &op) const {
  int nvalue = 0;
  for (int i=0; i<N; i++) {
      if (index(i).has_value()) nvalue++;
    }
  for (int i=0; i<N; i++) {
      if (o.index(i).has_value()) nvalue++;
    }
  if (nvalue) {
      do_binary_op_with_fixed_indices(f,o,initarray,op);
      return;
    }

  std::vector<int> indvec;
  for (int i=0; i<N; i++) {
      for (int j=0; j<N; j++) {
          if (o.index(i) == index(j)) {
              indvec.push_back(j);
              break;
            }
        }
    }
  if (indvec.size() != N) {
      throw std::invalid_argument("sma::ContractPart: nindex != N");
    }
  double alpha = f * o.factor() / factor();
  IndexList indlist(indvec);
  if (initarray) {
      array_.allocate_blocks();
      array_.zero();
      // the blocks and indices should already be initialized
      // so this shouldn't be needed.
      for (int i=0; i<N; i++)
          array_.set_index(indlist.i(i),o.array().index(i));
    }
  binary_op(array_, alpha, o.array(), indlist, skip_bounds_update_, op);
}

template <int N>
template <int Nl,int Nr>
void ContractPart<N>::doprod(double f, const ContractProd<Nl,Nr> &o,
                            bool initarray) const {
  std::vector<int> extCAv, extCBv, extAv, extBv, intAv, intBv;
  std::vector<int> fixextCAv, fixextCBv, fixextAv, fixextBv;

  for (int i=0; i<N; i++) {
      for (int j=0; j<Nl; j++) {
          if (indices_[i].symbolically_equivalent(o.l.index(j))) {
              extAv.push_back(j);
              extCAv.push_back(i);
              if (o.l.index(j).has_value() != indices_[i].has_value()
                  || (indices_[i].has_value()
                      && indices_[i].value() != o.l.index(j).value())) {
                  throw std::runtime_error("ContractPart<N>::doprod: two equivalent indices do not have the same value");
                }
              if (indices_[i].has_value()) {
                  fixextAv.push_back(j);
                  fixextCAv.push_back(i);
                }
            }
        }
    }

  for (int i=0; i<N; i++) {
      for (int j=0; j<Nr; j++) {
          if (indices_[i].symbolically_equivalent(o.r.index(j))) {
              extBv.push_back(j);
              extCBv.push_back(i);
              if (o.r.index(j).has_value() != indices_[i].has_value()
                  || (indices_[i].has_value()
                      && indices_[i].value() != o.r.index(j).value())) {
                  throw std::runtime_error("ContractPart::doprod: two equivalent indices do not have the same value");
                }
              if (indices_[i].has_value()) {
                  fixextBv.push_back(j);
                  fixextCBv.push_back(i);
                }
            }
        }
    }

  for (int i=0; i<Nl; i++) {
      for (int j=0; j<Nr; j++) {
          if (o.l.index(i).symbolically_equivalent(o.r.index(j))) {
              intAv.push_back(i);
              intBv.push_back(j);
            }
        }
    }

  IndexList extCA(extCAv), extCB(extCBv);
  IndexList intA(intAv), extA(extAv);
  IndexList intB(intBv), extB(extBv);
  IndexList fixextCA(fixextCAv), fixextCB(fixextCBv);
  IndexList fixextA(fixextAv), fixextB(fixextBv);

  // Find the values of the fixed indices.
  std::vector<int> fixCv, fixAv, fixBv;
  std::vector<sma2::bi_t> fixvalCv, fixvalAv, fixvalBv;
  for (int i=0; i<N; i++) {
      if (indices_[i].has_value()) {
          fixCv.push_back(i);
          fixvalCv.push_back(indices_[i].value());
        }
    }
  for (int i=0; i<Nl; i++) {
      if (o.l.index(i).has_value()) {
          fixAv.push_back(i);
          fixvalAv.push_back(o.l.index(i).value());
        }
    }
  for (int i=0; i<Nr; i++) {
      if (o.r.index(i).has_value()) {
          fixBv.push_back(i);
          fixvalBv.push_back(o.r.index(i).value());
        }
    }

  IndexList fixC(fixCv), fixA(fixAv), fixB(fixBv);
  BlockInfo<N> fixvalC(fixvalCv);
  BlockInfo<Nl> fixvalA(fixvalAv);
  BlockInfo<Nr> fixvalB(fixvalBv);

  double ABfactor = f * o.l.factor() * o.r.factor() / factor();

  if (initarray) {
      array_.allocate_blocks();
      array_.zero();
      for (int i=0; i<extCA.n(); i++) {
          array_.set_index(extCA.i(i), o.l.array().index(extA.i(i)));
        }
      for (int i=0; i<extCB.n(); i++) {
          array_.set_index(extCB.i(i), o.r.array().index(extB.i(i)));
        }
    }

  contract(array_, extCA, extCB, fixextCA, fixextCB, fixC, fixvalC,
           o.l.array(), extA, fixextA, intA, fixA, fixvalA, o.l.clear_after_use(),
           o.r.array(), extB, fixextB, intB, fixB, fixvalB, o.r.clear_after_use(),
           ABfactor, initarray, o.regtimer);
}

template <int N>
template <int Nl,int Nr>
void ContractPart<N>::dounion(const ContractProd<Nl,Nr> &o) const {
  std::vector<int> extCAv, extCBv, extAv, extBv, intAv, intBv;

  for (int i=0; i<N; i++) {
      for (int j=0; j<Nl; j++) {
          if (indices_[i].symbolically_equivalent(o.l.index(j))) {
              extAv.push_back(j);
              extCAv.push_back(i);
            }
        }
    }

  for (int i=0; i<N; i++) {
      for (int j=0; j<Nr; j++) {
          if (indices_[i].symbolically_equivalent(o.r.index(j))) {
              extBv.push_back(j);
              extCBv.push_back(i);
            }
        }
    }

  for (int i=0; i<Nl; i++) {
      for (int j=0; j<Nr; j++) {
          if (o.l.index(i).symbolically_equivalent(o.r.index(j))) {
              intAv.push_back(i);
              intBv.push_back(j);
            }
        }
    }

  IndexList extCA(extCAv), extCB(extCBv);
  IndexList intA(intAv), extA(extAv);
  IndexList intB(intBv), extB(extBv);

  std::vector<int> fixCv, fixAv, fixBv;
  std::vector<sma2::bi_t> fixvalCv, fixvalAv, fixvalBv;
  for (int i=0; i<N; i++) {
      if (indices_[i].has_value()) {
          if (fixCv.size() != i)
              throw std::invalid_argument("fixed indices must be first (in C)");
          fixCv.push_back(i);
          fixvalCv.push_back(indices_[i].value());
          if (fixvalCv[i] >= array().index(i).nblock()) {
              throw std::invalid_argument("fixed index out of range (in C)");
            }
        }
    }
  for (int i=0; i<Nl; i++) {
      if (o.l.index(i).has_value()) {
          if (fixAv.size() != i)
              throw std::invalid_argument("fixed indices must be first (in A)");
          fixAv.push_back(i);
          fixvalAv.push_back(o.l.index(i).value());
          if (fixvalAv[i] >= o.l.array().index(i).nblock()) {
              throw std::invalid_argument("fixed index out of range (in A)");
            }
        }
    }
  for (int i=0; i<Nr; i++) {
      if (o.r.index(i).has_value()) {
          fixBv.push_back(i);
          fixvalBv.push_back(o.r.index(i).value());
          if (fixvalBv[i] >= o.r.array().index(i).nblock()) {
              throw std::invalid_argument("fixed index out of range (in B)");
            }
        }
    }

  IndexList fixC(fixCv), fixA(fixAv), fixB(fixBv);
  BlockInfo<N> fixvalC(fixvalCv);
  BlockInfo<Nl> fixvalA(fixvalAv);
  BlockInfo<Nr> fixvalB(fixvalBv);

  bool bad_index = false;
  for (int i=0; i<extCA.n(); i++) {
      if (array_.index(extCA.i(i)) != o.l.array().index(extA.i(i)))
          bad_index = true;
      //array_.set_index(extCA.i(i), o.l.array().index(extA.i(i)));
    }
  for (int i=0; i<extCB.n(); i++) {
      if (array_.index(extCB.i(i)) != o.r.array().index(extB.i(i)))
          bad_index = true;
      //array_.set_index(extCB.i(i), o.r.array().index(extB.i(i)));
    }
  if (bad_index) {
      throw std::invalid_argument("sma2::contract_union: C range inconsistency");
    }

  contract_union(array_, extCA, extCB, fixC, fixvalC,
                 o.l.array(), extA, intA, fixA, fixvalA,
                 o.r.array(), extB, intB, fixB, fixvalB);
}

template <int N>
ContractPart<N>::ContractPart(Array<N> &array):
  array_(array), factor_(1.0), clear_after_use_(false),
  skip_bounds_update_(false) {
  if (N != 0) throw std::invalid_argument("sma::ContractPart: N != 0");
}

template <int N>
ContractPart<N>::ContractPart(Array<N> &array,
                              const Index &i1):
  array_(array), factor_(1.0), clear_after_use_(false),
  skip_bounds_update_(false) {
  if (N != 1) throw std::invalid_argument("sma::ContractPart: N != 1");
  indices_.push_back(i1);
}

template <int N>
ContractPart<N>::ContractPart(Array<N> &array,
                              const Index &i1,
                              const Index &i2):
  array_(array), factor_(1.0), clear_after_use_(false),
  skip_bounds_update_(false) {
  if (N != 2) throw std::invalid_argument("sma::ContractPart: N != 2");
  indices_.push_back(i1);
  indices_.push_back(i2);
}

template <int N>
ContractPart<N>::ContractPart(Array<N> &array,
                              const Index &i1,
                              const Index &i2,
                              const Index &i3):
  array_(array), factor_(1.0), clear_after_use_(false),
  skip_bounds_update_(false) {
  if (N != 3) throw std::invalid_argument("sma::ContractPart: N != 3");
  indices_.push_back(i1);
  indices_.push_back(i2);
  indices_.push_back(i3);
}

template <int N>
ContractPart<N>::ContractPart(Array<N> &array,
                              const Index &i1,
                              const Index &i2,
                              const Index &i3,
                              const Index &i4):
  array_(array), factor_(1.0), clear_after_use_(false),
  skip_bounds_update_(false) {
  if (N != 4) throw std::invalid_argument("sma::ContractPart: N != 4");
  indices_.push_back(i1);
  indices_.push_back(i2);
  indices_.push_back(i3);
  indices_.push_back(i4);
}

template <int N>
ContractPart<N>::ContractPart(Array<N> &array,
                              const Index &i1,
                              const Index &i2,
                              const Index &i3,
                              const Index &i4,
                              const Index &i5):
  array_(array), factor_(1.0), clear_after_use_(false),
  skip_bounds_update_(false) {
  if (N != 5) throw std::invalid_argument("sma::ContractPart: N != 5");
  indices_.push_back(i1);
  indices_.push_back(i2);
  indices_.push_back(i3);
  indices_.push_back(i4);
  indices_.push_back(i5);
}

template <int N>
ContractPart<N>::ContractPart(Array<N> &array,
                              const Index &i1,
                              const Index &i2,
                              const Index &i3,
                              const Index &i4,
                              const Index &i5,
                              const Index &i6):
  array_(array), factor_(1.0), clear_after_use_(false),
  skip_bounds_update_(false) {
  if (N != 6) throw std::invalid_argument("sma::ContractPart: N != 6");
  indices_.push_back(i1);
  indices_.push_back(i2);
  indices_.push_back(i3);
  indices_.push_back(i4);
  indices_.push_back(i5);
  indices_.push_back(i6);
}

template <int N>
void ContractPart<N>::apply_factor(double f)
{
  factor_ *= f;
}

template <int N>
double ContractPart<N>::factor() const
{
  return factor_;
}

template <int N>
Array<N>& ContractPart<N>::array() const { return array_; }

template <int N>
const Index &ContractPart<N>::index(int i) const { return indices_[i]; }

template <int N>
template <int N2>
void ContractPart<N>::operator |= (const ContractPart<N2> &o) const
{
  // compute the index arrays for C
  sma2::Array<N> &C = array_;
  std::vector<int> C_fixed, C_fixed_values;
  std::vector<int> C_external, C_external_on_A;
  get_all_indices(*this, o, C_fixed, C_fixed_values, C_external,
                  C_external_on_A, false);

  // compute the index arrays for A
  const sma2::Array<N2> &A = o.array();
  std::vector<int> A_fixed, A_fixed_values;
  std::vector<int> A_external, A_external_on_C;
  get_all_indices(o, *this, A_fixed, A_fixed_values, A_external,
                  A_external_on_C, false);

  for (int i=0; i<A_external.size(); i++) {
      if (A.index(A_external[i]) != C.index(A_external_on_C[i])) {
          throw std::invalid_argument("|=: Range conflict between A and C");
        }
    }

  // Fill in the blocks in C that are fixed.  The remaining
  // indices will be filled from A's external indices.
  sma2::BlockInfo<N> C_bi;
  for (int i=0; i<C_fixed.size(); i++) {
      C_bi.block(C_fixed[i]) = C_fixed_values[i];
    }

  // These index lists are needed to set up A's BlockInfo
  sma2::IndexList A_external_il(A_external);
  sma2::IndexList C_external_il(A_external_on_C);

  // find the range for the blocks in A that have the correct
  // fixed indices
  sma2::BlockInfo<N2> first_A_bi;
  sma2::BlockInfo<N2> last_A_bi;
  for (int i=0; i<N2; i++) {
      first_A_bi.block(i) = 0;
      last_A_bi.block(i) = A.index(i).nindex();
    }
  for (int i=0; i<A_fixed.size(); i++) {
      first_A_bi.block(A_fixed[i]) = A_fixed_values[i];
      last_A_bi.block(A_fixed[i]) = A_fixed_values[i];
    }
  typename sma2::Array<N2>::blockmap_t::const_iterator
      A_iter_begin = A.blockmap().lower_bound(first_A_bi);
  typename sma2::Array<N2>::blockmap_t::const_iterator
      A_iter_fence = A.blockmap().upper_bound(last_A_bi);

  // iterate through relevant blocks in A
  for (typename sma2::Array<N2>::blockmap_t::const_iterator A_iter = A_iter_begin;
       A_iter != A_iter_fence;
       A_iter++) {
      const typename sma2::BlockInfo<N2> &A_bi = A_iter->first;
      // find the contributing block in C
      C_bi.assign_blocks(C_external_il, A_bi, A_external_il);
      C.add_unallocated_block(C_bi);
    }
}

template <int N>
void ContractPart<N>::operator = (const ContractPart<N> &o) const {
  SumOperation op;
  do_binary_op(1.0, o, true, op);
}

template <int N>
template <int N2>
void ContractPart<N>::operator = (const ContractPart<N2> &o) const {
  SumOperation op;
  do_binary_op_with_fixed_indices(1.0,o,true,op);
}

template <int N>
void ContractPart<N>::operator += (const ContractPart<N> &o) const {
  SumOperation op;
  do_binary_op(1.0, o, false, op);
}

template <int N>
template <int N2>
void ContractPart<N>::operator += (const ContractPart<N2> &o) const {
  SumOperation op;
  do_binary_op_with_fixed_indices(1.0,o,false,op);
}

template <int N>
void ContractPart<N>::operator /= (const ContractPart<N> &o) const {
  DivOperation op;
  do_binary_op(1.0, o, false, op);
}

template <int N>
template <int N2>
void ContractPart<N>::operator /= (const ContractPart<N2> &o) const {
  DivOperation op;
  do_binary_op_with_fixed_indices(1.0,o,false,op);
}

template <int N>
void ContractPart<N>::operator -= (const ContractPart<N> &o) const {
  SumOperation op;
  do_binary_op(-1.0, o, false, op);
}

template <int N>
template <int Nl,int Nr>
void ContractPart<N>::operator = (const ContractProd<Nl,Nr> &o) const {
  doprod(1.0, o, true);
}

template <int N>
template <int Nl,int Nr>
void ContractPart<N>::operator += (const ContractProd<Nl,Nr> &o) const{
  doprod(1.0, o, false);
}

template <int N>
template <int Nl,int Nr>
void ContractPart<N>::operator -= (const ContractProd<Nl,Nr> &o) const{
  doprod(-1.0, o, false);
}

template <int N>
template <int Nl,int Nr>
void ContractPart<N>::operator |= (const ContractProd<Nl,Nr> &o) const {
  dounion(o);
}

template <int N>
ContractPart<N> operator *(double f, const ContractPart<N> &o)
{
  ContractPart<N> result(o);
  result.apply_factor(f);
  return result;
}

template <int N>
ContractPart<N> ContractPart<N>::operator ~ () const
{
  ContractPart<N> result(*this);
  result.clear_after_use_ = true;
  return result;
}

template <int N>
ContractPart<N> ContractPart<N>::skip_bounds_update () const
{
  ContractPart<N> result(*this);
  result.skip_bounds_update_ = true;
  return result;
}

template <int N>
double
ContractPart<N>::value()
{
  if (indices_.size() != N) {
      throw std::runtime_error("ContractPart::value: indices_.size() != N");
    }
  BlockInfo<N> bi;
  for (int i=0; i<N; i++) {
      if (!indices_[i].has_value()) {
          throw std::runtime_error("ContractPart::value: requires fixed indices");
        }
      bi.block(i) = indices_[i].value();
    }
  size_t blocksize = array_.block_size(bi);
  if (blocksize > 1) {
      throw std::runtime_error("ContractPart::value: blocksize must be 1");
    }
  typename Array<N>::blockmap_t::const_iterator
      biter = array_.blockmap().find(bi);
  if (biter == array_.blockmap().end()) {
      return 0.0;
    }
  if (biter->second == 0) {
      throw std::runtime_error("ContractPart::value: array not allocated");
    }
  return *biter->second;
}

}

}

#endif

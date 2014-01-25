
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

#ifndef _chemistry_qc_lmp2_contract_h
#define _chemistry_qc_lmp2_contract_h

#include <chemistry/qc/lmp2/dgemminfo.h>

#include <stdexcept>

#include <math/scmat/blas.h>

#define USE_BOUNDS_IN_CONTRACT_UNION 1

namespace sc {

namespace sma2 {

/// Returns true if an array needs to be repacked to perform a contraction.
/// @param A an Array object.
/// @param row_indices the indices in that form the rows of the contraction matrix.
/// @param col_indices the indices in that from the columns of the contraction matrix.
/// @param fixed a list of indices that will not vary
/// @param fixedvals the values of the fixed indices (must be have a block size of one)
template <int N>
inline bool
need_repack(const Array<N> *A, const IndexList &row_indices, const IndexList &col_indices,
            const IndexList &fixed, const BlockInfo<N> *fixedvals)
{
  // All that matters are indices that have block sizes greater than 1.
  // If those indices are ordered correctly, then a repack is not needed.

  std::vector<int> indices;

  for (int i=0; i<row_indices.n(); i++) {
      int index = row_indices.i(i);
      if (A->index(index).max_block_size() > 1) {
          indices.push_back(index);
        }
    }

  for (int i=0; i<col_indices.n(); i++) {
      int index = col_indices.i(i);
      if (A->index(index).max_block_size() > 1) {
          indices.push_back(index);
        }
    }

  for (int i=1; i<indices.size(); i++) {
      if (indices[i] < indices[i-1]) return true;
    }

  return false;
}

/// Returns an relative estimate of the cost to repack the matrix.
/// @param A an Array object.
/// @param row_indices the indices in that form the rows of the contraction matrix.
/// @param col_indices the indices in that from the columns of the contraction matrix.
/// @param fixed a list of indices that will not vary.
/// @param fixedvals the values of the fixed indices (must be have a block size of one).
template <int N>
inline double
repack_cost(const Array<N> *A,
            const IndexList &row_indices, const IndexList &col_indices,
            const IndexList &fixed, const BlockInfo<N> *fixedvals)
{
  // This gives a very rough estimate of the relative repack costs
  // nindex cost
  // 0      0
  // 1      1
  // 2      3
  // 3      7
  // 4      15
  double cost = A->n_element_allocated();
  for (int i=0; i<fixed.n(); i++) {
      cost = cost/A->index(fixed.i(i)).nblock();
    }
  return cost;
}

/// Determine the cost of repacking arrays for a contraction.
template <int NC, int NA, int NB>
class RepackScheme {
    double cost_;
    bool need_repack_A_, transpose_A_;
    bool need_repack_B_, transpose_B_;
    bool need_repack_C_, transpose_C_;
    IndexList extA_, intA_, fixA_;
    const BlockInfo<NA> *fixvalA_;
    IndexList extB_, intB_, fixB_;
    const BlockInfo<NB> *fixvalB_;
    IndexList extCA_, extCB_, fixC_;
    const BlockInfo<NC> *fixvalC_;
    const Array<NA> *A_;
    const Array<NB> *B_;
    const Array<NC> *C_;
    int n_C_repack_;

    void init() {
      cost_ = 0.0;
      need_repack_A_ = false;
      transpose_A_ = false;
      if (need_repack(A_, extA_, intA_, fixA_, fixvalA_)) {
          if (need_repack(A_, intA_, extA_, fixA_, fixvalA_)) {
              need_repack_A_ = true;
              cost_ += 2.0 * repack_cost(A_,extA_,intA_,fixA_,fixvalA_);
            }
          else {
              transpose_A_ = true;
            }
        }

      need_repack_B_ = false;
      transpose_B_ = false;
      if (need_repack(B_, intB_, extB_, fixB_, fixvalB_)) {
          if (need_repack(B_, extB_, intB_, fixB_, fixvalB_)) {
              need_repack_B_ = true;
              cost_ += 2.0 * repack_cost(B_,intB_,extB_,fixB_,fixvalB_);
            }
          else {
              transpose_B_ = true;
            }
        }

      need_repack_C_ = false;
      transpose_C_ = false;
      if (need_repack(C_, extCA_, extCB_, fixC_, fixvalC_)) {
          if (need_repack(C_, extCB_, extCA_, fixC_, fixvalC_)) {
              need_repack_C_ = true;
              cost_ += n_C_repack_
                     * repack_cost(C_,extCA_,extCB_,fixC_,fixvalC_);
            }
          else {
              transpose_C_ = true;
            }
        }
    }
    void reorder_indices(IndexList &i1,IndexList &i2) {
      std::map<int,int> index_map;
      for (int i=0; i<i1.n(); i++) {
          index_map[i1.i(i)] = i2.i(i);
        }
      int i=0;
      for (std::map<int,int>::iterator
               iter=index_map.begin();
           iter!=index_map.end(); i++,iter++) {
          i1.i(i) = iter->first;
          i2.i(i) = iter->second;
        }
    }
  public:
    /// Create the RepackScheme for a given contraction.
    RepackScheme(const Array<NA> *A,
                 const IndexList &extA, const IndexList &intA,
                 const IndexList &fixA, const BlockInfo<NA> *fixvalA,
                 const Array<NB> *B,
                 const IndexList &intB, const IndexList &extB,
                 const IndexList &fixB, const BlockInfo<NB> *fixvalB,
                 const Array<NC> *C,
                 const IndexList &extCA, const IndexList &extCB,
                 const IndexList &fixC, const BlockInfo<NC> *fixvalC,
                 int n_C_repack
        ):
      cost_(0.0),
      extA_(extA), intA_(intA), fixA_(fixA), fixvalA_(fixvalA),
      extB_(extB), intB_(intB), fixB_(fixB), fixvalB_(fixvalB),
      extCA_(extCA), extCB_(extCB), fixC_(fixC), fixvalC_(fixvalC),
      A_(A), B_(B), C_(C),
      n_C_repack_(n_C_repack) {
      init();
    }
    /// Set the array that determines the index ordering to C (i==0), A (i==1), or B (i==2)
    /// This will update the cost.
    void driver(int i) {
      if (i==0) {
          reorder_indices(extCA_,extA_);
          reorder_indices(extCB_,extB_);
        }
      else if (i==1) {
          reorder_indices(extA_,extCA_);
          reorder_indices(intA_,intB_);
        }
      else {
          reorder_indices(extB_,extCB_);
          reorder_indices(intB_,intA_);
        }
      init();
    }
    /// Returns true if A needs repacked in the current scheme.
    bool need_repack_A() const { return need_repack_A_; }
    /// Returns true if A needs transposed in the current scheme.
    bool transpose_A() const { return transpose_A_; }
    /// Returns true if B needs repacked in the current scheme.
    bool need_repack_B() const { return need_repack_B_; }
    /// Returns true if B needs transposed in the current scheme.
    bool transpose_B() const { return transpose_B_; }
    /// Returns true if C needs repacked in the current scheme.
    bool need_repack_C() const { return need_repack_C_; }
    /// Returns true if C needs transposed in the current scheme.
    bool transpose_C() const { return transpose_C_; }
    /// Returns the cost of the current scheme.
    double cost() const { return cost_; }

    /// Assign the contraction indices.
    /// This is not update the cost.
    void assign_indices(IndexList &extA, IndexList &intA,
                        IndexList &intB, IndexList &extB,
                        IndexList &extCA,IndexList &extCB) {
      extA  = extA_;
      intA  = intA_;
      intB  = intB_;
      extB  = extB_;
      extCA = extCA_;
      extCB = extCB_;
    }
};

/// Repack an array to prepare it for a contraction.
/// @param A an Array object.
/// @param row_indices the indices in that form the rows of the contraction matrix.
/// @param col_indices the indices in that from the columns of the contraction matrix.
/// @param fixed a list of indices that will not vary.
/// @param fixedvals the values of the fixed indices (must be have a block size of one).
/// @param reverse perform the reverse packing operation.
template <int N>
inline void
repack(Array<N> &A, const IndexList &row_indices, const IndexList &col_indices,
       const IndexList &fixed, const BlockInfo<N> &fixedvals,
       bool reverse = false)
{
  double *tmp_data = 0;
  int n_tmp_data = 0;

  const typename Array<N>::blockmap_t &amap = A.blockmap();

  typename Array<N>::blockmap_t::const_iterator begin, end;
  if (fixed.n() == 0) {
      begin = amap.begin();
      end = amap.end();
    }
  else {
      BlockInfo<N> sbi;

      sbi.zero();
      sbi.assign_blocks(fixed, fixedvals);
#ifdef USE_BOUND
      sbi.set_bound(DBL_MAX);
#endif
      begin = amap.lower_bound(sbi);

      for (int i=0; i<N; i++) sbi.block(i) = A.index(i).nblock();
      sbi.assign_blocks(fixed, fixedvals);
#ifdef USE_BOUND
      sbi.set_bound(0.0);
#endif
      end = amap.upper_bound(sbi);
    }

  for (typename Array<N>::blockmap_t::const_iterator aiter = begin;
       aiter != end;
       aiter++) {
      const BlockInfo<N> &bi = aiter->first;
      double *data = aiter->second;
      int ndata = bi.size(A.indices());
      if (n_tmp_data < ndata) {
          delete[] tmp_data;
          tmp_data = new double[ndata];
          n_tmp_data = ndata;
        }
      int nrow = bi.subset_size(A.indices(), row_indices);
      int ncol = bi.subset_size(A.indices(), col_indices);
      if (ndata != nrow * ncol) {
          throw std::length_error("sma::repack: ntotal != nrow * ncol");
        }
      memcpy(tmp_data,data,sizeof(double)*ndata);
      BlockIter<N> a_iter(A.indices(), bi);
      for (a_iter.start(); a_iter.ready(); a_iter++) {
          int row_index = a_iter.subset_offset(row_indices);
          int col_index = a_iter.subset_offset(col_indices);
          int index = a_iter.offset();
          int repacked_index = row_index*ncol + col_index;
          if (!reverse) {
              data[repacked_index] = tmp_data[index];
            }
          else {
              data[index] = tmp_data[repacked_index];
            }
        }
    }

  delete[] tmp_data;
}

/** Perform a contraction.
    The contraction C += f * A * B is computed, where C, A, and B
    are arrays and f is a scalar.  A and B can have fixed indices;
    however, the fixed indices must come before all other indices.
    @param C the C array.
    @param c_extCA the indices of C that are shared with A.
    @param c_extCB the indices of C that are shared with B.
    @param fixextCA the fixed indices of C that are shared with A.
    @param fixextCB the fixed indices of C that are shared with B.
    @param fixC all of the fixed indices of C.
    @param fixvalC the values of all the fixed indices of C.
    @param A the A array.
    @param c_extA the indices of A that are shared with C.
    @param fixextA the fixed indices of A that are shared with C.
    @param c_intA the internal indices (shared between A and B).
    @param fixA the fixed indices in A.
    @param fixvalA the values of the fixed indices in A.
    @param clear_A_after_use eliminate storage used by A as soon as possible.
    @param B the B array.
    @param c_extB the indices of B that are shared with C.
    @param fixextB the fixed indices of B that are shared with C.
    @param c_intB the internal indices (shared between A and B).
    @param fixB the fixed indices in B.
    @param fixvalB the values of the fixed indices in B.
    @param clear_B_after_use eliminate storage used by B as soon as possible.
    @param ABfactor a multiplicative factor.
    @param C_is_zero_on_entry if true, contract will assume that C is initalize zero.
 */
template <int NC, int NA, int NB>
inline void
contract(
    Array<NC> &C, const IndexList &c_extCA, const IndexList &c_extCB,
    const IndexList &fixextCA, const IndexList &fixextCB,
    const IndexList &fixC, const BlockInfo<NC> &fixvalC,
    Array<NA> &A, const IndexList &c_extA, const IndexList &fixextA,
    const IndexList &c_intA,
    const IndexList &fixA, const BlockInfo<NA> &fixvalA,
    bool clear_A_after_use,
    Array<NB> &B, const IndexList &c_extB, const IndexList &fixextB,
    const IndexList &c_intB,
    const IndexList &fixB, const BlockInfo<NB> &fixvalB,
    bool clear_B_after_use,
    double ABfactor,
    bool C_is_zero_on_entry = false,
    sc::Ref<sc::RegionTimer> timer = 0)
{
  // Some of the arguments are copied to local variables.  This
  // is to allow modification of those arguments locally to permit
  // optimization of the cost of repacking the arrays.
  IndexList extCA(c_extCA), extCB(c_extCB);
  IndexList extA(c_extA), extB(c_extB);
  IndexList intA(c_intA), intB(c_intB);

  if (C.n_element_allocated() == 0
      || A.n_element_allocated() == 0
      || B.n_element_allocated() == 0) {
      return;
    }

  for (int i=0; i<fixC.n(); i++)
    if (fixC.i(i) >= fixC.n())
      throw std::invalid_argument("contract: C's fixed indices must be first");

  // Consistency checks
  if (extA.n() + extB.n() != NC - fixC.n() + fixextCA.n() + fixextCB.n()) {
      throw std::invalid_argument("contract: Number of externals on A + B != C");
    }
  if (extCA.n() + extCB.n() != NC - fixC.n() + fixextCA.n() + fixextCB.n()) {
      throw std::invalid_argument("contract: Number of indices on C inconsistent");
    }
  if (intA.n() != intB.n()) {
      throw std::invalid_argument(
          "contract: Number of internals on A and B inconsistent");
    }
  if (intA.n() + extA.n() != NA - fixA.n() + fixextA.n()) {
      throw std::invalid_argument("contract: Number of indices on A inconsistent");
    }
  if (intB.n() + extB.n() != NB - fixB.n() + fixextB.n()) {
      throw std::invalid_argument("contract: Number of indices on B inconsistent");
    }
  for (int i=0; i<extA.n(); i++) {
      if (A.index(extA.i(i)) != C.index(extCA.i(i))) {
          throw std::invalid_argument("contract: Range conflict between A and C");
        }
    }
  for (int i=0; i<extB.n(); i++) {
      if (B.index(extB.i(i)) != C.index(extCB.i(i))) {
          throw std::invalid_argument("contract: Range conflict between B and C");
        }
    }
  for (int i=0; i<intA.n(); i++) {
      if (A.index(intA.i(i)) != B.index(intB.i(i))) {
          throw std::invalid_argument("contract: Range conflict between A and B");
        }
    }

  RepackScheme<NC,NA,NB> repack_scheme(&A, extA, intA, fixA, &fixvalA,
                                       &B, intB, extB, fixB, &fixvalB,
                                       &C, extCA,extCB,fixC, &fixvalC,
                                       C_is_zero_on_entry?1:2);

//   std::cout << "trying to repack: original cost = "
//             << repack_scheme.cost()
//             << ", rpk A = " << repack_scheme.need_repack_A()
//             << " (" << A.n_element_allocated() << ")"
//             << ", rpk B = " << repack_scheme.need_repack_B()
//             << " (" << B.n_element_allocated() << ")"
//             << ", rpk C = " << repack_scheme.need_repack_C()
//             << " (" << C.n_element_allocated() << ")"
//             << std::endl;

  if (repack_scheme.cost() != 0.0) {
      RepackScheme<NC,NA,NB> tmp_repack_scheme(repack_scheme);
      for (int driver1 = 0; driver1 < 3; driver1++) {
          for (int driver2 = 0; driver2 < 3; driver2++) {
              if (driver1 == driver2) continue;
              tmp_repack_scheme.driver(driver1);
              tmp_repack_scheme.driver(driver2);

//               std::cout << " repack scheme: "
//                         << " d1 = " << driver1
//                         << " d2 = " << driver2
//                         << ", cost = "
//                         << tmp_repack_scheme.cost()
//                         << ", rpk A = "
//                         << tmp_repack_scheme.need_repack_A()
//                         << " (" << A.n_element_allocated() << ")"
//                         << ", rpk B = "
//                         << tmp_repack_scheme.need_repack_B()
//                         << " (" << B.n_element_allocated() << ")"
//                         << ", rpk C = "
//                         << tmp_repack_scheme.need_repack_C()
//                         << " (" << C.n_element_allocated() << ")"
//                         << std::endl;

              if (tmp_repack_scheme.cost() < repack_scheme.cost()) {
                  repack_scheme = tmp_repack_scheme;
                  repack_scheme.assign_indices(extA, intA,
                                               intB, extB,
                                               extCA,extCB);
                  //std::cout << "found a more efficient scheme" << std::endl;
                }
              if (repack_scheme.cost() == 0.0) break;
            }
          if (repack_scheme.cost() == 0.0) break;
        }
    }

  // Remap the blocks of A so that it is sorted by its external indices
  // and any fixed indices.
  if (timer) timer->enter("remap A");
  IndexList cmpAlist(extA,fixA);
  IndexListLess<NA> cmpA(cmpAlist);
  typename Array<NA>::cached_blockmap_t remappedAbm_local(cmpA);
  typename Array<NA>::cached_blockmap_t *remappedAbm_ptr;
  if (fixA.n() == 0 && A.use_blockmap_cache()) {
      remappedAbm_ptr = &A.blockmap_cache_entry(cmpAlist);
    }
  else {
      remap(remappedAbm_local, A, fixA, fixvalA);
      remappedAbm_ptr = &remappedAbm_local;
    }
  typename Array<NA>::cached_blockmap_t &remappedAbm=*remappedAbm_ptr;
  if (timer) timer->exit();
//    std::cout << "A:" << std::endl << A;
//    std::cout << "remappedA:" << std::endl << remappedA;

  // Repack the data of A, B, and C so DGEMM can be used
  // note: if need_repack is true then the matrix has not been transposed
  if (timer) timer->enter("repack1");
  if (repack_scheme.need_repack_A()) repack(A, extA, intA, fixA, fixvalA);
  if (repack_scheme.need_repack_B()) repack(B, intB, extB, fixB, fixvalB);
  if (repack_scheme.need_repack_C() && !C_is_zero_on_entry) {
      repack(C, extCA, extCB, fixC, fixvalC);
    }
  if (timer) timer->exit();

//   std::cout << "tA:" << transpose_A
//             << " tB:" << transpose_B
//             << " rA:" << need_repack_A
//             << " rB:" << need_repack_B
//             << " rC:" << need_repack_C
//             << std::endl;

#if 0
  // Fixed indices imply that a loop over those indices is been done
  // external to this routine.  If one array doesn't have all of the fixed
  // indices that other arrays have, then repacking that array will result
  // in extra work.  The code below detects this case.
  std::set<int> fixed_all, fixed_A, fixed_B, fixed_C;
  for (int i=0; i<fixA.n(); i++) {
      fixed_all.insert(fixvalA.block(i));
      fixed_A.insert(fixvalA.block(i));
    }
  for (int i=0; i<fixB.n(); i++) {
      fixed_all.insert(fixvalB.block(i));
      fixed_B.insert(fixvalB.block(i));
    }
  for (int i=0; i<fixC.n(); i++) {
      fixed_all.insert(fixvalC.block(i));
      fixed_C.insert(fixvalC.block(i));
    }
  if (repack_scheme.need_repack_A()
      && fixed_A.size() > 0
      && fixed_A != fixed_all) {
      std::cout << "PERFORMANCE WARNING: contract needed to repack A"
                << " but B and/or C have different fixed indices"
                << std::endl;
      throw std::runtime_error("contract: performance exception");
    }
  if (repack_scheme.need_repack_B()
      && fixed_B.size() > 0
      && fixed_B != fixed_all) {
      std::cout << "PERFORMANCE WARNING: contract needed to repack B"
                << " but A and/or C have different fixed indices"
                << std::endl;
      throw std::runtime_error("contract: performance exception");
    }
  if (repack_scheme.need_repack_C()
      && fixed_C.size() > 0
      && fixed_C != fixed_all) {
      std::cout << "PERFORMANCE WARNING: contract needed to repack C"
                << " but A and/or B have different fixed indices"
                << std::endl;
      throw std::runtime_error("contract: performance exception");
    }
#endif

  const typename Array<NB>::blockmap_t &
      Bbm = B.blockmap();
#ifdef USE_HASH
  const typename Array<NB>::blockhash_t &
      Bbh = B.blockhash();
#endif
  const typename Array<NC>::blockmap_t &
      Cbm = C.blockmap();

  BlockInfo<NA> Abi;
  Abi.assign_blocks(fixA,fixvalA);

  BlockInfo<NB> Bbi;
  Bbi.zero();
  Bbi.assign_blocks(fixB,fixvalB);

#ifndef USE_HASH
  typename Array<NB>::blockmap_t::const_iterator B_fixed_hint
      = Bbm.lower_bound(Bbi);
#endif

  typename Array<NC>::blockmap_t::const_iterator C_begin, C_end;
  BlockInfo<NC> Cbi_lb;
  BlockInfo<NC> Cbi_ub;
  for (int i=0; i<NC; i++) {
      Cbi_lb.block(i) = 0;
      Cbi_ub.block(i) = C.index(i).nblock();
    }
  Cbi_lb.assign_blocks(fixC,fixvalC);
  Cbi_ub.assign_blocks(fixC,fixvalC);
  C_begin = Cbm.lower_bound(Cbi_lb);
  C_end = Cbm.upper_bound(Cbi_ub);

  if (timer) timer->enter("C loop");
  for (typename Array<NC>::blockmap_t::const_iterator
           Citer = C_begin;
       Citer != C_end;
       Citer++) {
      const BlockInfo<NC> &Cbi = Citer->first;
      double *Cdata = Citer->second;
      Abi.assign_blocks(extA, Cbi, extCA);
      std::pair<
          typename Array<NA>::cached_blockmap_t::const_iterator,
          typename Array<NA>::cached_blockmap_t::const_iterator >
          rangeA;
#ifdef USE_BOUND
      // cannot use equal range on remappedA because bound is used to sort
      // rangeA = remappedAbm.equal_range(Abi);
      Abi.set_bound(DBL_MAX);
      rangeA.first = remappedAbm.lower_bound(Abi);
      Abi.set_bound(0.0);
      rangeA.second = remappedAbm.upper_bound(Abi);
#else
      rangeA = remappedAbm.equal_range(Abi);
#endif
      typename Array<NA>::cached_blockmap_t::const_iterator
          firstA = rangeA.first,
          fenceA = rangeA.second;
      Bbi.assign_blocks(extB, Cbi, extCB);
      blasint n_extB = Cbi.subset_size(C.indices(), extCB);
      blasint n_extA = Cbi.subset_size(C.indices(), extCA);
#ifdef USE_HASH
      typename Array<NB>::blockhash_t::const_iterator Biter;
#else
      typename Array<NB>::blockmap_t::const_iterator Biter = Bbm.begin();
#endif
      if (timer) timer->enter("A loop");
      for (typename Array<NA>::cached_blockmap_t::const_iterator
               Aiter = firstA;
           Aiter != fenceA;
           Aiter++) {
          const BlockInfo<NA> &Abi = Aiter->first;
          double *Adata = Aiter->second;
          Bbi.assign_blocks(intB, Abi, intA);
#ifdef USE_HASH
          Biter = Bbh.find(Bbi);
          if (Biter == Bbh.end()) continue;
#else
          if (fixB.n() > 0) {
#if USE_STL_MULTIMAP
              Biter = Bbm.find(Bbi);
#else
              Biter = Bbm.find(B_fixed_hint, Bbi);
#endif
            }
          else {
              //blindly using a hint here makes this a bit slower
              //Biter = Bbm.find(Biter, Bbi);
              Biter = Bbm.find(Bbi);
            }
          if (Biter == Bbm.end()) continue;
#endif
          double *Bdata = Biter->second;
          blasint n_int = Abi.subset_size(A.indices(), intA);

          double one = 1.0;
          if (timer) timer->enter("dgemm");

          double t0 = cpu_walltime();

          if (n_extA == 1 && n_int == 1) {
              double tmp = ABfactor * Adata[0];
              for (int i=0; i<n_extB; i++) {
                  Cdata[i] += tmp*Bdata[i];
                }
            }
          else if (n_extA == 1 && n_extB == 1) {
              double tmp = 0.0;
              for (int i=0; i<n_int; i++) {
                  tmp += Adata[i]*Bdata[i];
                }
              Cdata[0] += ABfactor*tmp;
            }
          else if (n_int == 1 && n_extB == 1) {
              double tmp = ABfactor*Bdata[0];
              for (int i=0; i<n_extA; i++) {
                  Cdata[i] += Adata[i]*tmp;
                }
            }
          else if (n_int == 1) {
              if (repack_scheme.transpose_C()) {
                  for (int i=0,ij=0; i<n_extB; i++) {
                      for (int j=0; j<n_extA; j++,ij++) {
                          Cdata[ij] += ABfactor*Adata[j]*Bdata[i];
                        }
                    }
                }
              else {
                  for (int i=0,ij=0; i<n_extA; i++) {
                      for (int j=0; j<n_extB; j++,ij++) {
                          Cdata[ij] += ABfactor*Adata[i]*Bdata[j];
                        }
                    }
                }
            }
          else if (n_extA == 1) {
              if (repack_scheme.transpose_B()) {
                  for (int i=0,ij=0; i<n_extB; i++) {
                      double tmp = 0.0;
                      for (int j=0; j<n_int; j++,ij++) {
                          tmp += Adata[j]*Bdata[ij];
                        }
                      Cdata[i] += tmp * ABfactor;
                    }
                }
              else {
                  for (int i=0; i<n_extB; i++) {
                      double tmp = 0.0;
                      for (int j=0,ij=i; j<n_int; j++,ij+=n_extB) {
                          tmp += Adata[j]*Bdata[ij];
                        }
                      Cdata[i] += tmp * ABfactor;
                    }
                }
            }
          else if (n_extB == 1) {
              if (repack_scheme.transpose_A()) {
                  for (int i=0; i<n_extA; i++) {
                      double tmp = 0.0;
                      for (int j=0,ij=i; j<n_int; j++,ij+=n_extA) {
                          tmp += Bdata[j]*Adata[ij];
                        }
                      Cdata[i] += tmp * ABfactor;
                    }
                }
              else {
                  for (int i=0,ij=0; i<n_extA; i++) {
                      double tmp = 0.0;
                      for (int j=0; j<n_int; j++,ij++) {
                          tmp += Bdata[j]*Adata[ij];
                        }
                      Cdata[i] += tmp * ABfactor;
                    }
                }
            }
          else if (repack_scheme.transpose_C()) {
              const char *tA = "T";
              blasint lda = n_int;
              if (repack_scheme.transpose_A()) { tA = "N"; lda = n_extA; }

              const char *tB = "T";
              blasint ldb = n_extB;
              if (repack_scheme.transpose_B()) { tB = "N"; ldb = n_int; }

              blasint ldc = n_extA;

//               std::cout << " tA: " << tA
//                         << " tB: " << tB
//                         << " nr: " << n_extA
//                         << " nc: " << n_extB
//                         << " nl: " << n_int
//                         << " lda: " << lda
//                         << " ldb: " << ldb
//                         << " ldc: " << ldc
//                         << std::endl;

              F77_DGEMM(tA, tB, &n_extA, &n_extB, &n_int,
                        &ABfactor,Adata,&lda,Bdata,&ldb,
                        &one,Cdata,&ldc);
            }
          else {
              const char *tA = "N";
              blasint lda = n_int;
              if (repack_scheme.transpose_A()) { tA = "T"; lda = n_extA; }

              const char *tB = "N";
              blasint ldb = n_extB;
              if (repack_scheme.transpose_B()) { tB = "T"; ldb = n_int; }

              blasint ldc = n_extB;

              F77_DGEMM(tB, tA, &n_extB, &n_extA, &n_int,
                        &ABfactor,Bdata,&ldb,Adata,&lda,
                        &one,Cdata,&ldc);
            }
#ifdef USE_COUNT_DGEMM
          count_dgemm(n_extA, n_int, n_extB,
                      cpu_walltime()-t0);
#endif
          if (timer) timer->exit();
        }
      if (timer) timer->exit();
    }
  if (timer) timer->exit();
  
  // Repack the data of A, B, and C to the orginal data layout
  if (timer) timer->enter("repack2");
  if (clear_A_after_use) A.clear();
  else {
      if (repack_scheme.need_repack_A()) {
          repack(A, extA, intA, fixA, fixvalA, true);
        }
    }

  if (clear_B_after_use) B.clear();
  else {
      if (repack_scheme.need_repack_B()) {
          repack(B, intB, extB, fixB, fixvalB, true);
        }
    }

  if (repack_scheme.need_repack_C()) {
      repack(C, extCA, extCB, fixC, fixvalC, true);
    }
  if (timer) timer->exit();

  if (timer) timer->enter("bounds");
  C.compute_bounds();
  if (timer) timer->exit();
}

/// Contract two arrays to produce a scalar.
template <int N>
inline double
scalar_contract(
    Array<N> &c,
    Array<N> &a, const IndexList &alist)
{
  // Consistency checks
  if (alist.n() != N) {
      throw std::invalid_argument(
          "sma::scalar_contract: # of indices inconsistent");
    }
  for (int i=0; i<N; i++) {
      if (c.index(i) != a.index(alist.i(i)))
          throw std::invalid_argument(
              "sma::scalar_contract: indices don't agree");
    }

  bool same_index_order = alist.is_identity();

  double r = 0.0;
  const typename Array<N>::blockmap_t &amap = a.blockmap();
  const typename Array<N>::blockmap_t &cmap = c.blockmap();
  IndexList clist = alist.reverse_mapping();
  bool use_hint;
  if (clist.i(0) == 0) use_hint = true;
  else use_hint = false;
  typename Array<N>::blockmap_t::const_iterator citer = cmap.begin();
  for (typename Array<N>::blockmap_t::const_iterator aiter = amap.begin();
       aiter != amap.end();
       aiter++) {
      BlockInfo<N> cbi(aiter->first,clist);
#if USE_STL_MULTIMAP
      citer = cmap.find(cbi);
#else
      if (use_hint) citer = cmap.find(citer,cbi);
      else citer = cmap.find(cbi);
#endif
      if (citer == cmap.end()) continue;
      double *cdata = citer->second;
      double *adata = aiter->second;
      if (same_index_order) {
          int sz = c.block_size(cbi);
          for (int i=0; i<sz; i++) r += cdata[i] * adata[i];
        }
      else {
          BlockIter<N> cbiter(c.indices(),cbi);
          int coff = 0;
          for (cbiter.start(); cbiter.ready(); cbiter++,coff++) {
              r += cdata[coff] * adata[cbiter.subset_offset(alist)];
            }
        }
    }
  
  return r;
}

/// Determines the blocks that will be need to store a contraction result.
/// This will compute the blocks in an array, C, that will
/// be needed if the contraction A = C * B or B = C * A or
/// C = A * B is performed. The blocks in A and B must have
/// already been computed, and C is not cleared before its
/// blocks are allocated.
template <int NC, int NA, int NB>
inline void
contract_union(
    Array<NC> &C, const IndexList &extCA, const IndexList &extCB,
    const IndexList &fixC, const BlockInfo<NC> &fixvalC,
    Array<NA> &A, const IndexList &extA, const IndexList &intA,
    const IndexList &fixA, const BlockInfo<NA> &fixvalA,
    Array<NB> &B, const IndexList &extB, const IndexList &intB,
    const IndexList &fixB, const BlockInfo<NB> &fixvalB)
{
  // Consistency checks
  if (extA.n() + extB.n() != NC - fixC.n()) {
      std::cerr << "NA = " << NA << std::endl;
      std::cerr << "intA = " << intA << std::endl;
      std::cerr << "extA = " << extA << std::endl;
      std::cerr << "fixA = " << fixA << std::endl;
      std::cerr << "NB = " << NB << std::endl;
      std::cerr << "intB = " << intB << std::endl;
      std::cerr << "extB = " << extB << std::endl;
      std::cerr << "fixB = " << fixB << std::endl;
      std::cerr << "NC = " << NC << std::endl;
      std::cerr << "extCA = " << extCA << std::endl;
      std::cerr << "extCB = " << extCB << std::endl;
      std::cerr << "fixC = " << fixC << std::endl;
      throw std::invalid_argument("contract_union: Number of externals on A + B != C");
    }
  if (extCA.n() + extCB.n() != NC - fixC.n()) {
      throw std::invalid_argument("contract_union: Number of indices on C inconsistent");
    }
  if (intA.n() != intB.n()) {
      throw std::invalid_argument(
          "contract_union: Number of internals on A and B inconsistent");
    }
  if (intA.n() + extA.n() != NA - fixA.n()) {
      std::cerr << "NA = " << NA << std::endl;
      std::cerr << "extA.n() = " << extA.n() << std::endl;
      std::cerr << "fixA.n() = " << fixA.n() << std::endl;
      std::cerr << "intA.n() = " << intA.n() << " (";
      for (int i=0; i<intA.n(); i++) {
        std::cerr << " " << intA.i(i);
      }
      std::cerr << ")" << std::endl;
      throw std::invalid_argument("contract_union: Number of indices on A inconsistent");
    }
  if (intB.n() + extB.n() != NB - fixB.n()) {
      throw std::invalid_argument("contract_union: Number of indices on B inconsistent");
    }
  for (int i=0; i<extB.n(); i++) {
      if (B.index(extB.i(i)) != C.index(extCB.i(i))) {
          throw std::invalid_argument("contract_union: Range conflict between B and C");
        }
    }
  for (int i=0; i<intA.n(); i++) {
      if (A.index(intA.i(i)) != B.index(intB.i(i))) {
          throw std::invalid_argument("contract_union: Range conflict between A and B");
        }
    }


  // Remap the blocks of A so that it is sorted by its internal indices and
  // any fixed indices.  The fixed indices appear first (are most
  // significant wrt the ordering) so we can get the iterator bounds for
  // relevant internal indices more easily.  Data is not moved.
  IndexList cmpAlist(fixA, intA);
  IndexListLess<NA> cmpA(cmpAlist);
  typename Array<NA>::cached_blockmap_t remappedAbm_local(cmpA);
  typename Array<NA>::cached_blockmap_t *remappedAbm_ptr;
  if (fixA.n() == 0 && A.use_blockmap_cache()) {
      remappedAbm_ptr = &A.blockmap_cache_entry(cmpAlist);
    }
  else {
      remap(remappedAbm_local, A, fixA, fixvalA);
      remappedAbm_ptr = &remappedAbm_local;
    }
  typename Array<NA>::cached_blockmap_t &remappedAbm=*remappedAbm_ptr;

  // Remap the blocks of B so that it is sorted by its internal indices
  // and any fixed indices.  Data is not moved.
  IndexList cmpBlist(intB, fixB);
  IndexListLess<NB> cmpB(cmpBlist);
  typename Array<NB>::cached_blockmap_t remappedBbm_local(cmpB);
  typename Array<NB>::cached_blockmap_t *remappedBbm_ptr;
  if (fixB.n() == 0 && B.use_blockmap_cache()) {
      remappedBbm_ptr = &B.blockmap_cache_entry(cmpBlist);
    }
  else {
      remap(remappedBbm_local, B, fixB, fixvalB);
      remappedBbm_ptr = &remappedBbm_local;
    }
  typename Array<NB>::cached_blockmap_t &remappedBbm=*remappedBbm_ptr;
      

//    std::cout << "beginning loops" << std::endl;
//    std::cout << "extA = " << extA << std::endl;
//    std::cout << "extB = " << extB << std::endl;
//    std::cout << "extCA = " << extCA << std::endl;
//    std::cout << "extCB = " << extCB << std::endl;

  BlockInfo<NA> ablockinfo;
  for (int i=0; i<NA; i++) ablockinfo.block(i) = 0;
#ifdef USE_BOUND
  ablockinfo.set_bound(DBL_MAX);
#endif
  ablockinfo.assign_blocks(fixA, fixvalA);
  typename Array<NA>::cached_blockmap_t::const_iterator abegin;
  abegin = remappedAbm.lower_bound(ablockinfo);

  BlockInfo<NB> bblockinfo;
  bblockinfo.assign_blocks(fixB, fixvalB);

  BlockInfo<NC> cblockinfo;
  cblockinfo.assign_blocks(fixC, fixvalC);

  while (abegin != remappedAbm.end()) {
      ablockinfo = abegin->first;

      // if there are fixed indices, then abegin might be beyond
      // the fixed indices that we are interested in
      if (!ablockinfo.equiv_blocks(fixA, fixvalA)) break;

#ifdef USE_BOUND
#if 0 && USE_BOUNDS_IN_CONTRACT_UNION
      if (B.bound() < DBL_EPSILON) {
          ablockinfo.set_bound(C.tolerance()/DBL_EPSILON);
        }
      else {
          ablockinfo.set_bound(C.tolerance()/B.bound());
        }
#else
      ablockinfo.set_bound(0.0);
#endif
#endif
      typename Array<NA>::cached_blockmap_t::const_iterator
          afence = remappedAbm.upper_bound(ablockinfo);
      bblockinfo.assign_blocks(intB,ablockinfo,intA);
      std::pair<typename Array<NB>::cached_blockmap_t::const_iterator,
          typename Array<NB>::cached_blockmap_t::const_iterator>
          brange;
      // cannot use equal_range on remappedB since bounds are used to sort
      // brange = remappedBbm.equal_range(bblockinfo);
#ifdef USE_BOUND
      bblockinfo.set_bound(DBL_MAX);
#endif
      brange.first = remappedBbm.lower_bound(bblockinfo);
#ifdef USE_BOUND
#if 0 && USE_BOUNDS_IN_CONTRACT_UNION
      if (A.bound() < DBL_EPSILON) {
          bblockinfo.set_bound(C.tolerance()/DBL_EPSILON);
        }
      else {
          bblockinfo.set_bound(C.tolerance()/A.bound());
        }
#else
      bblockinfo.set_bound(0.0);
#endif
#endif
      brange.second = remappedBbm.upper_bound(bblockinfo);
      typename Array<NB>::cached_blockmap_t::const_iterator
          bbegin = brange.first,
          bfence = brange.second;
//        std::cout << "  in internal loop" << std::endl;
      for (typename Array<NA>::cached_blockmap_t::const_iterator
               aiter = abegin;
           aiter != afence;
           aiter++) {
#ifdef USE_BOUND
          double a_block_bound = aiter->first.bound();
#endif
          cblockinfo.assign_blocks(extCA,aiter->first,extA);
//            std::cout << "  A blocks: " << aiter->first << std::endl;
          for (typename Array<NB>::cached_blockmap_t::const_iterator
                   biter = bbegin;
               biter != bfence;
               biter++) {
#ifdef USE_BOUND
#if 0 && USE_BOUNDS_IN_CONTRACT_UNION
              if (a_block_bound * biter->first.bound() < C.tolerance()) {
                  continue;
                }
#endif
#endif
              cblockinfo.assign_blocks(extCB,biter->first,extB);
//                std::cout << "    B blocks: " << biter->first
//                          << " adding " << cblockinfo << std::endl;
              C.add_unallocated_block(cblockinfo);
            }
        }
#ifdef USE_BOUND
      ablockinfo.set_bound(0.0);
#endif
      abegin = remappedAbm.upper_bound(ablockinfo);
    }
}

}

}

#endif

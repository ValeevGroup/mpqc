//
// mtensor.h
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

#ifndef _ccr12_mtensor_h
#define _ccr12_mtensor_h

#include <numeric>
#include <chemistry/qc/ccr12/tensor.h>
#include <chemistry/qc/ccr12/ccr12_info.h>
#include <math/scmat/matrix.h>

namespace {
  void print_tile(double* data, int n01, int n23) {
    using namespace sc;
    RefSCMatrix mat = SCMatrixKit::default_matrixkit()->matrix(new SCDimension(n01),
                                                               new SCDimension(n23));
    mat.assign(data);
    mat.print("tile");
  }

}

namespace sc {

    /// return intersect of two ranges defined as pair<start,fence>, i.e. [start, fence)
    template <typename I> inline std::pair<I,I> intersect(const std::pair<I,I>& range1, const std::pair<I,I>& range2) {
      I start_max = std::max(range1.first, range2.first);
      I fence_min = std::min(range1.second, range2.second);
      if (start_max < fence_min) return make_pair(start_max, fence_min);
      return make_pair(I(0),I(0));
    }

    /// return true if r is contained in range defined as pair<start,fence>, i.e. [start, fence)
    template <typename I> inline bool in(const std::pair<I,I>& r, const std::pair<I,I>& range) {
      return r.first >= range.first && r.second <= range.second;
    }
    /// return true if i is in range defined as pair<start,fence>, i.e. [start, fence)
    template <typename I> inline bool in(I i, const std::pair<I,I>& range) {
      return i >= range.first && i < range.second;
    }

//////////

  /// Tensor metadata is implicit; MTensor is Tensor + metadata
  template <size_t NDIM>
  class MTensor {
    public:
      typedef CCR12_Info Info;
      typedef long tile_index;
      typedef long element_index;
      typedef std::pair<tile_index, tile_index> tile_range;  //< [start, fence)
      typedef std::vector<tile_range> tile_ranges;
      typedef std::pair<element_index, element_index> element_range;  //< [start, fence)
      typedef std::vector<element_range> element_ranges;
      typedef std::vector<element_index> element_index_map;
      static const size_t ndim = NDIM;

      /**
       * \param range is a vector of NDIM elements that contains tile range for each dimension [start, fence)
       */
      MTensor(Info const* info,
              Tensor* tensor,
              const tile_ranges& range) :
                tensor_(tensor), info_(info), range_(range)
      {
        MPQC_ASSERT(range.size() == NDIM);
      }

      /** copies contents of src to this.
       *  indexing of elements in src may be different from that in this.
       *  arrays eimap0 through eimap3 provide mapping
       *
       *  src integrals are assumed to not be antisymmetrized
       *
       *  mapped_element_ranges specifies the elements contained in src. Tiles that request elements
       *  outside these ranges will be skipped. Default to allow all elements.
       *
       *  \param src_is_2301 set to true if src contains <2 3|0 1> integrals
       */
      void convert(const Ref<DistArray4>& src, unsigned int tbtype,
                   const element_index_map& eimap0,
                   const element_index_map& eimap1,
                   const element_index_map& eimap2,
                   const element_index_map& eimap3,
                   element_ranges const* mapped_element_ranges = 0,
                   bool src_is_2301 = false);

      /** copies contents of src to this. src is assumed to be a RefSCMatrix or a RefSymmSCMatrix
       *  indexing of elements in src may be different from that in this.
       *  arrays eimap0 through eimap3 provide mapping.
       *
       *  src integrals are assumed to not be antisymmetrized.
       *
       *  to compute antisymmetrized integrals src must contain (01|23) and (01|32) or (10|23).
       *  if it doesn't, use convert() function that takes 2 src tensors.
       */
      template <typename SrcType> void convert(const SrcType& src,
                                               int n1,
                                               int n3,
                                               const element_index_map& eimap0,
                                               const element_index_map& eimap1,
                                               const element_index_map& eimap2,
                                               const element_index_map& eimap3,
                                               element_ranges const* mapped_element_ranges = 0) {

        MPQC_ASSERT(NDIM == 4);
        MPQC_ASSERT(mapped_element_ranges != 0);

        // determine which tiles map to the allowed element ranges in src
        const std::vector<bool> tile0_in_erange = tiles_in_erange(range_[0], eimap0, (*mapped_element_ranges)[0]);
        const std::vector<bool> tile1_in_erange = tiles_in_erange(range_[1], eimap1, (*mapped_element_ranges)[1]);
        const std::vector<bool> tile2_in_erange = tiles_in_erange(range_[2], eimap2, (*mapped_element_ranges)[2]);
        const std::vector<bool> tile3_in_erange = tiles_in_erange(range_[3], eimap3, (*mapped_element_ranges)[3]);

        // split work over tasks assuming that all tasks can access src
        const int nproc_with_ints = info_->mem()->n();
        const int me = info_->mem()->me();

        const size_t ntiles0 = range_[0].second - range_[0].first;
        const size_t ntiles1 = range_[1].second - range_[1].first;
        const size_t ntiles2 = range_[2].second - range_[2].first;
        const size_t ntiles3 = range_[3].second - range_[3].first;

        // if espace0 is equivalent to espace1, or espace2 equivalent to espace3,
        // can compute antisymmetric integrals using only (espace0 espace1| espace2 espace3)
        // else also need (espace0 espace1| espace3 espace2), or an equivalent
        const bool espace0_eq_espace1 = ((*mapped_element_ranges)[0] == (*mapped_element_ranges)[1]);
        const bool espace2_eq_espace3 = ((*mapped_element_ranges)[2] == (*mapped_element_ranges)[3]);
        const bool can_always_antisymmetrize = (espace0_eq_espace1 || espace2_eq_espace3);
        MPQC_ASSERT(can_always_antisymmetrize == true);  // this is likely to break the logic downstream
                                                    // I clearly don't understand enough at the moment

        // TODO check if equivalence of MTensor spaces is matched by their relationship in src

        // assume that src is accessible from all nodes
        if (1) {

          // how do I determine max size of the tiles?
          const size_t maxsize1 = info_->maxtilesize();
          const size_t maxtilesize = maxsize1 * maxsize1 * maxsize1 * maxsize1;
          double* data = info_->mem()->malloc_local_double(maxtilesize);

          for (long t0 = range_[0].first; t0 < range_[0].second; ++t0) {
            const long size0 = info_->get_range(t0);
            const long offset0 = info_->get_offset(t0);
            const long spin0 = info_->get_spin(t0);
            const bool in_erange0 = tile0_in_erange[t0];

            for (long t1 = std::max(range_[1].first, t0); t1 < range_[1].second; ++t1) {
              const long size1 = info_->get_range(t1);
              const long offset1 = info_->get_offset(t1);
              const long spin1 = info_->get_spin(t1);
              const bool in_erange1 = tile1_in_erange[t1];

              const bool aaaa = (spin0 == spin1);

              for (long t2 = range_[2].first; t2 < range_[2].second; ++t2) {
                const long size2 = info_->get_range(t2);
                const long offset2 = info_->get_offset(t2);
                const long spin2 = info_->get_spin(t2);
                const bool in_erange2 = tile2_in_erange[t2];

                const bool abab = ((spin0 != spin1) && (spin0 == spin2));
                const bool abba = ((spin0 != spin1) && (spin0 != spin2));

                for (long t3 = std::max(range_[3].first, t2); t3 < range_[3].second; ++t3) {
                  const long size3 = info_->get_range(t3);
                  const long offset3 = info_->get_offset(t3);
                  const long spin3 = info_->get_spin(t3);
                  const bool in_erange3 = tile3_in_erange[t3];

                  // cartesian ordinal is used as a key for hashing tiles
                  const long tile_key = (t3-range_[3].first) +
                                        ntiles3 * ( (t2-range_[2].first) +
                                                    ntiles2 * ( (t1-range_[1].first) +
                                                                ntiles1 * (t0-range_[0].first)
                                                              )
                                                  );

                  // since all procs can access integrals, use locality
                  if (tensor_->exists(tile_key) && tensor_->is_this_local(tile_key)) {

                    if ((!in_erange0 || !in_erange1 || !in_erange2 || !in_erange3))
                        continue;

                    // if erange[0] == erange[1], clearly can antisymmetrize 0 and 1
                    // if erange[2] == erange[3], clearly can antisymmetrize 2 and 3
                    bool antisymmetrize01 = espace0_eq_espace1;
                    bool antisymmetrize23 = espace2_eq_espace3;
                    // prefer to antisymmetrize 23 since can do that with 1 block of DistArray4
                    if (antisymmetrize23) antisymmetrize01 = false;

                    long size = size0 * size1 * size2 * size3;
                    std::fill(data, data+size, 0.0);

                      long i0123 = 0;
                      for (int i0 = 0; i0 < size0; ++i0) {
                        const int ii0 = eimap0[offset0 + i0];
                        for (int i1 = 0; i1 < size1; ++i1) {
                          const int ii1 = eimap1[offset1 + i1];

                          //if (ii0 == ii1 && aaaa) continue;

                          // split the work in round robin fashion by processors who have access to integrals
                          long i01;
                          if (antisymmetrize01) {
                            const long iii0 = std::max(ii0, ii1);
                            const long iii1 = std::min(ii0, ii1);
                            i01 = iii0 * (iii0+1)/2 + iii1;
                          } else {
                            i01 = ii0 * n1 + ii1;
                          }
                          const long i01_proc = i01 % nproc_with_ints;
                          if (i01_proc != me)
                            continue;

                          if (antisymmetrize23) {

                            for (int i2 = 0; i2 < size2; ++i2) {
                              const int ii2 = eimap2[offset2 + i2];
                              for (int i3 = 0; i3 < size3; ++i3, ++i0123) {
                                const int ii3 = eimap3[offset3 + i3];

                                const int i23 = ii2 * n3 + ii3;
                                const int i32 = ii3 * n3 + ii2;

                                if (abab) {
                                  const double integral_0123 = src(i01,i23);
                                  data[i0123] = integral_0123;
                                } else if (abba) {
                                  const double integral_0132 = src(i01,i32);
                                  data[i0123] = -integral_0132;
                                } else { // aaaa
                                  const double integral_0123 = src(i01,i23);
                                  const double integral_0132 = src(i01,i32);
                                  data[i0123] = integral_0123 - integral_0132;
                                }

                              }
                            }

                          } // antisymmetrize23

                          else if (antisymmetrize01) { // means can't antisymmetrize 2 and 3

                            const int i10 = i1 * n1 + i0;

                            for (int i2 = 0; i2 < size2; ++i2) {
                              const int ii2 = eimap2[offset2 + i2];
                              for (int i3 = 0; i3 < size3; ++i3, ++i0123) {
                                const int ii3 = eimap3[offset3 + i3];

                                const int i23 = ii2 * n3 + ii3;

                                if (abab) {
                                  const double integral_0123 = src(i01,i23);
                                  data[i0123] = integral_0123;
                                } else if (abba) {
                                  const double integral_1023 = src(i10,i23);
                                  data[i0123] = -integral_1023;
                                } else { // aaaa
                                  const double integral_0123 = src(i01,i23);
                                  const double integral_1023 = src(i10,i23);
                                  data[i0123] = integral_0123 - integral_1023;
                                }

                              }
                            }

                          } // antisymmetrize01

                          else { // do not antisymmetrize
                            MPQC_ASSERT(false); // should not happen, but only SMITH knows :-)
                          }

                        }
                      }

                      tensor_->put_block(tile_key, data);

#if 0
                const double sum = std::accumulate(data, data+size, 0.0);
                ExEnv::out0() << "tiles = (" << t0 << "," << t1 << "," << t2 << "," << t3
                              << ")  key = " << tile_key  << " sum = " << sum << std::endl;

                for(int i=0; i<size; ++i) {
                  ExEnv::out0() << "data[" << i << "] = " << data[i] << std::endl;
                }
#endif
                  } // if this task will process this tile

                } // t3
              } // t2
            } // t1
          } // t0

          info_->mem()->free_local_double(data);

        } // if this task can access the integrals

      } // MTensor<4>::convert() from RefSCMatrix or RefSymmSCMatrix

      /** copies contents of (src0, src1) to this, where src0 = (01|23) and src1 = (01|32), if
       *  src1_is_1023 == false, or src1 = (10|23), if src1_is_1023 == true.
       * .src arrays are assumed to be RefSCMatrix or RefSymmSCMatrix
       *  indexing of elements in src arrays may be different from that in this.
       *  arrays eimap0 through eimap3 provide mapping.
       *
       *  src integrals are assumed to not be antisymmetrized.
       *
       */
      template <typename SrcType> void convert(const SrcType& src0, const SrcType& src1,
                                               int n0,
                                               int n1,
                                               int n2,
                                               int n3,
                                               const element_index_map& eimap0,
                                               const element_index_map& eimap1,
                                               const element_index_map& eimap2,
                                               const element_index_map& eimap3,
                                               element_ranges const* mapped_element_ranges = 0,
                                               bool src1_is_1023 = false) {

        MPQC_ASSERT(NDIM == 4);
        MPQC_ASSERT(mapped_element_ranges != 0);

        // determine which tiles map to the allowed element ranges in src
        const std::vector<bool> tile0_in_erange = tiles_in_erange(range_[0], eimap0, (*mapped_element_ranges)[0]);
        const std::vector<bool> tile1_in_erange = tiles_in_erange(range_[1], eimap1, (*mapped_element_ranges)[1]);
        const std::vector<bool> tile2_in_erange = tiles_in_erange(range_[2], eimap2, (*mapped_element_ranges)[2]);
        const std::vector<bool> tile3_in_erange = tiles_in_erange(range_[3], eimap3, (*mapped_element_ranges)[3]);

        // split work over tasks assuming that all tasks can access src
        const int nproc_with_ints = info_->mem()->n();
        const int me = info_->mem()->me();

        const size_t ntiles0 = range_[0].second - range_[0].first;
        const size_t ntiles1 = range_[1].second - range_[1].first;
        const size_t ntiles2 = range_[2].second - range_[2].first;
        const size_t ntiles3 = range_[3].second - range_[3].first;

        // if espace0 is equivalent to espace1, or espace2 equivalent to espace3,
        // can compute antisymmetric integrals using only (espace0 espace1| espace2 espace3)
        // else also need (espace0 espace1| espace3 espace2), or an equivalent
        const bool espace0_eq_espace1 = ((*mapped_element_ranges)[0] == (*mapped_element_ranges)[1]);
        const bool espace2_eq_espace3 = ((*mapped_element_ranges)[2] == (*mapped_element_ranges)[3]);
        const bool can_always_antisymmetrize = (espace0_eq_espace1 || espace2_eq_espace3);
        MPQC_ASSERT(can_always_antisymmetrize == true);  // this is likely to break the logic downstream
                                                    // I clearly don't understand enough at the moment

        // TODO check if equivalence of MTensor spaces is matched by their relationship in src

        // assume that src is accessible from all nodes
        if (1) {

          // how do I determine max size of the tiles?
          const size_t maxsize1 = info_->maxtilesize();
          const size_t maxtilesize = maxsize1 * maxsize1 * maxsize1 * maxsize1;
          double* data = info_->mem()->malloc_local_double(maxtilesize);

          for (long t0 = range_[0].first; t0 < range_[0].second; ++t0) {
            const long size0 = info_->get_range(t0);
            const long offset0 = info_->get_offset(t0);
            const long spin0 = info_->get_spin(t0);
            const bool in_erange0 = tile0_in_erange[t0];

            for (long t1 = std::max(range_[1].first, t0); t1 < range_[1].second; ++t1) {
              const long size1 = info_->get_range(t1);
              const long offset1 = info_->get_offset(t1);
              const long spin1 = info_->get_spin(t1);
              const bool in_erange1 = tile1_in_erange[t1];

              const bool aaaa = (spin0 == spin1);

              for (long t2 = range_[2].first; t2 < range_[2].second; ++t2) {
                const long size2 = info_->get_range(t2);
                const long offset2 = info_->get_offset(t2);
                const long spin2 = info_->get_spin(t2);
                const bool in_erange2 = tile2_in_erange[t2];

                const bool abab = ((spin0 != spin1) && (spin0 == spin2));
                const bool abba = ((spin0 != spin1) && (spin0 != spin2));

                for (long t3 = std::max(range_[3].first, t2); t3 < range_[3].second; ++t3) {
                  const long size3 = info_->get_range(t3);
                  const long offset3 = info_->get_offset(t3);
                  const long spin3 = info_->get_spin(t3);
                  const bool in_erange3 = tile3_in_erange[t3];

                  // cartesian ordinal is used as a key for hashing tiles
                  const long tile_key = (t3-range_[3].first) +
                                        ntiles3 * ( (t2-range_[2].first) +
                                                    ntiles2 * ( (t1-range_[1].first) +
                                                                ntiles1 * (t0-range_[0].first)
                                                              )
                                                  );

                  // since all procs can access integrals, use locality
                  if (tensor_->exists(tile_key) && tensor_->is_this_local(tile_key)) {

                    if ((!in_erange0 || !in_erange1 || !in_erange2 || !in_erange3))
                        continue;

                    // if erange[0] == erange[1], clearly can antisymmetrize 0 and 1
                    // if erange[2] == erange[3], clearly can antisymmetrize 2 and 3
                    bool antisymmetrize01 = espace0_eq_espace1;
                    bool antisymmetrize23 = espace2_eq_espace3;
                    // prefer to antisymmetrize 23 since can do that with 1 block of DistArray4
                    if (antisymmetrize23) antisymmetrize01 = false;

                    long size = size0 * size1 * size2 * size3;
                    std::fill(data, data+size, 0.0);

                      long i0123 = 0;
                      for (int i0 = 0; i0 < size0; ++i0) {
                        const int ii0 = eimap0[offset0 + i0];
                        for (int i1 = 0; i1 < size1; ++i1) {
                          const int ii1 = eimap1[offset1 + i1];

                          //if (ii0 == ii1 && aaaa) continue;

                          // split the work in round robin fashion by processors who have access to integrals
                          long i01;
                          if (antisymmetrize01) {
                            const long iii0 = std::max(ii0, ii1);
                            const long iii1 = std::min(ii0, ii1);
                            i01 = iii0 * (iii0+1)/2 + iii1;
                          } else {
                            i01 = ii0 * n1 + ii1;
                          }
                          const long i01_proc = i01 % nproc_with_ints;
                          if (i01_proc != me)
                            continue;

                          if (antisymmetrize23) {

                            for (int i2 = 0; i2 < size2; ++i2) {
                              const int ii2 = eimap2[offset2 + i2];
                              for (int i3 = 0; i3 < size3; ++i3, ++i0123) {
                                const int ii3 = eimap3[offset3 + i3];

                                const int i23 = ii2 * n3 + ii3;
                                const int i32 = ii3 * n2 + ii2;

                                if (abab) {
                                  const double integral_0123 = src0(i01,i23);
                                  data[i0123] = integral_0123;
                                } else if (abba) {
                                  const double integral_0132 = src1(i01,i32);
                                  data[i0123] = -integral_0132;
                                } else { // aaaa
                                  const double integral_0123 = src0(i01,i23);
                                  const double integral_0132 = src1(i01,i32);
                                  data[i0123] = integral_0123 - integral_0132;
                                }

                              }
                            }

                          } // antisymmetrize23

                          else if (antisymmetrize01) { // means can't antisymmetrize 2 and 3

                            const int i10 = i1 * n0 + i0;

                            for (int i2 = 0; i2 < size2; ++i2) {
                              const int ii2 = eimap2[offset2 + i2];
                              for (int i3 = 0; i3 < size3; ++i3, ++i0123) {
                                const int ii3 = eimap3[offset3 + i3];

                                const int i23 = ii2 * n3 + ii3;

                                if (abab) {
                                  const double integral_0123 = src0(i01,i23);
                                  data[i0123] = integral_0123;
                                } else if (abba) {
                                  const double integral_1023 = src1(i10,i23);
                                  data[i0123] = -integral_1023;
                                } else { // aaaa
                                  const double integral_0123 = src0(i01,i23);
                                  const double integral_1023 = src1(i10,i23);
                                  data[i0123] = integral_0123 - integral_1023;
                                }

                              }
                            }

                          } // antisymmetrize01

                          else { // do not antisymmetrize
                            MPQC_ASSERT(false); // should not happen, but only SMITH knows :-)
                          }

                        }
                      }

                      tensor_->put_block(tile_key, data);

#if 0
                const double sum = std::accumulate(data, data+size, 0.0);
                ExEnv::out0() << "tiles = (" << t0 << "," << t1 << "," << t2 << "," << t3
                              << ")  key = " << tile_key  << " sum = " << sum << std::endl;

                for(int i=0; i<size; ++i) {
                  ExEnv::out0() << "data[" << i << "] = " << data[i] << std::endl;
                }
#endif
                  } // if this task will process this tile

                } // t3
              } // t2
            } // t1
          } // t0

          info_->mem()->free_local_double(data);

        } // if this task can access the integrals

      } // MTensor<4>::convert() from RefSCMatrix or RefSymmSCMatrix

      /** copies contents of src to this. src is assumed to be a RefSCMatrix or a RefSymmSCMatrix
       *  indexing of elements in src may be different from that in this.
       *  arrays eimap0 through eimap3 provide mapping
       *
       *  src integrals are assumed to not be antisymmetrized
       */
      template <typename SrcType> void convert(const SrcType& src,
                                               const element_index_map& eimap0,
                                               const element_index_map& eimap1,
                                               element_ranges const* mapped_element_ranges = 0) {

        MPQC_ASSERT(NDIM == 2);

        // determine which tiles map to the allowed element ranges in src
        const std::vector<bool> tile0_in_erange = tiles_in_erange(range_[0], eimap0, (*mapped_element_ranges)[0]);
        const std::vector<bool> tile1_in_erange = tiles_in_erange(range_[1], eimap1, (*mapped_element_ranges)[1]);

        // split work over tasks assuming that all tasks can access src
        const int nproc_with_ints = info_->mem()->n();
        const int me = info_->mem()->me();

        const size_t ntiles0 = range_[0].second - range_[0].first;
        const size_t ntiles1 = range_[1].second - range_[1].first;

        // check that the needed index ranges are present in src
      #if 0  // assume that the user knows what she's doing
        element_range src_range[2];
        src_range[0] = this->map_range(range_[0], eimap0);
        src_range[1] = this->map_range(range_[1], eimap1);
        MPQC_ASSERT(src_range[0].second <= src->ni());
        MPQC_ASSERT(src_range[1].second <= src->nj());
      #endif

        // TODO check if equivalence of MTensor spaces is matched by their relationship in src

        // assume that src is accessible from all nodes
        if (1) {

          // how do I determine max size of the tiles?
          const size_t maxsize1 = info_->maxtilesize();
          const size_t maxtilesize = maxsize1 * maxsize1;
          double* data = info_->mem()->malloc_local_double(maxtilesize);

          for (long t0 = range_[0].first; t0 < range_[0].second; ++t0) {
            const long size0 = info_->get_range(t0);
            const long offset0 = info_->get_offset(t0);
            const long spin0 = info_->get_spin(t0);
            const bool in_erange0 = tile0_in_erange[t0];

            for (long t1 = range_[1].first; t1 < range_[1].second; ++t1) {
              const long size1 = info_->get_range(t1);
              const long offset1 = info_->get_offset(t1);
              const long spin1 = info_->get_spin(t1);
              const bool in_erange1 = tile1_in_erange[t1];

              // cartesian ordinal is used as a key for hashing tiles
              const long tile_key = (t1-range_[1].first) + ntiles1 * (t0-range_[0].first);

              // since all procs can access integrals, use locality
              if (tensor_->exists(tile_key) && tensor_->is_this_local(tile_key)) {

                if ((!in_erange0 || !in_erange1))
                  continue;

                long size = size0 * size1;
                std::fill(data, data+size, 0.0);

                long i01 = 0;
                for (int i0 = 0; i0 < size0; ++i0) {
                  const int ii0 = eimap0[offset0 + i0];
                  for (int i1 = 0; i1 < size1; ++i1, ++i01) {
                    const int ii1 = eimap1[offset1 + i1];

                    data[i01] = src(ii0, ii1);
                  }
                }

                tensor_->add_block(tile_key, data);

#if 0
                const double sum = std::accumulate(data, data+size, 0.0);
                ExEnv::out0() << "tiles = (" << t0 << "," << t1
                << ")  key = " << tile_key  << " sum = " << sum << endl;
#endif
              } // if this task will process this tile

            } // t1
          } // t0

          info_->mem()->free_local_double(data);

        } // if this task can access the integrals

      } // MTensor<2>::convert() from RefSCMatrix or RefSymmSCMatrix

      void print(const std::string& label, std::ostream& os = ExEnv::out0()) const {
        tensor_->print(label, os);
      }

    private:
      Info const* info_;
      Tensor* tensor_;
      tile_ranges range_;  // range of tiles that defines this tensor. See info_ for tiling info

      // given a range of indices, return pair<mmin,mmax> where mmin and mmax are the minimum and maximum mapped indices
      element_range map_range(const element_range& rng,
                              const element_index_map& eimap) {
        element_range result(eimap[rng.first],
                             eimap[rng.first]);
        for(element_index i=rng.first; i != rng.second; ++i) {
          const element_index ii = eimap[i];
          result.first = std::min(result.first, ii);
          result.second = std::max(result.second, ii);
        }
        return result;
      }

      // given tile range [start, fence), element map, and element_range, result[t] is true if all elements within
      // the tile map to inside element_range
      std::vector<bool>
      tiles_in_erange(const tile_range& range,
                      const element_index_map& eimap,
                      const element_range& erange) {

        std::vector<bool> result(range.second, false);
        for(long t=range.first; t!=range.second; ++t) {
          const long start = info_->get_offset(t);
          const long size = info_->get_range(t);
          const long fence = start + size;
          const element_range tile(start, fence);
          // map elements
          const element_range mapped_tile = map_range(tile, eimap);
          result[t] = in(mapped_tile,erange);
        }
        return result;
      }

  };

} // end of namespace sc

#include <chemistry/qc/ccr12/mtensor.timpl.h>

#endif // end of header guard

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:

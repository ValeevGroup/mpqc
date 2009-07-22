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

#ifdef __GNUG__
#pragma interface
#endif

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
        assert(range.size() == NDIM);
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
      void convert(const Ref<DistArray4>& src, DistArray4::tbint_type type,
                   const element_index_map& eimap0,
                   const element_index_map& eimap1,
                   const element_index_map& eimap2,
                   const element_index_map& eimap3,
                   element_ranges const* mapped_element_ranges = 0,
                   bool src_is_2301 = false);

      /** copies contents of src to this. src is assumed to be a RefSCMatrix or a RefSymmSCMatrix
       *  indexing of elements in src may be different from that in this.
       *  arrays eimap0 through eimap3 provide mapping
       *
       *  src integrals are assumed to not be antisymmetrized
       */
      template <typename SrcType> void convert(const SrcType& src,
                                               int n1,
                                               int n3,
                                               bool antisymmetric01,
                                               bool antisymmetric23,
                                               const element_index_map& eimap0,
                                               const element_index_map& eimap1,
                                               const element_index_map& eimap2,
                                               const element_index_map& eimap3,
                                               element_ranges const* mapped_element_ranges = 0) {

        assert(NDIM == 4);

        // split work over tasks assuming that all tasks can access src
        const int nproc_with_ints = info_->mem()->n();
        const int me = info_->mem()->me();

        const size_t ntiles0 = range_[0].second - range_[0].first;
        const size_t ntiles1 = range_[1].second - range_[1].first;
        const size_t ntiles2 = range_[2].second - range_[2].first;
        const size_t ntiles3 = range_[3].second - range_[3].first;

        // if space0 is equivalent to space1, or space2 equivalent to space3,
        // can compute antisymmetric integrals using only (space0 space1| space2 space3)
        // else also need (space0 space1| space3 space2), or an equivalent
        const bool space0_eq_space1 = (range_[0] == range_[1]);
        const bool space2_eq_space3 = (range_[2] == range_[3]);
        const bool can_always_antisymmetrize = (space0_eq_space1 || space2_eq_space3);
        assert(can_always_antisymmetrize == true);  // this is likely to break the logic downstream
                                                    // I clearly don't understand enough at the moment
        const tile_range space0_n_space1 = intersect(range_[0], range_[1]);
        const tile_range space2_n_space3 = intersect(range_[2], range_[3]);

        // check that the needed index ranges are present in src
      #if 0  // assume that the user knows what she's doing
        element_range src_range[4];
        src_range[0] = this->map_range(range_[0], eimap0);
        src_range[1] = this->map_range(range_[1], eimap1);
        src_range[2] = this->map_range(range_[2], eimap2);
        src_range[3] = this->map_range(range_[3], eimap3);
        assert(src_range[0].second <= src->ni());
        assert(src_range[1].second <= src->nj());
        assert(src_range[2].second <= src->nx());
        assert(src_range[3].second <= src->ny());
      #endif

        // TODO check if equivalence of MTensor spaces is matched by their relationship in src

        // assume that src is accessible from all nodes
        if (1) {

          // how do I determine max size of the tiles?
          const size_t maxtilesize = 50 * 50 * 50 * 50;
          double* data = info_->mem()->malloc_local_double(maxtilesize);

          for (long t0 = range_[0].first; t0 < range_[0].second; ++t0) {
            const long size0 = info_->get_range(t0);
            const int offset0 = eimap0[info_->get_offset(t0)];
            const int fence0 = offset0 + size0;
            const long spin0 = info_->get_spin(t0);
            const bool in_erange0 = mapped_element_ranges ? in(element_range(offset0,fence0),(*mapped_element_ranges)[0])
                                                          : true;

            for (long t1 = std::max(range_[1].first, t0); t1 < range_[1].second; ++t1) {
              const long size1 = info_->get_range(t1);
              const int offset1 = eimap1[info_->get_offset(t1)];
              const int fence1 = offset1 + size1;
              const long spin1 = info_->get_spin(t1);
              const bool in_erange1 = mapped_element_ranges ? in(element_range(offset1,fence1),(*mapped_element_ranges)[1])
                                                            : true;

              const bool aaaa = (spin0 == spin1);

              for (long t2 = range_[2].first; t2 < range_[2].second; ++t2) {
                const long size2 = info_->get_range(t2);
                const int offset2 = eimap2[info_->get_offset(t2)];
                const int fence2 = offset2 + size2;
                const long spin2 = info_->get_spin(t2);
                const bool in_erange2 = mapped_element_ranges ? in(element_range(offset2,fence2),(*mapped_element_ranges)[2])
                                                              : true;

                const bool abab = ((spin0 != spin1) && (spin0 == spin2));
                const bool abba = ((spin0 != spin1) && (spin0 != spin2));

                for (long t3 = std::max(range_[3].first, t2); t3 < range_[3].second; ++t3) {
                  const long size3 = info_->get_range(t3);
                  const int offset3 = eimap3[info_->get_offset(t3)];
                  const int fence3 = offset3 + size3;
                  const long spin3 = info_->get_spin(t3);
                  const bool in_erange3 = mapped_element_ranges ? in(element_range(offset3,fence3),(*mapped_element_ranges)[3])
                                                                : true;

                  // cartesian ordinal is used as a key for hashing tiles
                  const long tile_key = (t3-range_[3].first) +
                                        ntiles3 * ( (t2-range_[2].first) +
                                                    ntiles2 * ( (t1-range_[1].first) +
                                                                ntiles1 * (t0-range_[0].first)
                                                              )
                                                  );

                  // since all procs can access integrals, use locality
                  if (tensor_->exists(tile_key) && tensor_->is_this_local(tile_key)) {

                    if (mapped_element_ranges &&
                        (!in_erange0 || !in_erange1 || !in_erange2 || !in_erange3))
                        continue;

                    // antisymmetrization is only possible if:
                    // 1) t0 and t1 belong to the intersection of range[0] and range[1], or
                    // 2) t2 and t3 belong to the intersection of range[2] and range[3]
                    // if range[0] == range[1], clearly can antisymmetrize 0 and 1
                    // if range[2] == range[3], clearly can antisymmetrize 2 and 3
                    bool antisymmetrize01 = space0_eq_space1;
                    if (!antisymmetrize01)
                      antisymmetrize01 = in(t0,space0_n_space1) && in(t1,space0_n_space1);
                    bool antisymmetrize23 = space2_eq_space3;
                    if (!antisymmetrize23)
                      antisymmetrize23 = in(t2,space2_n_space3) && in(t3,space2_n_space3);
                    // prefer to antisymmetrize 23 since can do that with 1 block of DistArray4
                    if (antisymmetrize23) antisymmetrize01 = false;

                    long size = size0 * size1 * size2 * size3;
                    std::fill(data, data+size, 0.0);

                      long i0123 = 0;
                      for (int i0 = 0; i0 < size0; ++i0) {
                        const int ii0 = i0 + offset0;
                        for (int i1 = 0; i1 < size1; ++i1) {
                          const int ii1 = i1 + offset1;

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
                              const int ii2 = i2 + offset2;
                              for (int i3 = 0; i3 < size3; ++i3, ++i0123) {
                                const int ii3 = i3 + offset3;

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
                              const int ii2 = i2 + offset2;
                              for (int i3 = 0; i3 < size3; ++i3, ++i0123) {
                                const int ii3 = i3 + offset3;

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
                            assert(false); // should not happen, but only SMITH knows :-)
                          }

                        }
                      }

                      tensor_->put_block(tile_key, data);

      #if 0
                      const double sum = std::accumulate(data, data+size, 0.0);
                      ExEnv::out0() << "tiles = (" << t0 << "," << t1 << "," << t2 << "," << t3
                                    << ")  key = " << tile_key  << " sum = " << sum << endl;
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

        assert(NDIM == 2);

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
        assert(src_range[0].second <= src->ni());
        assert(src_range[1].second <= src->nj());
      #endif

        // TODO check if equivalence of MTensor spaces is matched by their relationship in src

        // assume that src is accessible from all nodes
        if (1) {

          // how do I determine max size of the tiles?
          const size_t maxtilesize = 50 * 50;
          double* data = info_->mem()->malloc_local_double(maxtilesize);

          for (long t0 = range_[0].first; t0 < range_[0].second; ++t0) {
            const long size0 = info_->get_range(t0);
            const int offset0 = eimap0[info_->get_offset(t0)];
            const int fence0 = offset0 + size0;
            const long spin0 = info_->get_spin(t0);
            const bool in_erange0 = mapped_element_ranges ? in(element_range(offset0,fence0),(*mapped_element_ranges)[0])
                                                          : true;

            for (long t1 = range_[1].first; t1 < range_[1].second; ++t1) {
              const long size1 = info_->get_range(t1);
              const int offset1 = eimap1[info_->get_offset(t1)];
              const int fence1 = offset1 + size1;
              const long spin1 = info_->get_spin(t1);
              const bool in_erange1 = mapped_element_ranges ? in(element_range(offset1,fence1),(*mapped_element_ranges)[1])
                                                            : true;

              // cartesian ordinal is used as a key for hashing tiles
              const long tile_key = (t1-range_[1].first) + ntiles1 * (t0-range_[0].first);

              // since all procs can access integrals, use locality
              if (tensor_->exists(tile_key) && tensor_->is_this_local(tile_key)) {

                if (mapped_element_ranges && (!in_erange0 || !in_erange1))
                  continue;

                long size = size0 * size1;
                std::fill(data, data+size, 0.0);

                long i01 = 0;
                for (int i0 = 0; i0 < size0; ++i0) {
                  const int ii0 = i0 + offset0;
                  for (int i1 = 0; i1 < size1; ++i1, ++i01) {
                    const int ii1 = i1 + offset1;

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
          result.first = std::min(rng.first, ii);
          result.second = std::max(rng.second, ii);
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

//
// mtensor.cc
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

#include <cassert>
#include <chemistry/qc/ccr12/mtensor.h>

using namespace std;
using namespace sc;

namespace sc {

template<>
void MTensor<4ul>::convert(const Ref<DistArray4>& src,
                         unsigned int type,
                         const std::vector<element_index>& eimap0,
                         const std::vector<element_index>& eimap1,
                         const std::vector<element_index>& eimap2,
                         const std::vector<element_index>& eimap3,
                         element_ranges const* mapped_element_ranges,
                         bool src_is_2301) {

  MPQC_ASSERT(mapped_element_ranges != 0);

  // split work over tasks which have access to integrals
  vector<int> proc_with_ints;
  const int nproc_with_ints = src->tasks_with_access(proc_with_ints);
  const bool all_procs_have_access = (nproc_with_ints == info_->mem()->n());
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

  if (src->has_access(me)) {

    // how do I determine max size of the tiles?
    const size_t maxsize1 = info_->maxtilesize();
    const size_t maxsize2 = maxsize1 * maxsize1;
    const size_t maxsize4 = maxsize2 * maxsize2;
    double* data = info_->mem()->malloc_local_double(maxsize4);
    double* buf1 = new double[maxsize2];
    double* buf2 = new double[maxsize2];

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

            // if all procs can access integrals, use locality
            if ((all_procs_have_access && tensor_->is_this_local(tile_key))
                || (!all_procs_have_access && tensor_->exists(tile_key))) {

                if (mapped_element_ranges &&
                    (!in_erange0 || !in_erange1 || !in_erange2 || !in_erange3))
                    continue;

                // if all processors have access, parallelized by locality
                // else do round-robin distribution of tasks
                if (!all_procs_have_access) {
                  const long task_proc = tile_key % nproc_with_ints;
                  if (task_proc != proc_with_ints[me])
                    continue;
                }

                long size = size0 * size1 * size2 * size3;
                std::fill(data, data+size, 0.0);

#if 0
                std::cout << "MTensor::convert: " << t0 << " " << t1 << " " << t2 << " " << t3 << " " << tile_key << std::endl;
#endif
                // if erange[0] == erange[1], clearly can antisymmetrize 0 and 1
                // if erange[2] == erange[3], clearly can antisymmetrize 2 and 3
                bool antisymmetrize01 = espace0_eq_espace1;
                bool antisymmetrize23 = espace2_eq_espace3;
                // prefer to antisymmetrize 23 since can do that with 1 block of DistArray4
                if (antisymmetrize23 && !src_is_2301) antisymmetrize01 = false;
                // if src is transposed: prefer to antisymmetrize 01 since can do that with 1 block of DistArray4
                if (antisymmetrize01 && src_is_2301) antisymmetrize23 = false;

                if (!src_is_2301) { // src is 0123

                  double* buf_0123 = buf1;
                  double* buf_0132 = buf2;
                  double* buf_1023 = buf2;

                long i0123 = 0;
                for (int i0 = 0; i0 < size0; ++i0) {
                  const int ii0 = i0 + offset0;
                  for (int i1 = 0; i1 < size1; ++i1) {
                    const int ii1 = i1 + offset1;

                    if (antisymmetrize23) {

                      if (abab || aaaa)
                      src->retrieve_pair_subblock(ii0, ii1, type,
                                                  offset2, fence2,
                                                  offset3, fence3,
                                                  buf_0123);
                      if (abba || aaaa)
                      src->retrieve_pair_subblock(ii0, ii1, type,
                                                  offset3, fence3,
                                                  offset2, fence2,
                                                  buf_0132);

                      int i23 = 0;
                      for (int i2 = 0; i2 < size2; ++i2) {
                        int i32 = i2;
                        for (int i3 = 0; i3 < size3; ++i3, ++i23, i32+=size2, ++i0123) {

                          if (abab) {
                            const double integral_0123 = buf_0123[i23];
                            data[i0123] = integral_0123;
                          } else if (abba) {
                            const double integral_0132 = buf_0132[i32];
                            data[i0123] = -integral_0132;
                          } else { // aaaa
                            const double integral_0123 = buf_0123[i23];
                            const double integral_0132 = buf_0132[i32];
                            data[i0123] = integral_0123 - integral_0132;
                          }

                        }
                      }

                    } // antisymmetrize23

                    else if (antisymmetrize01) { // means can't antisymmetrize 2 and 3

                      if (abab || aaaa)
                      src->retrieve_pair_subblock(ii0, ii1, type,
                                                  offset2, fence2,
                                                  offset3, fence3,
                                                  buf_0123);
                      if (abba || aaaa)
                      src->retrieve_pair_subblock(ii1, ii0, type,
                                                  offset2, fence2,
                                                  offset3, fence3,
                                                  buf_1023);

                      int i23 = 0;
                      for (int i2 = 0; i2 < size2; ++i2) {
                        for (int i3 = 0; i3 < size3; ++i3, ++i0123, ++i23) {

                          if (abab) {
                            const double integral_0123 = buf_0123[i23];
                            data[i0123] = integral_0123;
                          } else if (abba) {
                            const double integral_1023 = buf_1023[i23];
                            data[i0123] = -integral_1023;
                          } else { // aaaa
                            const double integral_0123 = buf_0123[i23];
                            const double integral_1023 = buf_1023[i23];
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
                } // src is 0123

                else { // src is 2301

                  double* buf_2301 = buf1;
                  double* buf_2310 = buf2;
                  double* buf_3201 = buf2;

                  const int size23 = size2*size3;
                  for (int i2 = 0; i2 < size2; ++i2) {
                    const int ii2 = i2 + offset2;
                    for (int i3 = 0; i3 < size3; ++i3) {
                      const int ii3 = i3 + offset3;

                      if (antisymmetrize01) {

                        if (abab || aaaa)
                          src->retrieve_pair_subblock(ii2, ii3, type,
                                                      offset0, fence0,
                                                      offset1, fence1,
                                                      buf_2301);
                        if (abba || aaaa)
                          src->retrieve_pair_subblock(ii2, ii3, type,
                                                      offset1, fence1,
                                                      offset0, fence0,
                                                      buf_2310);

                        int i01 = 0;
                        long i0123 = i2*size3 + i3;
                        for (int i0 = 0; i0 < size0; ++i0) {
                          int i10 = i0;
                          for (int i1 = 0; i1 < size1; ++i1, ++i01, i10+=size0, i0123+=size23) {

                            if (abab) {
                              const double integral_0123 = buf_2301[i01];
                              data[i0123] = integral_0123;
                            } else if (abba) {
                              const double integral_1023 = buf_2310[i10];
                              data[i0123] = -integral_1023;
                            } else { // aaaa
                              const double integral_0123 = buf_2301[i01];
                              const double integral_1023 = buf_2310[i10];
                              data[i0123] = integral_0123 - integral_1023;
                            }

                          }
                        }

                      } // antisymmetrize01

                      else if (antisymmetrize23) { // means can't antisymmetrize 0 and 1

                        if (abab || aaaa)
                          src->retrieve_pair_subblock(ii2, ii3, type,
                                                      offset0, fence0,
                                                      offset1, fence1,
                                                      buf_2301);
                        if (abba || aaaa)
                          src->retrieve_pair_subblock(ii3, ii2, type,
                                                      offset0, fence0,
                                                      offset1, fence1,
                                                      buf_3201);

                        int i01 = 0;
                        long i0123 = i2*size3 + i3;
                        for (int i0 = 0; i0 < size0; ++i0) {
                          for (int i1 = 0; i1 < size1; ++i1, ++i01, i0123+=size23) {

                            if (abab) {
                              const double integral_0123 = buf_2301[i01];
                              data[i0123] = integral_0123;
                            } else if (abba) {
                              const double integral_0132 = buf_3201[i01];
                              data[i0123] = -integral_0132;
                            } else { // aaaa
                              const double integral_0123 = buf_2301[i01];
                              const double integral_0132 = buf_3201[i01];
                              data[i0123] = integral_0123 - integral_0132;
                            }

                          }
                        }

                      } // antisymmetrize23

                      else { // do not antisymmetrize
                        MPQC_ASSERT(false); // should not happen, but only SMITH knows :-)
                      }

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

    delete[] buf1;
    delete[] buf2;
    info_->mem()->free_local_double(data);

  } // if this task can access the integrals

}

} // end of namespace sc

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:

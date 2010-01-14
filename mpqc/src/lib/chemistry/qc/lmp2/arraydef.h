
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

#ifndef _chemistry_qc_lmp2_arraydef_h
#define _chemistry_qc_lmp2_arraydef_h

#include <iomanip>

namespace sc {

namespace sma2 {
    template <int N>
    inline void Array<N>::print_local(std::ostream&o) const {
      for (typename Array<N>::blockmap_t::const_iterator
               aiter = blocks_.begin();
           aiter != blocks_.end();
           aiter++) {
          const BlockInfo<N> &info = aiter->first;
          const double *data = aiter->second;
          BlockIter<N> iter(indices(),aiter->first);
          bool first = true;
          for (iter.start(); iter.ready(); iter++) {
              o << "  ";
              for (int i=0; i<N; i++) {
                  o << " " << iter.index(i)
                      + index(i).block_offset(info.block(i));
                }
              o << " = ";
              if (allocated_) {
                  o << std::showpoint << std::fixed
                    << std::setw(11) << std::setprecision(8)
                    << data[iter.offset()];
                }
              else {
                  o << "unallocated";
                }
              o << " offset: ";
              for (int i=0; i<N; i++) {
                  o << " " << iter.index(i);
                }
              if (first) {
                  o << " block: ";
                  info.print(o);
                  first = false;
                }
              o << std::endl;
            }
        }
    }
    template <int N>
    inline void Array<N>::print(const sc::Ref<sc::MessageGrp> &grp,
                                bool distributed, std::ostream&o) const {
//       if (grp.null() || grp->me() == 0) {
//           for (int i = 0; i < N; i++) {
//               indices_[i].print(o);
//             }
//         }
      if (grp.nonnull() && distributed) {
          if (grp->me() == 0) o << "array " << name_ << ":" << std::endl;
          for (int i=0; i<grp->n(); i++) {
              if (grp->me() == i) {
                  o << "node " << i << ":" << std::endl;
                  print_local(o);
                  o << std::flush;
                }
              grp->sync();
              sleep(1);
            }
        }
      else if (grp.null() || (grp.nonnull() && grp->n() == 1)) {
          print_local(o);
        }
    }

    template <int N, class V>
    void
    apply_denominator(Array<N> &array, double denominator_offset,
                      const std::vector<V> &eigenvalues)
    {
      if (N != eigenvalues.size()) {
          throw std::runtime_error("apply_denominator: eigenvalues.size() != N");
        }
      bool inited = false;
      const bool print_factor_info = false;
      double smallest_factor = 0.0;
      double largest_factor = 0.0;
      const typename Array<N>::blockmap_t &bmap = array.blockmap();
      for (typename Array<N>::blockmap_t::const_iterator bmapiter = bmap.begin();
           bmapiter != bmap.end();
           bmapiter++) {
          const BlockInfo<N> &binfo = bmapiter->first;
          double *data = bmapiter->second;
          int off[N];
          for (int i=0; i<N; i++) {
              off[i] = array.index(i).block_offset(binfo.block(i));
            }
          BlockIter<N> biter(array.indices(),binfo);
          for (biter.start(); biter.ready(); biter++) {
              int offset = biter.offset();
              double denom = denominator_offset;
              for (int i=0; i<N; i++) {
                  denom += (*eigenvalues[i])[biter.index(i) + off[i]];
                }
              data[offset] /= denom;
              if (print_factor_info) {
                  if (!inited) {
                      largest_factor = smallest_factor = denom;
                      inited = true;
                    }
                  if (smallest_factor > denom) smallest_factor = denom;
                  if (largest_factor < denom) largest_factor = denom;
                }
            }
        }

      if (inited && print_factor_info) {
          sc::ExEnv::out0() << sc::indent
                            << "The denominators are in the range ["
                            << smallest_factor
                            << ", "
                            << largest_factor
                            << "]"
                            << std::endl;
        }
    }

}

}

#endif

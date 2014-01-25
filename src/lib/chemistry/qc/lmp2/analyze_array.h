
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

#ifndef _chemistry_qc_analyze_array_h
#define _chemistry_qc_analyze_array_h

#include <chemistry/qc/lmp2/sma.h>
#include <util/group/message.h>

namespace sc {

/** \brief Prints out information about an Array.
 * @param array the array to be analyzed.
 * @param name the name of the array to be printed in the output.
 * @param grp the MessageGrp.
 * @param distrubed whether or not the array's data is distributed.
 */
template <int N>
void
analyze_array(typename sma2::Array<N> &array, const char *name,
              const sc::Ref<sc::MessageGrp> &grp = 0, bool distributed  = false)
{
  if ((distributed && grp && grp->n() == 1)
      || grp == 0) {
      distributed = false;
    }
  double n_block_max_total = 1;
  double n_element_max_total = 1;
  for (int i=0; i<N; i++) {
      n_block_max_total *= array.index(i).nblock();
      n_element_max_total *= array.index(i).nindex();
    }
  double n_block_max = n_block_max_total;
  double n_element_max = n_element_max_total;
  if (distributed) {
      n_block_max /= grp->n();
      n_element_max /= grp->n();
    }
  double nblock = array.n_block();
  double nelement = array.n_element();
  double nblock_total = nblock;
  double nelement_total = nelement;
  if (distributed) {
      grp->sum(nblock_total);
      grp->sum(nelement_total);
      grp->max(nblock);
      grp->max(nelement);
    }
  sc::ExEnv::out0() << sc::indent
                    << sc::scprintf(
                        "  %10s %5.0f of %5.0f blocks %10.0f of %10.0f elements locally",
                        name, nblock, n_block_max, nelement, n_element_max)
                    << std::endl;
  if (distributed) {
      sc::ExEnv::out0() << sc::indent
                        << sc::scprintf(
                            "  %10s %5.0f of %5.0f blocks %10.0f of %10.0f elements globally",                            name, nblock_total, n_block_max_total, nelement_total,
                            n_element_max_total)
                        << std::endl;
    }
#if USE_STL_MULTIMAP
  sc::ExEnv::out0() << sc::indent
                << sc::scprintf(
                    "    uses %11.3f MB for data locally",
                    nelement*sizeof(double)*1e-6)
                << std::endl;
#else
  sc::ExEnv::out0() << sc::indent
                << sc::scprintf(
                    "    uses %11.3f MB for data and %11.3f MB for the map locally",
                    nelement*sizeof(double)*1e-6,
                    nblock*sizeof(AVLMMapNode<sma2::BlockInfo<N>,double*>)*1e-6)
                << std::endl;
#endif
}

}

#endif


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

#include <math.h>

#include <chemistry/qc/lmp2/sma.h>

namespace sc {
namespace sma2 {

static double *
get_element_pointer(Array<4> & ao_integrals,
                      int p, int q, int r, int s)
{
  int indices[4];
  indices[0] = p;
  indices[1] = q;
  indices[2] = r;
  indices[3] = s;

  BlockInfo<4> binfo;

  for (int i=0; i<4; i++) {
      binfo.block(i) = ao_integrals.index(i).index_to_block(indices[i]);
    }

  Array<4>::blockmap_t::const_iterator block
      = ao_integrals.blockmap().find(binfo);
  if (block == ao_integrals.blockmap().end()) return 0;
  double *data = block->second;

  BlockIter<4> biter(ao_integrals.indices(), binfo);

  for (int i=0; i<4; i++) {
      biter.index(i) = ao_integrals.index(i).index_to_offset(indices[i]);
    }

  return &data[biter.offset()];
}

void
pack_2e_integrals_into_array(Array<4> & ao_integrals,
                             const sc::Ref<sc::TwoBodyInt> &tbint)
{
  const double *buffer = tbint->buffer();
  sc::Ref<sc::GaussianBasisSet> basis = tbint->basis();
  double log2bound = log(DBL_EPSILON)/log(2.0);

  for (int i=0; i<4; i++) {
      if (ao_integrals.index(i).nindex() != basis->nbasis()) {
          throw std::runtime_error("bad range in 2e array");
        }
    }

  ao_integrals.allocate_blocks();

  int nshell = basis->nshell();
  for(int P = 0; P < nshell; P++) {
    int nump = basis->shell(P).nfunction();
    for(int Q = 0; Q < nshell; Q++) {
      int numq = basis->shell(Q).nfunction();
      if (tbint->log2_shell_bound(P,Q) < log2bound) continue;
      for(int R = 0; R < nshell; R++) {
        int numr = basis->shell(R).nfunction();
        if (tbint->log2_shell_bound(P,Q,R) < log2bound) continue;
        for(int S = 0; S < nshell; S++) {
          int nums = basis->shell(S).nfunction();
          if (tbint->log2_shell_bound(P,Q,R,S) < log2bound) continue;

          bool computed = false;

          int index = 0;
          for(int p=0; p < nump; p++) {
            int op = basis->shell_to_function(P)+p;
            int mop = ao_integrals
                     .index(0).basis_index_to_range_index(op);
            for(int q = 0; q < numq; q++) {
              int oq = basis->shell_to_function(Q)+q;
              int moq = ao_integrals
                       .index(1).basis_index_to_range_index(oq);
              for(int r = 0; r < numr; r++) {
                int oor = basis->shell_to_function(R)+r;
                int mor = ao_integrals
                         .index(2).basis_index_to_range_index(oor);
                for(int s = 0; s < nums; s++,index++) {
                  int os = basis->shell_to_function(S)+s;
                  int mos = ao_integrals
                           .index(3).basis_index_to_range_index(os);

                  double *dat
                    = get_element_pointer(ao_integrals,mop,moq,mor,mos);
                  if (dat != 0) {
                      if (!computed) {
                          tbint->compute_shell(P,Q,R,S);
                          computed = true;
                        }
                      *dat = buffer[index];
                    }
                }
              }
            }
          }

        }
      }
    }
  }

  ao_integrals.compute_bounds();
}

void
pack_2e_integrals_into_shell_blocked_array(Array<4> & ao_integrals,
                             const sc::Ref<sc::TwoBodyInt> &tbint)
{
  const double *buffer = tbint->buffer();
  sc::Ref<sc::GaussianBasisSet> basis = tbint->basis();
  double log2bound = log(DBL_EPSILON)/log(2.0);

  for (int i=0; i<4; i++) {
      if (ao_integrals.index(i).nindex() != basis->nbasis()) {
          throw std::runtime_error("bad range in 2e array (nindex)");
        }
      if (ao_integrals.index(i).nblock() != basis->nshell()) {
          throw std::runtime_error("bad range in 2e array (nshell)");
        }
      for (int j=0; j<basis->nshell(); j++) {
          if (ao_integrals.index(i).block_size(j) != basis->shell(j).nfunction()) {
              throw std::runtime_error("bad range in 2e array (nfunction_in_shell)");
            }
        }
    }

  ao_integrals.allocate_blocks();

  const Array<4>::blockmap_t &ao_blockmap = ao_integrals.blockmap();

  for (Array<4>::blockmap_t::const_iterator ia = ao_blockmap.begin();
       ia != ao_blockmap.end();
       ia++) {
      const BlockInfo<4> &iablockinfo = ia->first;
      int P = iablockinfo.block(0);
      int Q = iablockinfo.block(1);
      int R = iablockinfo.block(2);
      int S = iablockinfo.block(3);
      int n = ao_integrals.block_size(iablockinfo);
      if (tbint->log2_shell_bound(P,Q,R,S) >= log2bound) {
          tbint->compute_shell(P,Q,R,S);
          memcpy(ia->second, buffer, n*sizeof(*buffer));
        }
      else {
          memset(ia->second, 0, n*sizeof(*buffer));
        }
    }

  ao_integrals.compute_bounds();
}

}
}

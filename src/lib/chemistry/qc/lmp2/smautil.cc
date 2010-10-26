
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

#include <stdexcept>

#include <chemistry/qc/lmp2/sma.h>

namespace sc {

namespace sma2 {

void
extract_diagonal(Array<2> &array, std::vector<double> &diag)
{
  if (array.index(0) != array.index(1)) {
      throw std::invalid_argument("sma::extract_diagonal: bad array");
    }
  const Range &r = array.index(0);
  diag.resize(r.nindex());
  BlockInfo<2> bi;
  int nblock = r.nblock();
  int ij=0;
  for (int i=0; i<nblock; i++) {
      int blocksize = r.block_size(i);
      bi.block(0) = i;
      bi.block(1) = i;
      Array<2>::blockmap_t::const_iterator blockii
          = array.blockmap().find(bi);
      if (blockii == array.blockmap().end()) continue;
      double *data = blockii->second;
      int offset = 0;
      for (int j=0; j<blocksize; j++,ij++,offset+=blocksize+1) {
          diag[ij] = data[offset];
        }
    }
}

void
scale_diagonal(Array<2> &array, double factor)
{
  if (array.index(0) != array.index(1)) {
      throw std::invalid_argument("sma::scale_diagonal: bad array");
    }
  const Range &r = array.index(0);
  BlockInfo<2> bi;
  int nblock = r.nblock();
  int ij=0;
  for (int i=0; i<nblock; i++) {
      int blocksize = r.block_size(i);
      bi.block(0) = i;
      bi.block(1) = i;
      Array<2>::blockmap_t::const_iterator blockii
          = array.blockmap().find(bi);
      if (blockii == array.blockmap().end()) continue;
      double *data = blockii->second;
      int offset = 0;
      for (int j=0; j<blocksize; j++,ij++,offset+=blocksize+1) {
          data[offset] *= factor;
        }
    }
}

}

}

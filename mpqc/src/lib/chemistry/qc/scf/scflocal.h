//
// scflocal.h --- local inline functions
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Edward Seidl <seidl@janed.com>
// Maintainer: LPS
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

#ifndef _chemistry_qc_scf_scflocal_h
#define _chemistry_qc_scf_scflocal_h

#include <math/scmat/block.h>

namespace sc {

static inline double *
get_tri_block(SCMatrixBlock* blk,
              int& istart, int& iend, int& jstart, int& jend, int& sub)
{
  double *data=0;
  
  if (dynamic_cast<SCMatrixLTriBlock*>(blk)) {
    SCMatrixLTriBlock *lblk = dynamic_cast<SCMatrixLTriBlock*>(blk);
    istart = lblk->start; iend=lblk->end;
    jstart = lblk->start; jend=lblk->end;
    data = lblk->data;
    sub=0;
  } else if (dynamic_cast<SCMatrixLTriSubBlock*>(blk)) {
    SCMatrixLTriSubBlock *lblk = dynamic_cast<SCMatrixLTriSubBlock*>(blk);
    istart = lblk->istart; iend=lblk->iend;
    jstart = lblk->jstart; jend=lblk->jend;
    data = lblk->data;
    sub=1;
  }

  return data;
}

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:


#ifndef _chemistry_qc_scf_scflocal_h
#define _chemistry_qc_scf_scflocal_h

#include <math/scmat/block.h>

static inline double *
get_tri_block(SCMatrixBlock* blk,
              int& istart, int& iend, int& jstart, int& jend, int& sub)
{
  double *data=0;
  
  if (SCMatrixLTriBlock::castdown(blk)) {
    SCMatrixLTriBlock *lblk = SCMatrixLTriBlock::castdown(blk);
    istart = lblk->start; iend=lblk->end;
    jstart = lblk->start; jend=lblk->end;
    data = lblk->data;
    sub=0;
  } else if (SCMatrixLTriSubBlock::castdown(blk)) {
    SCMatrixLTriSubBlock *lblk = SCMatrixLTriSubBlock::castdown(blk);
    istart = lblk->istart; iend=lblk->iend;
    jstart = lblk->jstart; jend=lblk->jend;
    data = lblk->data;
    sub=1;
  }

  return data;
}

#endif

// Local Variables:
// mode: c++
// eval: (c-set-style "ETS")


#ifndef _chemistry_qc_scf_scflocal_h
#define _chemistry_qc_scf_scflocal_h

#include <math/scmat/block.h>

static inline double *
get_tri_block(SCMatrixBlock* blk,
              int& istart, int& iend, int& jstart, int& jend, int&tri)
{
  double *data=0;
  
  if (SCMatrixLTriBlock::castdown(blk)) {
    SCMatrixLTriBlock *lblk = SCMatrixLTriBlock::castdown(blk);
    istart = lblk->start; iend=lblk->end;
    jstart = lblk->start; jend=lblk->end;
    data = lblk->data;
    tri=1;
  } else if (SCMatrixLTriSubBlock::castdown(blk)) {
    SCMatrixLTriSubBlock *lblk = SCMatrixLTriSubBlock::castdown(blk);
    istart = lblk->istart; iend=lblk->iend;
    jstart = lblk->jstart; jend=lblk->jend;
    data = lblk->data;
    tri=1;
  } else if (SCMatrixRectSubBlock::castdown(blk)) {
    SCMatrixRectSubBlock *lblk = SCMatrixRectSubBlock::castdown(blk);
    istart = lblk->istart; iend=lblk->iend;
    jstart = lblk->jstart; jend=lblk->jend;
    data = lblk->data;
    tri=0;
  }

  return data;
}

#endif

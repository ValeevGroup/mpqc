
#ifndef _math_scmat_util_h
#define _math_scmat_util_h

#include <math/scmat/elemop.h>
#include <math/scmat/block.h>

namespace sc {

  void
  scmat_perform_op_on_blocks(const Ref<SCElementOp>& op,
                             const Ref<SCMatrixBlockList> &blocklist);

};

#endif

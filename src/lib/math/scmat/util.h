
#ifndef _math_scmat_util_h
#define _math_scmat_util_h

#include <math/scmat/elemop.h>
#include <math/scmat/block.h>

namespace sc {

  void
  scmat_perform_op_on_blocks(const Ref<SCElementOp>& op,
                             const Ref<SCMatrixBlockList> &blocklist);

  /// Canonicalize phases of SCMatrix A
  /// phases are canonical when the largest-magnitude coefficient in each column is positive
  void canonicalize_column_phases(RefSCMatrix& A);

};

#endif

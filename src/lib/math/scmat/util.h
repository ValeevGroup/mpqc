
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

  /// Compute and print out neatly various matrix norms of A
  template <class RefSCMat>
    void print_scmat_norms(const RefSCMat& A, const std::string& label, std::ostream& os = ExEnv::out0())
    {
      Ref<SCElementMaxAbs> maxabs_op(new SCElementMaxAbs);
      A.element_op(maxabs_op);
      const double maxabs = maxabs_op->result();

      Ref<SCElementKNorm> onenorm_op(new SCElementKNorm(1));
      A.element_op(onenorm_op);
      const double onenorm = onenorm_op->result();

      Ref<SCElementKNorm> twonorm_op(new SCElementKNorm(2));
      A.element_op(twonorm_op);
      const double twonorm = twonorm_op->result();

      os << indent << "Norms of " << label << std::endl;
      os << indent << "------------------------" << std::endl;
      os << indent << "||A||_{\\infty} = " << scprintf("%10.5lf",maxabs) << std::endl;
      os << indent << "||A||_1        = " << scprintf("%10.5lf",onenorm) << std::endl;
      os << indent << "||A||_2        = " << scprintf("%10.5lf",twonorm) << std::endl << std::endl;
    }

};

#endif

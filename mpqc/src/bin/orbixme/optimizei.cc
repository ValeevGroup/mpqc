
#include <optimizei.h>

#include <math/optimize/opt.h>

// force linkage of QNewtonOpt
#include <math/optimize/qnewton.h>
const ClassDesc &fl0 = QNewtonOpt::class_desc_;

C_OptimizeImpl::C_OptimizeImpl()
{
}

C_OptimizeImpl::~C_OptimizeImpl()
{
}

Optimize *
C_OptimizeImpl::opt()
{
  Optimize *ret;
  ret = Optimize::castdown(dc_);
  return ret;
}

long
C_OptimizeImpl::optimize(CORBA_Environment &)
{
  if (!opt()) return 0;
  return opt()->optimize();
}


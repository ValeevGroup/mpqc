
#include <optimizei.h>

#include <math/optimize/opt.h>

// force linkage of QNewtonOpt
#include <math/optimize/qnewton.h>
static ForceLink<QNewtonOpt> fl0;

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
  ret = dynamic_cast<Optimize*>(dc_);
  return ret;
}

long
C_OptimizeImpl::optimize(CORBA_Environment &)
{
  if (!opt()) return 0;
  return opt()->optimize();
}


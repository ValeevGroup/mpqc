
#include <functioni.h>

#include <math/optimize/nlp.h>

C_FunctionImpl::C_FunctionImpl()
{
}

C_FunctionImpl::~C_FunctionImpl()
{
}

NLP2 *
C_FunctionImpl::func()
{
  NLP2 *ret;
  ret = dynamic_cast<NLP2*>(dc_);
  return ret;
}

double
C_FunctionImpl::value(CORBA_Environment &)
{
  if (!func()) return 0.0;
  return func()->value();
}



#ifdef __GNUC__
#pragma implementation
#endif

#include <chemistry/qc/wfn/hcore.h>
#include <chemistry/qc/intv2/integralv2.h>

#define CLASSNAME AccumHCore
#define PARENTS public AccumDIH
#define HAVE_CTOR
#define HAVE_STATEIN_CTOR
#define HAVE_KEYVAL_CTOR
#include <util/class/classi.h>

void *
AccumHCore::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = AccumDIH::_castdown(cd);
  return do_castdowns(casts,cd);
}

AccumHCore::AccumHCore()
{
}

AccumHCore::AccumHCore(StateIn&s) :
  AccumDIH(s)
{
}

AccumHCore::AccumHCore(const RefKeyVal& keyval) :
  AccumDIH(keyval)
{
}

AccumHCore::~AccumHCore()
{
}

void
AccumHCore::save_data_state(StateOut& s)
{
  AccumDIH::save_data_state(s);
}

void
AccumHCore::accum(const RefSymmSCMatrix& h)
{
  h.assign(0.0);
  RefSCElementOp op = new GaussianKineticIntv2(_basis_set);
  h.element_op(op);
  op=0;
  
  op = new GaussianNuclearIntv2(_basis_set);
  h.element_op(op);
  op=0;
}

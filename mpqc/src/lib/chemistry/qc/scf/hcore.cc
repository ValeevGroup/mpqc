
#ifdef __GNUC__
#pragma implementation
#endif

#include <math/optimize/diis.h>
#include <chemistry/qc/scf/hcore.h>
#include <chemistry/qc/integral/integralv2.h>

#define CLASSNAME AccumHCore
#define PARENTS public AccumDIH
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

AccumHCore::AccumHCore(StateIn&s) :
  AccumDIH(s)
{
}

AccumHCore::AccumHCore(const RefKeyVal& keyval) :
  AccumDIH(keyval)
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
  RefSCElementOp op = new GaussianKineticIntv2(_basis_set, _molecule);
  h.element_op(op);
  op = new GaussianNuclearIntv2(_basis_set, _molecule);
  h.element_op(op);
}

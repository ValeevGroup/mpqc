
#ifdef __GNUC__
#pragma implementation
#endif

extern "C" {
#  include <math.h>
}

#include "update.h"
#include <util/keyval/keyval.h>

SavableState_REF_def(HessianUpdate);
#define CLASSNAME HessianUpdate
#define PARENTS virtual_base public SavableState
#include <util/state/statei.h>
#include <util/class/classia.h>

/////////////////////////////////////////////////////////////////////////
// HessianUpdate

void *
HessianUpdate::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SavableState::_castdown(cd);
  return do_castdowns(casts,cd);
}

HessianUpdate::HessianUpdate() : inverse_hessian_(0)
{
}

HessianUpdate::HessianUpdate(StateIn&s):
  SavableState(s)
{
  s.get(inverse_hessian_);
}

HessianUpdate::HessianUpdate(const RefKeyVal&keyval) :
  inverse_hessian_(0)
{
}

HessianUpdate::~HessianUpdate()
{
}

void
HessianUpdate::save_data_state(StateOut&s)
{
  s.put(inverse_hessian_);
}

void
HessianUpdate::set_inverse(void)
{
  inverse_hessian_ = 1;
}

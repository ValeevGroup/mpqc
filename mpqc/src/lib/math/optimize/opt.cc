
#ifdef __GNUC__
#pragma implementation
#endif

extern "C" {
#  include <math.h>
}

#include "opt.h"
#include <util/keyval/keyval.h>

SavableState_REF_def(Optimize);
#define CLASSNAME Optimize
#define PARENTS virtual public SavableState
#include <util/state/statei.h>
#include <util/class/classia.h>

SavableState_REF_def(LineOpt);
#define CLASSNAME LineOpt
#define PARENTS public Optimize
#include <util/state/statei.h>
#include <util/class/classia.h>

/////////////////////////////////////////////////////////////////////////
// Optimize

void *
Optimize::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SavableState::_castdown(cd);
  return do_castdowns(casts,cd);
}

Optimize::Optimize()
{
}

Optimize::Optimize(StateIn&s):
  SavableState(s,Optimize::class_desc_)
{
  s.get(max_iterations_);
  n_iterations_ = 0;
}

Optimize::Optimize(const RefKeyVal&keyval)
{
  max_iterations_ = keyval->intvalue("max_iterations");
  if (keyval->error() != KeyVal::OK) max_iterations_ = 10;
  n_iterations_ = 0;
}

Optimize::~Optimize()
{
}

void
Optimize::save_data_state(StateOut&s)
{
  s.put(max_iterations_);
}

void
Optimize::init()
{
  n_iterations_ = 0;
}

int
Optimize::optimize()
{
  int result;
  while((!(result = update())) && (n_iterations_ < max_iterations_))
      n_iterations_++;
  return result;
}

/////////////////////////////////////////////////////////////////////////
// LineOpt

void *
LineOpt::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = Optimize::_castdown(cd);
  return do_castdowns(casts,cd);
}

LineOpt::LineOpt()
{
}

LineOpt::LineOpt(StateIn&s):
  SavableState(s,LineOpt::class_desc_),
  Optimize(s)
{
  search_direction_.restore_state(s);
}

LineOpt::LineOpt(const RefKeyVal&keyval):
  Optimize(keyval)
{
}

LineOpt::~LineOpt()
{
}

void
LineOpt::save_data_state(StateOut&s)
{
  Optimize::save_data_state(s);
  search_direction_.save_state(s);
}

void
LineOpt::set_search_direction(RefSCVector&s)
{
  search_direction_ = s.copy();
}


#ifdef __GNUC__
#pragma implementation
#endif

extern "C" {
#  include <math.h>
}

#include <math/optimize/opt.h>
#include <util/keyval/keyval.h>

SavableState_REF_def(Optimize);
#define CLASSNAME Optimize
#define PARENTS virtual_base public SavableState
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

Optimize::Optimize() :
  ckpt_(0), ckpt_file(0)
{
}

Optimize::Optimize(StateIn&s):
  SavableState(s)
{
  s.get(ckpt_);
  s.getstring(ckpt_file);
  s.get(max_iterations_);
  n_iterations_ = 0;
}

Optimize::Optimize(const RefKeyVal&keyval)
{
  ckpt_ = keyval->booleanvalue("checkpoint");
  if (keyval->error() != KeyVal::OK) ckpt_ = 0;
  ckpt_file = keyval->pcharvalue("checkpoint_file");
  if (keyval->error() != KeyVal::OK) {
    ckpt_file = new char[13];
    strcat(ckpt_file,"opt_ckpt.dat");
  }
  max_iterations_ = keyval->intvalue("max_iterations");
  if (keyval->error() != KeyVal::OK) max_iterations_ = 10;
  n_iterations_ = 0;
}

Optimize::~Optimize()
{
  if (ckpt_file) delete[] ckpt_file;
  ckpt_file=0;
}

void
Optimize::save_data_state(StateOut&s)
{
  s.put(ckpt_);
  s.putstring(ckpt_file);
  s.put(max_iterations_);
}

void
Optimize::init()
{
  n_iterations_ = 0;
}

void
Optimize::set_checkpoint()
{
  ckpt_=1;
}

void
Optimize::set_max_iterations(int mi)
{
  max_iterations_ = mi;
}

void
Optimize::set_checkpoint_file(const char *path)
{
  if (ckpt_file) delete[] ckpt_file;
  if (path) {
    ckpt_file = new char[strlen(path)+1];
    strcpy(ckpt_file,path);
  } else
    ckpt_file=0;
}
  

#ifndef OPTSTATEOUT
#define OPTSTATEOUT StateOutBinXDR
#endif

int
Optimize::optimize()
{
  int result;
  while((n_iterations_ < max_iterations_) && (!(result = update()))) {
      n_iterations_++;
      if (ckpt_) {
        OPTSTATEOUT so(ckpt_file);
        this->save_state(so);
      }
    }
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
  Optimize(s)
  maybe_SavableState(s)
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

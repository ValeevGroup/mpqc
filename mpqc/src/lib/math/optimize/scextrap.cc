
#ifdef __GNUC__
#pragma implementation
#endif

#include <math/optimize/scextrap.h>

SavableState_REF_def(SCExtrapData);

#define CLASSNAME SCExtrapData
#define PARENTS public SavableState
#include <util/class/classia.h>

void *
SCExtrapData::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SavableState::_castdown(cd);
  return do_castdowns(casts,cd);
}

SCExtrapData::SCExtrapData()
{
}

SCExtrapData::SCExtrapData(StateIn& s) :
  SavableState(s)
{
}

SCExtrapData::~SCExtrapData()
{
}

void
SCExtrapData::save_data_state(StateOut& s)
{
}

////////////////////////////////////////////////////////////////////////////

SavableState_REF_def(SCExtrapError);

#define CLASSNAME SCExtrapError
#define PARENTS public SavableState
#include <util/class/classia.h>
void *
SCExtrapError::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SavableState::_castdown(cd);
  return do_castdowns(casts,cd);
}

SCExtrapError::SCExtrapError()
{
}

SCExtrapError::SCExtrapError(StateIn& s) :
  SavableState(s)
{
}

SCExtrapError::~SCExtrapError()
{
}

void
SCExtrapError::save_data_state(StateOut& s)
{
}

////////////////////////////////////////////////////////////////////////////

SavableState_REF_def(SelfConsistentExtrapolation);

#define CLASSNAME SelfConsistentExtrapolation
#define PARENTS public SavableState
#include <util/class/classia.h>
void *
SelfConsistentExtrapolation::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SavableState::_castdown(cd);
  return do_castdowns(casts,cd);
}

SelfConsistentExtrapolation::SelfConsistentExtrapolation()
{
  errorset_ = 0;
  tolerance_ = 1.0e-8;
}

SelfConsistentExtrapolation::SelfConsistentExtrapolation(StateIn& s) :
  SavableState(s)
{
  s.get(error_);
  s.get(errorset_);
  s.get(tolerance_);
}

SelfConsistentExtrapolation::SelfConsistentExtrapolation(
    const RefKeyVal&keyval)
{
  errorset_ = 0;
  tolerance_ = keyval->doublevalue("tolerance");
  if (keyval->error() != KeyVal::OK) tolerance_ = 1.0e-8;
}

SelfConsistentExtrapolation::~SelfConsistentExtrapolation()
{
}

void
SelfConsistentExtrapolation::save_data_state(StateOut& s)
{
  s.put(error_);
  s.put(errorset_);
  s.put(tolerance_);
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// eval: (c-set-style "CLJ")

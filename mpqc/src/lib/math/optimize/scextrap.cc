
#ifdef __GNUC__
#pragma implementation
#endif

#include <math/optimize/scextrap.h>

#define CLASSNAME SCExtrapData
#define PARENTS public DescribedClass
#include <util/class/classia.h>
void *
SCExtrapData::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = DescribedClass::_castdown(cd);
  return do_castdowns(casts,cd);
}

SCExtrapData::~SCExtrapData()
{
}

#define CLASSNAME SCExtrapError
#define PARENTS public DescribedClass
#include <util/class/classia.h>
void *
SCExtrapError::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = DescribedClass::_castdown(cd);
  return do_castdowns(casts,cd);
}

SCExtrapError::~SCExtrapError()
{
}

#define CLASSNAME SelfConsistentExtrapolation
#define PARENTS public DescribedClass
#include <util/class/classia.h>
void *
SelfConsistentExtrapolation::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = DescribedClass::_castdown(cd);
  return do_castdowns(casts,cd);
}

SelfConsistentExtrapolation::SelfConsistentExtrapolation()
{
  errorset_ = 0;
  tolerance_ = 1.0e-8;
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

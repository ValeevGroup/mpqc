
#ifdef __GNUC__
#pragma implementation
#endif

#include <chemistry/qc/basis/symmint.h>

////////////////////////////////////////////////////////////////////////////
// SymmOneBodyIntIter

SymmOneBodyIntIter::SymmOneBodyIntIter(const RefOneBodyInt& ints,
                                       const RefPetiteList& p) :
  OneBodyIntIter(ints), pl(p)
{
}

SymmOneBodyIntIter::~SymmOneBodyIntIter()
{
}

void
SymmOneBodyIntIter::next()
{
  OneBodyIntIter::next();
  while (!pl->lambda(icur,jcur) && OneBodyIntIter::ready())
    OneBodyIntIter::next();
}

void
SymmOneBodyIntIter::start(int ist, int jst, int ien, int jen)
{
  OneBodyIntIter::start(ist,jst,ien,jen);
  while (!pl->lambda(icur,jcur) && OneBodyIntIter::ready())
    OneBodyIntIter::next();
}

double
SymmOneBodyIntIter::scale() const
{
  return (double) pl->lambda(icur,jcur) / (double) pl->order();
}

////////////////////////////////////////////////////////////////////////////
// SymmTwoBodyIntIter

SymmTwoBodyIntIter::SymmTwoBodyIntIter(const RefTwoBodyInt& ints,
                                       const RefPetiteList& p) :
  TwoBodyIntIter(ints), pl(p)
{
}

SymmTwoBodyIntIter::~SymmTwoBodyIntIter()
{
}

// very inefficient...fix later
void
SymmTwoBodyIntIter::next()
{
  TwoBodyIntIter::next();
  while (TwoBodyIntIter::ready() &&
         !pl->in_p4(i_offset(icur)+jcur, i_offset(kcur)+lcur,
                    icur,jcur,kcur,lcur))
    TwoBodyIntIter::next();
}

// very inefficient...fix later
void
SymmTwoBodyIntIter::start()
{
  TwoBodyIntIter::start();
  while (TwoBodyIntIter::ready() &&
         !pl->in_p4(i_offset(icur)+jcur, i_offset(kcur)+lcur,
                    icur,jcur,kcur,lcur))
    TwoBodyIntIter::next();
}

// very inefficient...fix later
double
SymmTwoBodyIntIter::scale() const
{
  return (double)
    pl->in_p4(i_offset(icur)+jcur, i_offset(kcur)+lcur,icur,jcur,kcur,lcur);
}

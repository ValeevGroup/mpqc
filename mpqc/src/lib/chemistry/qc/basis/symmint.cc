
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
  while (!pl->lambda(icur,jcur))
    OneBodyIntIter::next();
}

void
SymmOneBodyIntIter::start(int ist, int jst, int ien, int jen)
{
  OneBodyIntIter::start(ist,jst,ien,jen);
  while (!pl->lambda(icur,jcur))
    OneBodyIntIter::next();
}

double
SymmOneBodyIntIter::scale() const
{
  return (double) pl->lambda(icur,jcur) / (double) pl->order();
}

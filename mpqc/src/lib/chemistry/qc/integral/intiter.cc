
#ifdef __GNUC__
#pragma implementation
#endif

#include <chemistry/qc/integral/intiter.h>

OneBodyIntIter::OneBodyIntIter()
{
  reset(0,0,0,0);
}

OneBodyIntIter::OneBodyIntIter(int is, int ie, int js, int je)
{
  reset(is,ie,js,je);
}

OneBodyIntIter::~OneBodyIntIter()
{
}

void
OneBodyIntIter::reset(int is, int ie, int js, int je)
{
  istart=is;
  jstart=js;
  iend=ie;
  jend=je;
}

void
OneBodyIntIter::start()
{
  icur=istart;
  jcur=jstart;
}

int
OneBodyIntIter::ready()
{
  return (icur < iend);
}

void
OneBodyIntIter::next()
{
  if (jcur < jend-1) {
    jcur++;
  } else {
    jcur=jstart;
    icur++;
  }
}

void
OneBodyIntIter::start_ltri()
{
  icur=istart;
  jcur=istart;
}

int
OneBodyIntIter::ready_ltri()
{
  return (icur < iend);
}

void
OneBodyIntIter::next_ltri()
{
  if (jcur < icur) {
    jcur++;
  } else {
    jcur=istart;
    icur++;
  }
}

int
OneBodyIntIter::ishell() const
{
  return icur;
}

int
OneBodyIntIter::jshell() const
{
  return jcur;
}

double
OneBodyIntIter::scale() const
{
  return 1.0;
}

////////////////////////////////////////////////////////////////////////////
// SymmOneBodyIntIter

SymmOneBodyIntIter::SymmOneBodyIntIter(const RefSymmGaussianBasisSet& gbs) :
  pl(gbs->petite_list())
{
}

SymmOneBodyIntIter::~SymmOneBodyIntIter()
{
}

void
SymmOneBodyIntIter::next()
{
  OneBodyIntIter::next();
  while (!pl.lambda(icur,jcur))
    OneBodyIntIter::next();
}

void
SymmOneBodyIntIter::next_ltri()
{
  OneBodyIntIter::next_ltri();
  while (!pl.lambda(icur,jcur))
    OneBodyIntIter::next_ltri();
}

void
SymmOneBodyIntIter::start()
{
  OneBodyIntIter::start();
  while (!pl.lambda(icur,jcur))
    OneBodyIntIter::next();
}

void
SymmOneBodyIntIter::start_ltri()
{
  OneBodyIntIter::start_ltri();
  while (!pl.lambda(icur,jcur))
    OneBodyIntIter::next_ltri();
}

double
SymmOneBodyIntIter::scale() const
{
  return (double) pl.lambda(icur,jcur) / (double) pl.order();
}

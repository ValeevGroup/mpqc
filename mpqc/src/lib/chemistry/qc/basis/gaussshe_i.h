
#ifdef INLINE_FUNCTIONS
#define INLINE inline
#else
#define INLINE
#endif

INLINE
CartesianIter::CartesianIter(int l):l_(l)
{
}

INLINE
CartesianIter::~CartesianIter()
{
}

INLINE void
CartesianIter::start()
{
  bfn_=a_=c_=0;
}

INLINE void
CartesianIter::next()
{
  if (c_<l_-a_) c_++; else {c_=0; a_++;} bfn_++;
}

INLINE
CartesianIter::operator int()
{
  return a_<=l_;
}

INLINE int
CartesianIter::a()
{
  return a_;
}

INLINE int
CartesianIter::b()
{
  return l_-a_-c_;
}

INLINE int
CartesianIter::c()
{
  return c_;
}

INLINE int
CartesianIter::l()
{
  return l_;
}

INLINE int
CartesianIter::bfn()
{
  return bfn_;
}

INLINE int
GaussianShell::nprimitive()
{
  return nprim;
}

INLINE int
GaussianShell::ncontraction()
{
  return ncon;
}

INLINE int
GaussianShell::nfunction()
{
  return nfunc;
}

INLINE int
GaussianShell::am(int con)
{
  return l[con];
}

INLINE char
GaussianShell::amchar(int con)
{
  return amtypes[l[con]];
}

INLINE int
GaussianShell::is_cartesian(int con)
{
  return !puream[con];
}

INLINE int
GaussianShell::is_pure(int con)
{
  return puream[con];
}

INLINE double
GaussianShell::coefficient_unnorm(int con,int prim)
{
  return coef[con][prim];
}

INLINE double
GaussianShell::normalization(int con,int bfn)
{
  return norm[con][bfn];
}

INLINE double
GaussianShell::exponent(int iprim)
{
  return exp[iprim];
}

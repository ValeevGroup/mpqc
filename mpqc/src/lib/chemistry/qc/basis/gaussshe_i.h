
#ifdef __GNUC__
#pragma interface
#endif

#ifdef INLINE_FUNCTIONS
#define INLINE inline
#else
#define INLINE
#endif

INLINE int
GaussianShell::nprimitive() const
{
  return nprim;
}

INLINE int
GaussianShell::ncontraction() const
{
  return ncon;
}

INLINE int
GaussianShell::nfunction() const
{
  return nfunc;
}

INLINE int
GaussianShell::am(int con) const
{
  return l[con];
}

INLINE char
GaussianShell::amchar(int con) const
{
  return amtypes[l[con]];
}

INLINE int
GaussianShell::is_cartesian(int con) const
{
  return !puream[con];
}

INLINE int
GaussianShell::is_pure(int con) const
{
  return puream[con];
}

INLINE double
GaussianShell::coefficient_unnorm(int con,int prim) const
{
  return coef[con][prim];
}

INLINE double
GaussianShell::exponent(int iprim) const
{
  return exp[iprim];
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// eval: (c-set-style "CLJ")

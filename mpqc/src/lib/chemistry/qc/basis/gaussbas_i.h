
#ifdef INLINE_FUNCTIONS
#define INLINE inline
#else
#define INLINE
#endif

INLINE const char*
GaussianBasisSet::name() const
{
  return name_;
}

INLINE int
GaussianBasisSet::nprimitive() const
{
  return nprim_;
}

INLINE int
GaussianBasisSet::nshell() const
{
  return nshell_;
}

INLINE int
GaussianBasisSet::nbasis() const
{
  return nbasis_;
}

INLINE int
GaussianBasisSet::shell_to_function(int i) const
{
  return shell_to_function_(i);
}

INLINE const GaussianShell&
GaussianBasisSet::operator()(int i) const
{
  return *shell[i];
}

INLINE const GaussianShell&
GaussianBasisSet::operator[](int i) const
{
  return *shell[i];
}

INLINE const double&
GaussianBasisSet::r(int icenter,int xyz) const
{
  return this->center_to_r_(icenter,xyz);
}

INLINE GaussianShell&
GaussianBasisSet::operator()(int i)
{
  return *shell[i];
}

INLINE GaussianShell&
GaussianBasisSet::operator[](int i)
{
  return *shell[i];
}

INLINE double&
GaussianBasisSet::r(int i,int xyz)
{
  return this->center_to_r_(i,xyz);
}

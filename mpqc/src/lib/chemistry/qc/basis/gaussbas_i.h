
#ifdef __GNUC__
#pragma interface
#endif

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
  return *shell_[i];
}

INLINE const GaussianShell&
GaussianBasisSet::operator[](int i) const
{
  return *shell_[i];
}

INLINE GaussianShell&
GaussianBasisSet::operator()(int i)
{
  return *shell_[i];
}

INLINE GaussianShell&
GaussianBasisSet::operator[](int i)
{
  return *shell_[i];
}

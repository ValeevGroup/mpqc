
#ifdef __GNUC__
#pragma implementation
#endif

#include <chemistry/qc/basis/dercent.h>

DerivCenters::DerivCenters()
{
  clear();
}

void
DerivCenters::clear()
{
  ncenter_ = 0;
  omitted_center_ = -1;
  omitted_atom_ = -1;
}

void
DerivCenters::add_center(int center, int atom)
{
  center_[ncenter_] = center;
  atom_[ncenter_] = atom;
  ncenter_++;
}

void
DerivCenters::add_omitted(int center, int atom)
{
  omitted_center_ = center;
  omitted_atom_ = atom;
  ncenter_++;
}

void
DerivCenters::add_center(int center, const RefGaussianBasisSet &b, int shell)
{
  add_center(center, b->shell_to_center(shell));
}

void
DerivCenters::add_omitted(int center, const RefGaussianBasisSet &b, int shell)
{
  add_omitted(center, b->shell_to_center(shell));
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// eval: (c-set-style "CLJ")

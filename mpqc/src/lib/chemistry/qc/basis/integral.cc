
#ifdef __GNUC__
#pragma implementation
#endif

#include <chemistry/qc/basis/integral.h>
#include <chemistry/qc/basis/shellrot.h>
#include <chemistry/qc/basis/petite.h>

SavableState_REF_def(Integral);

#define CLASSNAME Integral
#define PARENTS public SavableState
#include <util/state/statei.h>
#include <util/class/classia.h>
void *
Integral::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SavableState::_castdown(cd);
  return do_castdowns(casts,cd);
}

Integral::Integral()
{
}

Integral::Integral(StateIn& s) :
  SavableState(s)
{
}

Integral::Integral(const RefKeyVal&)
{
}

void
Integral::save_data_state(StateOut&)
{
}

RefPetiteList
Integral::petite_list(const RefGaussianBasisSet& gbs)
{
  return new PetiteList(gbs, *this);
}

ShellRotation
Integral::shell_rotation(int am, SymmetryOperation& so, int pure)
{
  ShellRotation r(am, so, *this, pure);
  return r;
}


#ifdef __GNUC__
#pragma implementation
#endif

#include <chemistry/qc/basis/integral.h>
#include <chemistry/qc/basis/shellrot.h>
#include <chemistry/qc/basis/petite.h>
#include <chemistry/qc/basis/obint.h>

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

Integral::Integral(const RefGaussianBasisSet &b1,
                   const RefGaussianBasisSet &b2,
                   const RefGaussianBasisSet &b3,
                   const RefGaussianBasisSet &b4)
{
  grp_ = MessageGrp::get_default_messagegrp();
  set_basis(b1,b2,b3,b4);
}

Integral::Integral(StateIn& s) :
  SavableState(s)
{
  bs1_.restore_state(s);
  bs2_.restore_state(s);
  bs3_.restore_state(s);
  bs4_.restore_state(s);
  s.get(storage_);
  grp_ = MessageGrp::get_default_messagegrp();
}

Integral::Integral(const RefKeyVal&)
{
  grp_ = MessageGrp::get_default_messagegrp();
}

void
Integral::save_data_state(StateOut&o)
{
  bs1_.save_state(o);
  bs2_.save_state(o);
  bs3_.save_state(o);
  bs4_.save_state(o);
  o.put(storage_);
}

RefPetiteList
Integral::petite_list()
{
  return new PetiteList(bs1_, this);
}

RefPetiteList
Integral::petite_list(const RefGaussianBasisSet& gbs)
{
  return new PetiteList(gbs, this);
}

ShellRotation
Integral::shell_rotation(int am, SymmetryOperation& so, int pure)
{
  this->reference();
  ShellRotation r(am, so, this, pure);
  this->dereference();
  return r;
}

void
Integral::set_basis(const RefGaussianBasisSet &b1,
                    const RefGaussianBasisSet &b2,
                    const RefGaussianBasisSet &b3,
                    const RefGaussianBasisSet &b4)
{
  bs1_ = b1;
  bs2_ = b2;
  bs3_ = b3;
  bs4_ = b4;
  if (bs2_.null()) bs2_ = bs1_;
  if (bs3_.null()) bs3_ = bs2_;
  if (bs4_.null()) bs4_ = bs3_;
}

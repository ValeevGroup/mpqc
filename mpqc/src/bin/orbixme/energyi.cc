
#include <strstream.h>
#include <iostream.h>
#include <energyi.h>

#include <chemistry/molecule/energy.h>

// force linkage of MPSCF
#include <chemistry/qc/mpqc/mpqc.h>
const ClassDesc &fl0 = MPSCF::class_desc_;

C_MolecularEnergyImpl::C_MolecularEnergyImpl()
{
}

C_MolecularEnergyImpl::~C_MolecularEnergyImpl()
{
}

MolecularEnergy *
C_MolecularEnergyImpl::mole()
{
  MolecularEnergy *ret;
  ret = MolecularEnergy::castdown(dc_);
  return ret;
}

double
C_MolecularEnergyImpl::energy(CORBA_Environment &)
{
  if (!mole()) return 0.0;
  return mole()->energy();
}



#include <energyi.h>

#include <chemistry/molecule/energy.h>

// force linkages

#include <chemistry/qc/mpqc/mpqc.h>
const ClassDesc &fl0 = MPSCF::class_desc_;
#include <chemistry/molecule/coor.h>
const ClassDesc &fl1 = SymmMolecularCoor::class_desc_;
const ClassDesc &fl2 = RedundMolecularCoor::class_desc_;

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

C_Molecule *
C_MolecularEnergyImpl::molecule(CORBA_Environment &e)
{
  Molecule *m = 0;
  if (mole()) m = mole()->molecule().pointer();
  TIE_C_Molecule(C_MoleculeImpl) *c_m = 0;

  if (m) {
      c_m = new TIE_C_Molecule(C_MoleculeImpl)
            (new C_MoleculeImpl(m), /*marker*/ 0, /*loader*/ 0);
    }

  return c_m;
}

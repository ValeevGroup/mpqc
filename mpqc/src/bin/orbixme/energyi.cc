
#include <energyi.h>

#include <chemistry/molecule/energy.h>

// force linkages

#include <chemistry/qc/mpqc/mpqc.h>
static ForceLink<MPSCF> fl0;
#include <chemistry/molecule/coor.h>
static ForceLink<SymmMolecularCoor> fl1;
static ForceLink<RedundMolecularCoor> fl2;

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
  ret = dynamic_cast<MolecularEnergy*>(dc_);
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

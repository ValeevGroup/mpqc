
#ifndef energyi_h
#define energyi_h

#include "keyvali.h"
#include "energy.h"

class MolecularEnergy;

class C_MolecularEnergyImpl: public C_KeyValCreatableImpl {
  protected:
    MolecularEnergy *mole();
  public:
    C_MolecularEnergyImpl();
    ~C_MolecularEnergyImpl();

    double energy(CORBA_Environment &IT_env);
};

DEF_TIE_C_MolecularEnergy(C_MolecularEnergyImpl);

#endif

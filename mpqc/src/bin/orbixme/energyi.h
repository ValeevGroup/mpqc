
#ifndef energyi_h
#define energyi_h

#include "functioni.h"
#include "moleculei.h"
#include "energy.h"

class MolecularEnergy;

class C_MolecularEnergyImpl: public C_FunctionImpl {
  protected:
    MolecularEnergy *mole();
  public:
    C_MolecularEnergyImpl();
    ~C_MolecularEnergyImpl();

    double energy(CORBA_Environment &IT_env);
    C_Molecule *molecule(CORBA_Environment &e);

    unsigned char molecularenergy_has_object(CORBA_Environment &e)
    { return has_object(e); }
};

DEF_TIE_C_MolecularEnergy(C_MolecularEnergyImpl);

#endif

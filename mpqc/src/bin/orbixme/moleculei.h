
#ifndef moleculei_h
#define moleculei_h

#include "keyvali.h"
#include "molecule.h"

class Molecule;

class C_MoleculeImpl: public C_KeyValCreatableImpl {
  protected:
    Molecule *mol();
  public:
    C_MoleculeImpl();
    ~C_MoleculeImpl();

    long natom(CORBA_Environment &IT_env);
};

DEF_TIE_C_Molecule(C_MoleculeImpl);

#endif

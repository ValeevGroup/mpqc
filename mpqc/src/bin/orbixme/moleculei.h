
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
    C_MoleculeImpl(Molecule*);
    ~C_MoleculeImpl();

    long natom(CORBA_Environment &IT_env);
    long atomic_number(long atomnum, CORBA_Environment &);
    double x(long atomnum, CORBA_Environment &);
    double y(long atomnum, CORBA_Environment &);
    double z(long atomnum, CORBA_Environment &);
    double r(long atom1, long atom2, CORBA_Environment &);

    unsigned char molecule_has_object(CORBA_Environment &e)
    { return has_object(e); }
};

DEF_TIE_C_Molecule(C_MoleculeImpl);

#endif

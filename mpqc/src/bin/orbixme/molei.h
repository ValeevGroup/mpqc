
#ifndef molei_h
#define molei_h

#include "mole.h"

class MolecularEnergy;

class MolEImpl: public MolEBOAImpl {
    MolecularEnergy *mole_;
public:
    MolEImpl(const char *, CORBA(LoaderClass)*);
    ~MolEImpl();

    void create(const char *, CORBA_Environment &IT_env);
    double energy(CORBA_Environment &IT_env);
};

#endif

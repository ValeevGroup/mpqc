
#ifndef optimizei_h
#define optimizei_h

#include "keyvali.h"
#include "optimize.h"

class Optimize;

class C_OptimizeImpl: public C_KeyValCreatableImpl {
  protected:
    Optimize *opt();
  public:
    C_OptimizeImpl();
    ~C_OptimizeImpl();

    long optimize(CORBA_Environment &IT_env);
};

DEF_TIE_C_Optimize(C_OptimizeImpl);

#endif

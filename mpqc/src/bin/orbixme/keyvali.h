
#ifndef keyvali_h
#define keyvali_h

#include "keyval.h"

class DescribedClass;

class C_KeyValCreatableImpl {
  protected:
    DescribedClass *dc_;
  public:
    C_KeyValCreatableImpl();
    virtual ~C_KeyValCreatableImpl();

    virtual void keyval_create(const char *, CORBA_Environment &IT_env);
};

DEF_TIE_C_KeyValCreatable(C_KeyValCreatableImpl);

#endif


#ifndef keyvali_h
#define keyvali_h

#include "keyval.h"

class DescribedClass;

class C_KeyValCreatableImpl {
  protected:
    DescribedClass *dc_;

    void clear_dc();
    void set_dc(DescribedClass *);
  public:
    C_KeyValCreatableImpl();
    C_KeyValCreatableImpl(DescribedClass *);
    virtual ~C_KeyValCreatableImpl();

    virtual void keyval_create(const char *, CORBA_Environment &IT_env);
    unsigned char has_object(CORBA_Environment &IT_env);

    DescribedClass *object() { return dc_; }
};

DEF_TIE_C_KeyValCreatable(C_KeyValCreatableImpl);

#endif

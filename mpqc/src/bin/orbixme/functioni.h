
#ifndef functioni_h
#define functioni_h

#include "keyvali.h"
#include "function.h"

class NLP2;

class C_FunctionImpl: public C_KeyValCreatableImpl {
  protected:
    NLP2 *func();
  public:
    C_FunctionImpl();
    ~C_FunctionImpl();

    double value(CORBA_Environment &IT_env);

    unsigned char function_has_object(CORBA_Environment &e)
    { return has_object(e); }
};

DEF_TIE_C_Function(C_FunctionImpl);

#endif

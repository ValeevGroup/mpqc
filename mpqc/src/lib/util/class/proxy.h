
#ifndef _util_keyval_proxy_h
#define _util_keyval_proxy_h

#include <util/class/class.h>

class DescribedClassProxy: public DescribedClass {
#   define CLASSNAME DescribedClassProxy
#   include <util/class/classda.h>
  public:
    virtual RefDescribedClass object() = 0;
};

#endif

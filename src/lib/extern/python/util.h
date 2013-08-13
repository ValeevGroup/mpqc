
#ifndef _util_h
#define _util_h

#include <util/ref/ref.h>

namespace sc {
  template <class T>
    inline T*
    get_pointer(const sc::Ref<T>& ref) {
    return ref.pointer();
  }
}

#endif

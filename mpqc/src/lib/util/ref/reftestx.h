
#ifdef __GNUG__
#pragma interface
#endif

#include <util/ref/ref.h>

class X: public VRefCount {
  private:
    int x;
  public:
    static int nx;
    X();
    ~X();
};

Ref_declare(X);


#include <util/misc/scexception.h>

using namespace sc;

class X: public DescribedClass {
  public:
    X() {};
    void x() {
      throw AlgorithmException("algorithm exception in X::x()",
                               __FILE__,
                               __LINE__,
                               this->class_desc());
    }
};

static ClassDesc X_cd(typeid(X),"X",1,"public DescribedClass");

void
f()
{
  throw AlgorithmException("algorithm exception in f()",
                           __FILE__,
                           __LINE__
                           );
}

void
g()
{
  throw MaxIterExceeded("maxiter exception in g()",
                        __FILE__,
                        __LINE__,
                        10);
}

void
h()
{
  throw ToleranceExceeded("tolerance exception in g()",
                          __FILE__,
                          __LINE__,
                          1.0, 2.0);
}

main()
{
  try {
      f();
      std::cout << "f() ran OK" << std::endl;
    }
  catch (SCException &e) {
      std::cout << "got an f() exception" << std::endl;
      std::cout << e.what() << std::endl;
    }

  try {
      g();
      std::cout << "g() ran OK" << std::endl;
    }
  catch (SCException &e) {
      std::cout << "got an g() exception" << std::endl;
      std::cout << e.what() << std::endl;
    }

  try {
      h();
      std::cout << "h() ran OK" << std::endl;
    }
  catch (SCException &e) {
      std::cout << "got an h() exception" << std::endl;
      std::cout << e.what() << std::endl;
    }

  try {
      X x;
      x.x();
      std::cout << "x.x() ran OK" << std::endl;
    }
  catch (SCException &e) {
      std::cout << "got an x.x() exception" << std::endl;
      std::cout << e.what() << std::endl;
    }

  return 0;
}



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

void
i()
{
  throw SystemException("system exception in i()",
                           __FILE__,
                           __LINE__
                           );
}

void
j()
{
  throw MemAllocFailed("memalloc exception in j()",
                        __FILE__,
                        __LINE__,
                        100000);
}

void
k()
{
  throw ProgrammingError("programming error in k()",
                         __FILE__,
                         __LINE__);
}

void
l()
{
  throw InputError("input error in l()",
                   __FILE__,
                   __LINE__,
                   "the_keyword");
}

void
m()
{
  throw LimitExceeded<int>("limit exceeded in m()",
                           __FILE__,
                           __LINE__,
                           10, 11);
}

void
n()
{
  throw FileOperationFailed("file op failed in n()",
                            __FILE__,
                            __LINE__,
                            "outfile",
                            FileOperationFailed::OpenRW);
}

void
o()
{
  throw SyscallFailed("syscall failed in o()",
                      __FILE__,
                      __LINE__,
                      "xyz",
                      10);
}

main()
{
  try {
      f();
      std::cout << "ERROR: f() ran OK" << std::endl;
    }
  catch (SCException &e) {
      std::cout << "EXPECTED: got an f() exception" << std::endl;
      std::cout << e.what() << std::endl;
    }

  try {
      g();
      std::cout << "ERROR: g() ran OK" << std::endl;
    }
  catch (SCException &e) {
      std::cout << "EXPECTED: got an g() exception" << std::endl;
      std::cout << e.what() << std::endl;
    }

  try {
      h();
      std::cout << "ERROR: h() ran OK" << std::endl;
    }
  catch (SCException &e) {
      std::cout << "EXPECTED: got an h() exception" << std::endl;
      std::cout << e.what() << std::endl;
    }

  try {
      i();
      std::cout << "ERROR: i() ran OK" << std::endl;
    }
  catch (SCException &e) {
      std::cout << "EXPECTED: got an i() exception" << std::endl;
      std::cout << e.what() << std::endl;
    }

  try {
      j();
      std::cout << "ERROR: j() ran OK" << std::endl;
    }
  catch (SCException &e) {
      std::cout << "EXPECTED: got an j() exception" << std::endl;
      std::cout << e.what() << std::endl;
    }

  try {
      k();
      std::cout << "ERROR: k() ran OK" << std::endl;
    }
  catch (SCException &e) {
      std::cout << "EXPECTED: got an k() exception" << std::endl;
      std::cout << e.what() << std::endl;
    }

  try {
      l();
      std::cout << "ERROR: l() ran OK" << std::endl;
    }
  catch (SCException &e) {
      std::cout << "EXPECTED: got an l() exception" << std::endl;
      std::cout << e.what() << std::endl;
    }

  try {
      m();
      std::cout << "ERROR: m() ran OK" << std::endl;
    }
  catch (SCException &e) {
      std::cout << "EXPECTED: got an m() exception" << std::endl;
      std::cout << e.what() << std::endl;
    }

  try {
      n();
      std::cout << "ERROR: n() ran OK" << std::endl;
    }
  catch (SCException &e) {
      std::cout << "EXPECTED: got an n() exception" << std::endl;
      std::cout << e.what() << std::endl;
    }

  try {
      o();
      std::cout << "ERROR: o() ran OK" << std::endl;
    }
  catch (SCException &e) {
      std::cout << "EXPECTED: got an o() exception" << std::endl;
      std::cout << e.what() << std::endl;
    }

  try {
      X x;
      x.x();
      std::cout << "x.x() ran OK" << std::endl;
    }
  catch (SCException &e) {
      std::cout << "EXPECTED: got an x.x() exception" << std::endl;
      std::cout << e.what() << std::endl;
    }

  return 0;
}


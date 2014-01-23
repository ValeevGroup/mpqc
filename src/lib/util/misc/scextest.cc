
#include <util/misc/scexception.h>
#include <stdexcept>

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

class NestedException: public SCException {
  public:
    NestedException(
        const char *description = 0,
        const char *file = 0,
        int line = 0,
        const ClassDesc *class_desc = 0,
        const char *exception_type = "NestedException") MPQC__NOEXCEPT:
      SCException(description, file, line, class_desc, exception_type)
        {
          try {
              if (description) throw NestedException();
              throw std::runtime_error("ctor");
            }
          catch (...) {}
        }
  
    NestedException(const NestedException& ref) MPQC__NOEXCEPT:
      SCException(ref)
        {
          try {
              if (description()) throw NestedException();
              throw std::runtime_error("ctor");
            }
          catch (...) {}
        }

    ~NestedException() MPQC__NOEXCEPT
        {
          try {
              if (description()) throw NestedException();
              throw std::runtime_error("ctor");
            }
          catch (...) {}
        }
};

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
                   "the_keyword",
                   "the_value");
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

void
p()
{
  throw FeatureNotImplemented("p() not implemented",
                              __FILE__,
                              __LINE__);
}

void
ex_on_stack()
{
  ProgrammingError ex("programming error in ex_on_stack()",
                      __FILE__, __LINE__);
  try {
      ex.elaborate() << "more info about the problem" << std::endl;
      throw std::runtime_error("whoops");
    }
  catch (...) {}
  throw ex;
}

void
nested()
{
  throw NestedException("nested exception test",
                        __FILE__, __LINE__);
}

int main(int argc, char* argv[])
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
      p();
      std::cout << "ERROR: p() ran OK" << std::endl;
    }
  catch (SCException &e) {
      std::cout << "EXPECTED: got an p() exception" << std::endl;
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

  try {
      ex_on_stack();
      std::cout << "ERROR: ex_on_stack() ran OK" << std::endl;
    }
  catch (SCException &e) {
      std::cout << "EXPECTED: got an ex_on_stack() exception" << std::endl;
      std::cout << e.what() << std::endl;
    }

  try {
      nested();
      std::cout << "ERROR: nested() ran OK" << std::endl;
    }
  catch (SCException &e) {
      std::cout << "EXPECTED: got an nested() exception" << std::endl;
      std::cout << e.what() << std::endl;
    }

  return 0;
}


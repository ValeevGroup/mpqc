
#include "catch.hpp"
#include "mpqc/util/core/exception.h"

using namespace mpqc;

class X {
  public:
    X() {};
    void x() {
      throw AlgorithmException("algorithm exception in X::x()",
                               __FILE__,
                               __LINE__);
    }
};

class NestedException: public Exception {
  public:
    NestedException(
        const char *description = 0,
        const char *file = 0,
        int line = 0,
        const char *exception_type = "NestedException") MPQC__NOEXCEPT:
      Exception(description, file, line, exception_type)
        {
          try {
              if (description) throw NestedException();
              throw std::runtime_error("ctor");
            }
          catch (...) {}
        }
  
    NestedException(const NestedException& ref) MPQC__NOEXCEPT:
      Exception(ref)
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
q()
{
  throw FeatureDisabled("q() disabled",
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

TEST_CASE("Exception", "[exception]") {
  REQUIRE_THROWS_AS(f(), mpqc::Exception);
  REQUIRE_THROWS_AS(g(), mpqc::Exception);
  REQUIRE_THROWS_AS(h(), mpqc::Exception);
  REQUIRE_THROWS_AS(i(), mpqc::Exception);
  REQUIRE_THROWS_AS(j(), mpqc::Exception);
  REQUIRE_THROWS_AS(k(), mpqc::Exception);
  REQUIRE_THROWS_AS(l(), mpqc::Exception);
  REQUIRE_THROWS_AS(m(), mpqc::Exception);
  REQUIRE_THROWS_AS(n(), mpqc::Exception);
  REQUIRE_THROWS_AS(o(), mpqc::Exception);
  REQUIRE_THROWS_AS(p(), mpqc::Exception);
  REQUIRE_THROWS_AS(q(), mpqc::Exception);

  X x;
  REQUIRE_THROWS_AS(x.x(), mpqc::Exception);

  REQUIRE_THROWS_AS(ex_on_stack(), mpqc::Exception);

  REQUIRE_THROWS_AS(nested(), mpqc::Exception);

};  // TEST_CASE

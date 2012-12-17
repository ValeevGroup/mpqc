
#include <iostream>
#include <util/misc/autovec.h>

using namespace sc;

class A {
    int me;
    static int nA;
  public:
    A() {
      me = nA++;
      std::cout << "A::A (" << me << ")" << std::endl;
    }
    ~A() {
      std::cout << "A::~A (" << me << ")" << std::endl;
    }
};

int A::nA = 0;

void
expect(const std::string &m)
{
  std::cout << "expect " << m << ": ";
}

int
main(int argc, char* argv[])
{
  auto_vec<double> double_data(new double[4]);
  double_data.release();

  {
    expect("A::A (0)");
    auto_vec<A> a(new A[1]);
    auto_vec<A> b(a);
    expect("A::~A (0)");
    b.reset();
  }

  {
    expect("A::A (1)");
    auto_vec<A> a(new A[1]);
    expect("A::A (2)");
    auto_vec<A> b(new A[1]);
    expect("A::~A (2)");
    b = a;
    expect("A::~A (1)");
    b.reset();
  }

  {
    expect("A::A (3)");
    auto_vec<A> a(new A[1]);
    expect("A::~A (3)");
  }

  {
    expect("A::A (4)");
    auto_vec<A> a(new A[1]);
    auto_vec<A> b(a.release());
    expect("A::~A (4)");
  }

  {
    expect("A::A (5)");
    auto_vec<A> a(new A[1]);
    expect("A::A (6)");
    A *ap = new A[1];
    expect("A::~A (5)");
    a.reset(ap);
    expect("A::~A (6)");
  }

  return 0;
}

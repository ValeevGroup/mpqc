
//compilation:
// g++ -g    bug.cc   -o bug
//output on SGI IRIX 4.0.5F with gcc 2.5.{0,2}:
// Aref(B&) OK
// Aref(C&) OK
// Aval(B) FAILED
// Aval(C) OK
//works on SGI IRIX 4.0.5F with gcc 2.4.5
//works on sparc SunOS 4.1.1 with gcc 2.5.2

extern "C" {
    int printf(const char*, ...);
}

class B {
  public:
    B() {}
    B(const B&) {}
    ~B() {}
};

// C is the same as B, except it has no user defined dtor
class C {
  public:
    C() {}
    C(const C&) {}
};

class Aval {
  public:
    Aval(B) {};
    Aval(C) {};
};

class Aref {
  public:
    Aref(B&) {};
    Aref(C&) {};
};

main()
{
  if (new Aref(B())==0) printf("Aref(B&) FAILED\n");
  else printf("Aref(B&) OK\n");

  if (new Aref(C())==0) printf("Aref(C&) FAILED\n");
  else printf("Aref(C&) OK\n");

  if (new Aval(B())==0) printf("Aval(B) FAILED\n");
  else printf("Aval(B) OK\n");

  if (new Aval(C())==0) printf("Aval(C) FAILED\n");
  else printf("Aval(C) OK\n");
}

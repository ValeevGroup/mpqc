
#include <stdio.h>
#include <iostream.h>
#include <util/container/ref.h>
#include <util/keyval/keyval.h>
#include <util/class/class.h>
#include <util/state/state.h>
#include <util/group/message.h>
#include <util/group/mstate.h>

#define A_parents virtual public SavableState
class A: A_parents {
#   define CLASSNAME A
#   define HAVE_CTOR
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  private:
    int ia;
    int n;
    int* array;
    double d;
  public:
    A(int size);
    A(const RefKeyVal&);
    A(StateIn&);
    ~A();
    void save_data_state(StateOut&);
    inline int& a() { return ia; };
    virtual void print (ostream&s = cout)
    {
      s << "A::a = " << a() << '\n';
      s << "A::array = {";
      for (int i=0; i<n; i++) s << array[i] << ' ';
      s << "}\n";
    }
};
SavableState_REF_dec(A);
SavableState_REF_def(A);
A::A(int size):
  ia(1),
  n(size),
  array(new int[size]),
  d(-1.24)
{
  for (int i=0; i<size; i++) array[i] = size - i - 1;
}
A::A(const RefKeyVal&keyval):
  ia(keyval->intvalue("a")),
  n(keyval->intvalue("n")),
  d(-1.24)

{
  array = new int[n];
  for (int i=0; i<n; i++) array[i] = i + 10000;
}
A::A(StateIn&s):
  SavableState(s,A::class_desc_)
{
  printf("getting d\n"); fflush(stdout);
  s.get(d);
  printf("getting ia\n"); fflush(stdout);
  s.get(ia);
  printf("getting array\n"); fflush(stdout);
  s.get(n);
  s.get(array);
  printf("got everything\n"); fflush(stdout);
}
A::~A()
{
  delete[] array;
}
void
A::save_data_state(StateOut&s)
{
  printf("putting d\n"); fflush(stdout);
  s.put(d);
  printf("putting ia\n"); fflush(stdout);
  s.put(ia);
  printf("putting array\n"); fflush(stdout);
  s.put(n);
  s.put(array,n);
  printf("put everything\n"); fflush(stdout);
}

#define CLASSNAME A
#define PARENTS A_parents
#define HAVE_CTOR
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
A::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SavableState::_castdown(cd);
  return do_castdowns(casts,cd);
}

void test(MessageGrp&, int source, int target);

int
main()
{
  //ProcMessageGrp grp;
  //test(grp, 0, 0);

  ShmMessageGrp grp(4);
  test(grp, 2, 1);

  return 0;
}

void
test(MessageGrp& grp, int source, int target)
{
  RefA a,b;
  
  if (grp.me() == source) {
      StateSend so(grp);
      so.set_buffer_size(5);
      so.target(target);
      a = new A(10);
      a.save_state(so);
      so.flush();
    }

  if (grp.me() == target) {
      StateRecv si(grp);
      si.set_buffer_size(5);
      si.source(source);
      b.restore_state(si);
    }

  if (grp.me() == target) {
      printf("target:\n");
      b->print();
    }

  grp.sync();

  if (grp.me() == source) {
      printf("source:\n");
      a->print();
    }

  ///////////////////////////////////////////////////
  // Test broadcast
  grp.sync();

  b = 0;
  
  if (grp.me() == source) {
      BcastStateSend so(grp);
      a.save_state(so);
    }
  else {
      BcastStateRecv si(grp,source);
      b.restore_state(si);
    }

  if (grp.me() == target) {
      printf("bcast target:\n");
      b->print();
    }

}

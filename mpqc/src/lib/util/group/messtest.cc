
#include <stdio.h>
#include <iostream.h>
#include <util/container/ref.h>
#include <util/keyval/keyval.h>
#include <util/class/class.h>
#include <util/state/state.h>
#include <util/group/message.h>
#include <util/group/mstate.h>
#include <util/group/hcube.h>

// Force linkages:
#ifndef __PIC__
# ifndef PARAGON
#   include <util/group/messshm.h>
    const ClassDesc &fl0 = ShmMessageGrp::class_desc_;
# endif
# ifdef PARAGON
#   include <util/group/messpgon.h>
    const ClassDesc &fl1 = ParagonMessageGrp::class_desc_;
# endif
#endif

#define A_parents virtual_base public SavableState
class A: A_parents {
#   define CLASSNAME A
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
  SavableState(s)
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

void test(const RefMessageGrp&, int source, int target);
void test_hcube(int nproc, int root, int fwd);

int
main(int argc, char**argv)
{
  RefMessageGrp grp;

  const char* input = SRCDIR "/messtest.in";
  const char* keyword = "message";

  if (argc >= 2) input = argv[1];
  if (argc >= 3) keyword = argv[2];

  RefKeyVal keyval = new ParsedKeyVal(input);

  grp = keyval->describedclassvalue(keyword);

  if (grp.null()) {
      fprintf(stderr,"Couldn't initialize MessageGrp\n");
      abort();
    }

  grp->sync();

  if (grp->me() == 0) {
      test_hcube(39, 0, 1);
      test_hcube(16, 0, 1);
      test_hcube(17, 4, 1);
      test_hcube(17, 4, 0);
      test_hcube(1, 0, 0);
    }

  grp->sync();

  if (grp->n() >= 3) {
      test(grp, 2, 1);
    }
  else {
      test(grp, 0, 0);
    }

  grp->sync();
  return 0;
}

void
test_hcube(int nproc, int root, int fwd)
{
  int i, j;
  RefGlobalMsgIter *gmi = new RefGlobalMsgIter[nproc];
  for (i=0; i<nproc; i++) {
      gmi[i] =  new HypercubeGMI(nproc, i, root);
    }
  int iter = 1;
  for (j=0; j<nproc; j++) {
      if (fwd) {
          gmi[j]->forwards();
        }
      else {
          gmi[j]->backwards();
        }
    }
  while (!gmi[0]->done()) {
      printf("------ step %d of %d ------\n", iter, gmi[0]->n());
      for (j=0; j<nproc; j++) {
          if (gmi[j]->send()) {
              if (0 <= gmi[j]->sendto() && gmi[j]->sendto() < nproc) {
                  if (gmi[gmi[j]->sendto()]->recvfrom() == j) {
                      printf(" %d -> %d\n", j, gmi[j]->sendto());
                    }
                  else {
                      printf(" %d -> (%d)\n", j, gmi[j]->sendto());
                    }
                }
              else {
                  printf(" %d -> %d?\n", j, gmi[j]->sendto());
                }
            }
          else if (gmi[j]->recv()) {
              if (0 <= gmi[j]->recvfrom() && gmi[j]->recvfrom() < nproc) {
                  if (gmi[gmi[j]->recvfrom()]->sendto() == j) {
                      // to be printed by sender
                    }
                  else {
                      printf(" (%d) -> %d\n", gmi[j]->recvfrom(), j);
                    }
                }
              else {
                  printf(" %d? -> %d\n", gmi[j]->recvfrom(), j);
                }
            }
        }
      for (j=0; j<nproc; j++) gmi[j]->next();
      iter++;
    }
  fflush(stdout);
}

void
test(const RefMessageGrp& grp, int source, int target)
{
  RefA a,b;
  
  if (grp->me() == source) {
      StateSend so(grp);
      so.set_buffer_size(5);
      so.target(target);
      a = new A(10);
      a.save_state(so);
      so.flush();
    }

  if (grp->me() == target) {
      StateRecv si(grp);
      si.set_buffer_size(5);
      si.source(source);
      b.restore_state(si);
    }

  if (grp->me() == target) {
      printf("target:\n");
      b->print();
    }

  grp->sync();

  if (grp->me() == source) {
      printf("source:\n");
      a->print();
    }

  ///////////////////////////////////////////////////
  // Test broadcast

  if (source != target) {
      grp->sync();

      b = 0;
  
      if (grp->me() == source) {
          BcastStateSend so(grp);
          a.save_state(so);
        }
      else {
          BcastStateRecv si(grp,source);
          b.restore_state(si);
        }

      if (grp->me() == target) {
          printf("bcast target:\n");
          b->print();
        }
    }

}

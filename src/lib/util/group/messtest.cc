//
// messtest.cc
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@limitpt.com>
// Maintainer: LPS
//
// This file is part of the SC Toolkit.
//
// The SC Toolkit is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The SC Toolkit is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the SC Toolkit; see the file COPYING.LIB.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//
// The U.S. Government is granted a limited license as per AL 91-7.
//

#include <util/misc/formio.h>
#include <util/keyval/keyval.h>
#include <util/class/class.h>
#include <util/state/state.h>
#include <util/misc/bug.h>
#include <util/group/message.h>
#include <util/group/mstate.h>
#include <util/group/hcube.h>

using namespace std;
using namespace sc;

// Force linkages:
//#ifndef __PIC__
# ifdef HAVE_MPI
#   include <util/group/messmpi.h>
    static ForceLink<MPIMessageGrp> fl2;
# endif
//#endif

class A: virtual public SavableState {
  private:
    int ia;
    int n;
    int* array;
    double d;
  public:
    A(int size);
    A(const Ref<KeyVal>&);
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

A::A(int size):
  ia(1),
  n(size),
  array(new int[size]),
  d(-1.24)
{
  for (int i=0; i<size; i++) array[i] = size - i - 1;
}
A::A(const Ref<KeyVal>&keyval):
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
  cout << "getting d" << endl;
  s.get(d);
  cout << "getting ia" << endl;
  s.get(ia);
  cout << "getting array" << endl;
  s.get(n);
  s.get(array);
  cout << "got everything" << endl;
}
A::~A()
{
  delete[] array;
}
void
A::save_data_state(StateOut&s)
{
  cout << "putting d" << endl;
  s.put(d);
  cout << "putting ia" << endl;
  s.put(ia);
  cout << "putting array" << endl;
  s.put(n);
  s.put(array,n);
  cout << "put everything" << endl;
}

static ClassDesc A_cd(
  typeid(A),"A",1,"virtual public SavableState",
  0, create<A>, create<A>);

void test(const Ref<MessageGrp>&, int source, int target);
void test_hcube(int nproc, int root, int fwd);

int
main(int argc, char**argv)
{
  Ref<MessageGrp> grp = MessageGrp::initial_messagegrp(argc, argv);

  Ref<Debugger> debugger;

  if (grp == 0) {
      const char* input = SRCDIR "/messtest.in";
      const char* keyword = "message";

      if (argc >= 2) input = argv[1];
      if (argc >= 3) keyword = argv[2];

      Ref<KeyVal> keyval = new ParsedKeyVal(input);

      grp << keyval->describedclassvalue(keyword);

      debugger << keyval->describedclassvalue(":debug");

      if (grp == 0) {
          cerr << scprintf("Couldn't initialize MessageGrp\n");
          abort();
        }
    }

  if (debugger) {
      debugger->set_exec(argv[0]);
      debugger->set_prefix(grp->me());
    }

  Debugger::set_default_debugger(debugger);

  grp->sync();
  if (grp->n() > 1) {
      BcastState bc(grp,1);
      bc.bcast(debugger);
      bc.flush();
    }
  grp->sync();
  if (debugger) {
      debugger->set_exec(argv[0]);
      debugger->set_prefix(grp->me());
      debugger->traceback();
    }
  grp->sync();

  if (0 && grp->me() == 0) {
      test_hcube(3, 0, 1);
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
  else if (grp->n() >= 2) {
      test(grp, 1, 0);
    }
  else {
      test(grp, 0, 0);
    }

  int testsum = 1;
  grp->sum(&testsum,1);
  if (testsum != grp->n()) {
      cerr << scprintf("WARNING: sum wrong\n");
    }

  double testdsum = 1.0;
  grp->sum(&testdsum,1);
  cout << scprintf("on %d testdsum = %4.1f\n", grp->me(), testdsum);

  grp->sync();
  grp = 0;
  MessageGrp::set_default_messagegrp(0);
  return 0;
}

void
test_hcube(int nproc, int root, int fwd)
{
  int i, j;
  Ref<GlobalMsgIter> *gmi = new Ref<GlobalMsgIter>[nproc];
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
      cout << scprintf("------ step %d of %d ------\n", iter, gmi[0]->n());
      for (j=0; j<nproc; j++) {
          if (gmi[j]->send()) {
              if (0 <= gmi[j]->sendto() && gmi[j]->sendto() < nproc) {
                  if (gmi[gmi[j]->sendto()]->recvfrom() == j) {
                      cout << scprintf(" %d -> %d\n", j, gmi[j]->sendto());
                    }
                  else {
                      cout << scprintf(" %d -> (%d)\n", j, gmi[j]->sendto());
                    }
                }
              else {
                  cout << scprintf(" %d -> %d?\n", j, gmi[j]->sendto());
                }
            }
          else if (gmi[j]->recv()) {
              if (0 <= gmi[j]->recvfrom() && gmi[j]->recvfrom() < nproc) {
                  if (gmi[gmi[j]->recvfrom()]->sendto() == j) {
                      // to be printed by sender
                    }
                  else {
                      cout << scprintf(" (%d) -> %d\n", gmi[j]->recvfrom(), j);
                    }
                }
              else {
                  cout << scprintf(" %d? -> %d\n", gmi[j]->recvfrom(), j);
                }
            }
        }
      for (j=0; j<nproc; j++) gmi[j]->next();
      iter++;
    }
  cout.flush();
}

void
test(const Ref<MessageGrp>& grp, int source, int target)
{
  Ref<A> a,b;
  const int nca = 1000000;
  char ca[nca];
  
  if (grp->me() == source) {
      StateSend so(grp);
      //so.set_buffer_size(5);
      so.target(target);
      a = new A(10);
      SavableState::save_state(a,so);
      so.flush();
      grp->send(target, ca, nca);
      grp->recv(target, ca, nca);
    }

  if (grp->me() == target) {
      StateRecv si(grp);
      //si.set_buffer_size(5);
      si.source(source);
      b << SavableState::restore_state(si);
      grp->recv(source, ca, nca);
      grp->send(source, ca, nca);
    }

  if (grp->me() == target) {
      cout << "target:" << endl;
      b->print();
    }

  grp->sync();

  if (grp->me() == source) {
      cout << "source:" << endl;
      a->print();
    }

  ///////////////////////////////////////////////////
  // Test broadcast

  if (source != target) {
      grp->sync();

      b = 0;
  
      if (grp->me() == source) {
          BcastStateSend so(grp);
          SavableState::save_state(a,so);
        }
      else {
          BcastStateRecv si(grp,source);
          b << SavableState::restore_state(si);
        }

      if (grp->me() == target) {
          cout << "bcast target:" << endl;
          b->print();
        }
    }

}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:

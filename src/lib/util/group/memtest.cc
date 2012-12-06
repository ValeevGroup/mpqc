//
// memtest.cc
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

#include <math.h>
#include <util/misc/formio.h>
#include <util/misc/bug.h>
#include <util/group/message.h>
#include <util/group/mstate.h>
#include <util/group/hcube.h>
#include <util/group/memshm.h>
#include <util/group/memregion.h>
#ifdef HAVE_NX
#  include <util/group/memipgon.h>
#endif

using namespace std;
using namespace sc;

// Force linkages:
//#ifndef __PIC__
#ifdef HAVE_MPI
#   include <util/group/messmpi.h>
#   include <util/group/memmtmpi.h>
    static ForceLink<MPIMessageGrp> fl2;
    static ForceLink<MTMPIMemoryGrp> fl3;
#endif
//#endif

#include <util/misc/scexception.h>
static const char * (sc::SCException::*force_except_link)() const
    = &sc::SCException::description;

// this is needed for debugging
#ifdef HAVE_NX
extern "C" {
#include <nx.h>
}
#endif // HAVE_NX

#ifdef HAVE_HRECV
#  define DISABLE do { masktrap(1); cout.flush(); } while(0)
#  define ENABLE do { cout.flush(); masktrap(0); } while(0)
   extern "C" {
       long masktrap(long state);
     }
#else
#  define DISABLE
#  define ENABLE
#endif

#define PRINTF(args) do { DISABLE; \
                          cout << scprintf args; \
                          cout.flush(); \
                          ENABLE; \
                         } while(0)

void do_simple_tests(const Ref<MessageGrp>&,const Ref<MemoryGrp>&);
void do_int_tests(const Ref<MessageGrp>&,const Ref<MemoryGrp>&);
void do_double_tests(const Ref<MessageGrp>&,const Ref<MemoryGrp>&);
void do_double2_tests(const Ref<MessageGrp>&,const Ref<MemoryGrp>&);
void do_memgrp_tests(const Ref<MessageGrp>&,const Ref<MemoryGrp>&);
void do_memgrpregion_tests(const Ref<MessageGrp>&,const Ref<MemoryGrp>&);

int
main(int argc, char**argv)
{
  Ref<MessageGrp> msg = MessageGrp::initial_messagegrp(argc, argv);

  const char* input = SRCDIR "/memtest.in";
  Ref<KeyVal> keyval = new ParsedKeyVal(input);

  if (msg.null()) {
      const char* keyword = "message";

      if (argc >= 2) input = argv[1];
      if (argc >= 3) keyword = argv[2];

      msg << keyval->describedclassvalue(keyword);

      if (msg.null()) {
          cerr << scprintf("Couldn't initialize MessageGrp\n");
          abort();
        }
    }

// This causes problems for automated testing:
//   // now set up the debugger
//   Ref<Debugger> debugger;
//   debugger << keyval->describedclassvalue(":debug");
//   if (debugger.nonnull()) {
//     debugger->set_exec(argv[0]);
//     debugger->set_prefix(msg->me());
//   }

  keyval = 0;

  msg->sync();

  MessageGrp::set_default_messagegrp(msg);

  Ref<MemoryGrp> mem = MemoryGrp::initial_memorygrp(argc, argv);
  if (mem.nonnull())
    MemoryGrp::set_default_memorygrp(mem);
  else
    mem = MemoryGrp::get_default_memorygrp();

  do_memgrp_tests(msg, mem);
  do_memgrpregion_tests(msg, mem);
  
  MemoryGrp::set_default_memorygrp(0);
  MessageGrp::set_default_messagegrp(0);

  return 0;
}

void
do_memgrp_tests(const Ref<MessageGrp>&msg,
             const Ref<MemoryGrp>&mem)
{
  do_simple_tests(msg, mem);

  do_double_tests(msg, mem);
  do_double2_tests(msg, mem);

  do_int_tests(msg, mem);
}

void
do_simple_tests(const Ref<MessageGrp>&msg,
                const Ref<MemoryGrp>&mem)
{
  mem->set_localsize(8);

  cout << scprintf("Using memory group \"%s\".\n", mem->class_name());

  mem->sync();
  mem->set_localsize(0);
}

void
do_int_tests(const Ref<MessageGrp>&msg,
             const Ref<MemoryGrp>&mem)
{
  const int intbufsize = 10;

  mem->set_localsize(intbufsize*sizeof(int));

  cout << scprintf("Using memory group \"%s\".\n", mem->class_name());

  //sleep(1);
  cout.flush();
  cout << scprintf("111111111111111111111111111111111\n");
  cout.flush();
  //sleep(1);

  mem->sync();

  //sleep(1);
  cout.flush();
  cout << scprintf("222222222222222222222222222222222\n");
  cout.flush();
  //sleep(1);

  //mem->deactivate();
  //sleep(1);
  cout.flush();
  cout << scprintf("333333333333333333333333333333333\n");
  cout.flush();
  //sleep(1);
  //mem = 0;
  //return;

  PRINTF(("creating MemoryGrpBuf\n"));

  MemoryGrpBuf<int> buf(mem);

  PRINTF(("%d: obtaining writelock\n", mem->me()));
  int *data = buf.writeonly(0, intbufsize);
  PRINTF(("%d: releasing writelock\n", mem->me()));
  int i;
  for (i=0; i<intbufsize; i++) data[i] = 0;
  buf.release();

  PRINTF(("%d: obtaining readlock\n", mem->me()));
  const int *cdata = buf.readonly(0, intbufsize);
  PRINTF(("%d: releasing readlock\n", mem->me()));
  buf.release();

//   if (mem->me() == 0) {
//       cdata = buf.readonly(0, intbufsize);
//       for (i=0; i<intbufsize; i++) {
//           cout << scprintf("data[%3d] = %4d\n", i, cdata[i]);
//         }
//       buf.release();
//     }

  PRINTF(("%d: syncing\n", mem->me()));
  mem->sync();
  PRINTF(("---------------------------------------------------------\n"));
  mem->sync();
  PRINTF(("%d: done syncing\n", mem->me()));

  if (mem->me() == 0 || mem->me() == 1) {
      int start[2];
      int length[2];
      start[0] = 0;
      length[0] = intbufsize/2;
      start[1] = length[0];
      length[1] = intbufsize - length[0];
      data = buf.readwrite(start[mem->me()], length[mem->me()]);
      PRINTF(("%d: adding %d to [%d, %d)\n",
              mem->me(), 10 * (mem->me()+1),
              start[mem->me()], start[mem->me()]+length[mem->me()]));
      for (i=0; i<buf.length(); i++) {
          data[i] += 10 * (mem->me() + 1);
        }
      buf.release();
      mem->sync();
      PRINTF(("------------------------------------------------------\n"));
      mem->sync();
      data = buf.readwrite(start[1-mem->me()], length[1-mem->me()]);
      PRINTF(("%d: adding %d to [%d, %d)\n",
              mem->me(), 100 * (mem->me()+1),
              start[1-mem->me()], start[1-mem->me()]+length[1-mem->me()]));
      for (i=0; i<buf.length(); i++) {
          data[i] += 100 * (mem->me() + 1);
        }
      buf.release();
      mem->sync();
      PRINTF(("------------------------------------------------------\n"));
      mem->sync();
    }
  else {
      mem->sync();
      PRINTF(("------------------------------------------------------\n"));
      mem->sync();
      mem->sync();
      PRINTF(("------------------------------------------------------\n"));
      mem->sync();
   }

  if (mem->me() == 0) {
      cdata = buf.readonly(0, intbufsize);
      for (i=0; i<intbufsize; i++) {
          PRINTF(("data[%3d] = %4d\n", i, cdata[i]));
        }
      buf.release();
    }

  PRINTF(("%d: syncing\n", mem->me()));
  mem->sync();
  PRINTF(("---------------------------------------------------------\n"));
  mem->sync();

  PRINTF(("%d: exiting\n", mem->me()));
  PRINTF(("%d: syncing\n", mem->me()));
  mem->sync();
  PRINTF(("==========================================================\n"));
  mem->sync();

  mem->set_localsize(0);
}

void
do_double_tests(const Ref<MessageGrp>&msg,
                const Ref<MemoryGrp>&mem)
{
  PRINTF(("double tests entered\n"));

  int i,j;

  const int doublebufsize = 4;

  mem->set_localsize(doublebufsize*sizeof(double));

  cout << scprintf("Using memory group \"%s\".\n", mem->class_name());

  mem->sync();

  MemoryGrpBuf<double> dbuf(mem);

  PRINTF(("%d: double tests mem = 0x%x\n", mem->me(), mem.pointer()));

  double factor = 1.0;
  for (i=0; i<mem->me(); i++) factor *= 10.0;

  PRINTF(("%d: setting 5th digit\n", mem->me()));

  double *data = dbuf.writeonly_on_node(0, doublebufsize);
  for (i=0; i<doublebufsize; i++) {
      data[i] = 10000.0 * (mem->me()+1) + 100000.0 * i;
    }
  dbuf.release();

  PRINTF(("%d: 5th digit set\n", mem->me()));

  double contrib[doublebufsize];
  for (i=0; i<doublebufsize; i++) {
      contrib[i] = (mem->me()+1) * factor;
    }

  mem->sync();
  PRINTF(("---------------------------------------------------------\n"));
  mem->sync();

  PRINTF(("%d: starting sum reduction\n", mem->me()));

  for (i=0; i<mem->n(); i++) {
      mem->sum_reduction_on_node(contrib, mem->me(),
                                 doublebufsize-mem->me(),
                                 i);
    }

  PRINTF(("%d: done with sum reduction\n", mem->me()));

  mem->sync();
  PRINTF(("---------------------------------------------------------\n"));
  mem->sync();

  for (i=0; i<mem->n(); i++) {
      mem->sync();
      if (i==mem->me()) {
          const double *cdata = dbuf.readonly_on_node(0, doublebufsize);
          for (j=0; j<doublebufsize; j++) {
              PRINTF(("%2d: %12.1f\n", i, cdata[j]));
            }
          dbuf.release();
        }
    }

  mem->sync();
  PRINTF(("---------------------------------------------------------\n"));
  mem->sync();

  for (i=0; i<mem->n(); i++) {
      mem->sync();
      if (i==mem->me()) {
          const double *cdata = dbuf.readonly(0, doublebufsize*mem->n());
          for (j=0; j<doublebufsize*mem->n(); j++) {
              PRINTF(("%2d: data[%2d] = %12.1f\n", i, j, cdata[j]));
            }
          dbuf.release();
        }
    }

  mem->sync();
  PRINTF(("==========================================================\n"));
  mem->sync();

  mem->set_localsize(0);
}

void
do_double2_tests(const Ref<MessageGrp>&msg,
                 const Ref<MemoryGrp>&mem)
{
  PRINTF(("double2 tests entered\n"));

  int i,j;

  const int doublebufsize = 4;
  mem->set_localsize(doublebufsize*sizeof(double));
  cout << scprintf("Using memory group \"%s\".\n", mem->class_name());

  mem->sync();

  MemoryGrpBuf<double> dbuf(mem);

  double *data = dbuf.writeonly_on_node(0, doublebufsize);
  for (i=0; i<doublebufsize; i++) {
      data[i] = 0.0;
    }
  dbuf.release();

  int target = (mem->me() + 1)%mem->n();

  double contrib[doublebufsize];
  for (i=0; i<doublebufsize; i++) {
      contrib[i] = mem->me()/200.0;
    }

  mem->sync();
  PRINTF(("---------------------------------------------------------\n"));
  mem->sync();

  PRINTF(("%d: starting sum reduction\n", mem->me()));

  for (i=0; i<200; i++) {
      mem->sum_reduction_on_node(contrib, 0,
                                 doublebufsize,
                                 target);
    }

  PRINTF(("%d: done with sum reduction\n", mem->me()));

  mem->sync();
  PRINTF(("---------------------------------------------------------\n"));
  mem->sync();

  for (i=0; i<mem->n(); i++) {
      mem->sync();
      if (i==mem->me()) {
          const double *cdata = dbuf.readonly(0, doublebufsize*mem->n());
          for (j=0; j<doublebufsize*mem->n(); j++) {
              PRINTF(("%2d: data[%2d] = %12.1f\n", i, j, cdata[j]));
            }
          dbuf.release();
        }
    }

  mem->sync();
  PRINTF(("==========================================================\n"));
  mem->sync();

  mem->set_localsize(0);
}

void
do_memgrpregion_tests(const Ref<MessageGrp>&msg,
                      const Ref<MemoryGrp>&mem)
{
  PRINTF(("MemoryGrpRegion tests entered\n"));
  cout << scprintf("Using host memory group \"%s\".\n", mem->class_name());  
  // grab enough memory to fit in 2 regions
  mem->set_localsize(400);

  Ref<MemoryGrpRegion> rg1 = new MemoryGrpRegion(mem,0,200);
  do_memgrp_tests(msg, rg1);
  Ref<MemoryGrpRegion> rg2 = new MemoryGrpRegion(mem,200,200);
  do_memgrp_tests(msg, rg2);

  mem->sync();
  PRINTF(("==========================================================\n"));
  mem->sync();
  mem->set_localsize(0);
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:

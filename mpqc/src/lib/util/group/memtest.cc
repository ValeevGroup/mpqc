
#include <stdio.h>
#include <math.h>
#include <util/group/message.h>
#include <util/group/mstate.h>
#include <util/group/hcube.h>
#include <util/group/memshm.h>
#ifdef HAVE_NX
#  include <util/group/memipgon.h>
#endif

// Force linkages:
//#ifndef __PIC__
#ifdef HAVE_SYSV_IPC
#   include <util/group/messshm.h>
    const ClassDesc &fl0 = ShmMessageGrp::class_desc_;
#endif
#ifdef HAVE_PVM
#   include <util/group/messpvm.h>
    const ClassDesc &fl2 = PVMMessageGrp::class_desc_;
#endif
#ifdef HAVE_MPI
#   include <util/group/messmpi.h>
    const ClassDesc &fl2 = MPIMessageGrp::class_desc_;
#endif
//#endif

// this is needed for debugging
#ifdef HAVE_NX
extern "C" {
#include <nx.h>
}
#endif // HAVE_NX

#ifdef HAVE_HRECV
#  define DISABLE do { masktrap(1); fflush(stdout); } while(0)
#  define ENABLE do { fflush(stdout); masktrap(0); } while(0)
   extern "C" {
       long masktrap(long state);
     }
#else
#  define DISABLE
#  define ENABLE
#endif

#define PRINTF(args) do { DISABLE; \
                          printf args; \
                          fflush(stdout); \
                          ENABLE; \
                         } while(0)

void do_simple_tests(const RefMessageGrp&);
void do_int_tests(const RefMessageGrp&);
void do_double_tests(const RefMessageGrp&);
void do_double2_tests(const RefMessageGrp&);

#ifdef HAVE_HRECV
#  define MemoryGrp_CTOR(msg) new IParagonMemoryGrp(msg)
#else
#  define MemoryGrp_CTOR(msg) MemoryGrp::create_memorygrp()
#endif

int
main(int argc, char**argv)
{
  RefMessageGrp msg = MessageGrp::initial_messagegrp(argc, argv);

  if (msg.null()) {
      const char* input = SRCDIR "/memtest.in";
      const char* keyword = "message";

      if (argc >= 2) input = argv[1];
      if (argc >= 3) keyword = argv[2];

      RefKeyVal keyval = new ParsedKeyVal(input);

      msg = keyval->describedclassvalue(keyword);

      if (msg.null()) {
          fprintf(stderr,"Couldn't initialize MessageGrp\n");
          abort();
        }
    }

  msg->sync();

  MessageGrp::set_default_messagegrp(msg);

  do_simple_tests(msg);

  do_double_tests(msg);
  do_double2_tests(msg);

  do_int_tests(msg);

  return 0;
}

void
do_simple_tests(const RefMessageGrp&msg)
{
  RefMemoryGrp mem = MemoryGrp_CTOR(msg);

  mem->set_localsize(8);

  printf("Using memory group \"%s\".\n", mem->class_name());

  mem->sync();
  mem = 0;
}

void
do_int_tests(const RefMessageGrp&msg)
{
  const int intbufsize = 10;
  RefMemoryGrp mem = MemoryGrp_CTOR(msg);

  mem->set_localsize(intbufsize*sizeof(int));

  printf("Using memory group \"%s\".\n", mem->class_name());

  //sleep(1);
  fflush(stdout);
  printf("111111111111111111111111111111111\n");
  fflush(stdout);
  //sleep(1);

  mem->sync();

  //sleep(1);
  fflush(stdout);
  printf("222222222222222222222222222222222\n");
  fflush(stdout);
  //sleep(1);

  //mem->deactivate();
  //sleep(1);
  fflush(stdout);
  printf("333333333333333333333333333333333\n");
  fflush(stdout);
  //sleep(1);
  //mem = 0;
  //return;

  int uselocks = 0;

  mem->lock(uselocks);

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
//           printf("data[%3d] = %4d\n", i, cdata[i]);
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
      if (!uselocks) {
          mem->sync();
          PRINTF(("------------------------------------------------------\n"));
          mem->sync();
        }
      data = buf.readwrite(start[1-mem->me()], length[1-mem->me()]);
      PRINTF(("%d: adding %d to [%d, %d)\n",
              mem->me(), 100 * (mem->me()+1),
              start[1-mem->me()], start[1-mem->me()]+length[1-mem->me()]));
      for (i=0; i<buf.length(); i++) {
          data[i] += 100 * (mem->me() + 1);
        }
      buf.release();
      if (!uselocks) {
          mem->sync();
          PRINTF(("------------------------------------------------------\n"));
          mem->sync();
        }
    }
  else {
      if (!uselocks) {
          mem->sync();
          PRINTF(("------------------------------------------------------\n"));
          mem->sync();
          mem->sync();
          PRINTF(("------------------------------------------------------\n"));
          mem->sync();
        }
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

  //mem->deactivate();
  mem = 0;
}

void
do_double_tests(const RefMessageGrp&msg)
{
  PRINTF(("double tests entered\n"));

  int i,j;

  RefMemoryGrp mem;

  const int doublebufsize = 4;
  mem = MemoryGrp_CTOR(msg);

  mem->set_localsize(doublebufsize*sizeof(double));

  printf("Using memory group \"%s\".\n", mem->class_name());

  mem->sync();
  mem->lock(0);

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

  mem->lock(1);
  for (i=0; i<mem->n(); i++) {
      mem->sum_reduction_on_node(contrib, mem->me(),
                                 doublebufsize-mem->me(),
                                 i);
    }
  mem->lock(0);

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

  //mem->deactivate();
}

void
do_double2_tests(const RefMessageGrp&msg)
{
  PRINTF(("double2 tests entered\n"));

  int i,j;

  RefMemoryGrp mem;

  const int doublebufsize = 4;
  mem = MemoryGrp_CTOR(msg);
  mem->set_localsize(doublebufsize*sizeof(double));
  printf("Using memory group \"%s\".\n", mem->class_name());

  mem->sync();
  mem->lock(0);

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

  mem->lock(1);
  for (i=0; i<200; i++) {
      mem->sum_reduction_on_node(contrib, 0,
                                 doublebufsize,
                                 target);
    }
  mem->lock(0);

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

  //mem->deactivate();
}

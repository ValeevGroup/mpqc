
#include <stdio.h>
#include <math.h>
#include <util/misc/bug.h>
#include <util/group/messshm.h>
#include <math/scmat/dist.h>

ClassDesc* f0 = &ShmMessageGrp::class_desc_;
#ifdef HAVE_MPI
#include <util/group/messmpi.h>
ClassDesc* f1 = &MPIMessageGrp::class_desc_;
#endif

void matrixtest(RefSCMatrixKit kit, RefKeyVal keyval,
                RefSCDimension d1,RefSCDimension d2,RefSCDimension d3);

main(int argc, char** argv)
{
  RefKeyVal keyval = new ParsedKeyVal(SRCDIR "/matrixtest.in");

  RefMessageGrp msg = MessageGrp::initial_messagegrp(argc, argv);

  if (msg.null()) {
      msg = keyval->describedclassvalue("messagegrp");

      if (msg.null()) {
          fprintf(stderr,"Couldn't initialize MessageGrp\n");
          abort();
        }
    }

  MessageGrp::set_default_messagegrp(msg);

  RefDebugger d = keyval->describedclassvalue("debugger");
  if (d.nonnull()) {
      d->set_prefix(msg->me());
      d->set_exec(argv[0]);
    }

  // test the blocklist send and receive
  if (msg->n() > 1) {
      RefSCMatrixBlockList l = new SCMatrixBlockList;
      l->append(new SCMatrixRectBlock(0,5,0,2));
      l->append(new SCMatrixRectBlock(0,3,0,11));
      l->append(new SCMatrixLTriBlock(7,13));
      if (msg->me() == 0) {
          StateSend out(msg);
          out.target(1);
          out.copy_references();
          l.save_state(out);
          out.flush();
        }
      else if (msg->me() == 1) {
          StateRecv in(msg);
          in.source(0);
          in.copy_references();
          l.restore_state(in);
        }
    }

  RefSCMatrixKit kit = new DistSCMatrixKit;
  RefSCDimension d1(keyval->describedclassvalue("d1"));
  RefSCDimension d2(keyval->describedclassvalue("d2"));
  RefSCDimension d3(keyval->describedclassvalue("d3"));

  int nblocks = (int)sqrt(msg->n());

  // replace dimensions with dimensions that have subblocks
  d1 = new SCDimension(d1.n(), nblocks);
  d2 = new SCDimension(d2.n(), nblocks);
  d3 = new SCDimension(d3.n(), nblocks);

  matrixtest(kit,keyval,d1,d2,d3);

  return 0;
}

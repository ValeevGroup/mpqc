
#include <util/misc/formio.h>
#include <util/group/messshm.h>
#include <math/scmat/repl.h>

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
          cerr << indent << "Couldn't initialize MessageGrp\n";
          abort();
        }
    }

  MessageGrp::set_default_messagegrp(msg);

  RefSCMatrixKit kit = new ReplSCMatrixKit;
  RefSCDimension d1(keyval->describedclassvalue("d1"));
  RefSCDimension d2(keyval->describedclassvalue("d2"));
  RefSCDimension d3(keyval->describedclassvalue("d3"));

  matrixtest(kit,keyval,d1,d2,d3);

  return 0;
}

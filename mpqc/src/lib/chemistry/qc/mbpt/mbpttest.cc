
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string.h>

#include <sys/stat.h>
#include <unistd.h>
#include <new.h>

#include <util/keyval/keyval.h>
#include <util/group/message.h>
#include <util/group/pregtime.h>
#include <util/misc/bug.h>
#include <util/misc/formio.h>

#include <math/scmat/repl.h>
#include <math/scmat/dist.h>

#include <math/optimize/qnewton.h>
#include <math/optimize/gdiis.h>
#include <math/optimize/efc.h>
#include <math/optimize/update.h>

#include <chemistry/molecule/coor.h>
#include <chemistry/molecule/energy.h>

#include <chemistry/qc/scf/clscf.h>
#include <chemistry/qc/scf/hsosscf.h>

#include <chemistry/qc/mbpt/mbpt.h>

// Force linkages:
#ifndef __PIC__
const ClassDesc &fl0a = CLSCF::class_desc_;
const ClassDesc &fl0b = HSOSSCF::class_desc_;

const ClassDesc &fl0e = MBPT2::class_desc_;

const ClassDesc &fl1a = RedundMolecularCoor::class_desc_;
const ClassDesc &fl1b = CartMolecularCoor::class_desc_;
const ClassDesc &fl1c = SymmMolecularCoor::class_desc_;

const ClassDesc &fl2 = QNewtonOpt::class_desc_;
const ClassDesc &fl3 = GDIISOpt::class_desc_;
const ClassDesc &fl4 = EFCOpt::class_desc_;
const ClassDesc &fl5 = BFGSUpdate::class_desc_;

const ClassDesc &fl6 = ReplSCMatrixKit::class_desc_;
const ClassDesc &fl7 = DistSCMatrixKit::class_desc_;

# ifdef HAVE_SYSV_IPC
#   include <util/group/messshm.h>
    const ClassDesc &fl8 = ShmMessageGrp::class_desc_;
# endif
const ClassDesc &fl9 = ProcMessageGrp::class_desc_;
# ifdef HAVE_NX_H
#  include <util/group/messpgon.h>
    const ClassDesc &fl10 = ParagonMessageGrp::class_desc_;
# endif
#endif

RefRegionTimer tim;
RefMessageGrp grp;

static RefMessageGrp
init_mp(const RefKeyVal& keyval)
{
  // if we are on a paragon then use a ParagonMessageGrp
  // otherwise read the message group from the input file
#ifdef HAVE_NX_H
  grp = new ParagonMessageGrp;
#else
  grp = keyval->describedclassvalue("message");
#endif

  if (grp.nonnull()) MessageGrp::set_default_messagegrp(grp);
  else grp = MessageGrp::get_default_messagegrp();

  RefDebugger debugger = keyval->describedclassvalue(":debug");
  // Let the debugger know the name of the executable and the node
  if (debugger.nonnull()) {
    debugger->set_exec("mbpttest");
    debugger->set_prefix(grp->me());
    debugger->debug("curt is a hog");
  }
  
  tim = new ParallelRegionTimer(grp,"mbpttest",1,0);
  RegionTimer::set_default_regiontimer(tim);

  SCFormIO::set_printnode(0);
  SCFormIO::set_messagegrp(grp);
  //SCFormIO::set_debug(1);

  SCFormIO::setindent(cout, 2);
  SCFormIO::setindent(cerr, 2);
  
  return grp;
}

main(int argc, char**argv)
{
  char *input =      (argc > 1)? argv[1] : SRCDIR "/mbpttest.in";
  char *keyword =    (argc > 2)? argv[2] : "mole";
  char *optkeyword = (argc > 3)? argv[3] : "opt";

  // open keyval input
  RefKeyVal rpkv(new ParsedKeyVal(input));

  init_mp(rpkv);

  tim->enter("input");
  
  if (rpkv->exists("matrixkit"))
    SCMatrixKit::set_default_matrixkit(rpkv->describedclassvalue("matrixkit"));
  
  struct stat sb;
  RefMolecularEnergy mole;
  RefOptimize opt;

  if (stat("mbpttest.ckpt",&sb)==0 && sb.st_size) {
    StateInBinXDR si("mbpttest.ckpt");
    opt.restore_state(si);
    mole = opt->function();
  } else {
    mole = rpkv->describedclassvalue(keyword);
    opt = rpkv->describedclassvalue(optkeyword);
    if (opt.nonnull()) {
      opt->set_checkpoint();
      opt->set_checkpoint_file("mbpttest.ckpt");
    }
  }

  tim->exit("input");

  if (mole.nonnull()) {
    if (mole->gradient_implemented()) {
      if (opt.nonnull()) {
        opt->optimize();
      } else {
        mole->gradient().print("gradient");
      }
    } else if (mole->value_implemented()) {
      cout << node0 << indent
           << scprintf("value of mole is %20.15f\n\n", mole->energy());
    }

  mole->print(cout);
  }

  StateOutBinXDR so("mbpttest.wfn");
  mole.save_state(so);
  
  tim->print(cout);

  return 0;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// eval: (c-set-style "CLJ")

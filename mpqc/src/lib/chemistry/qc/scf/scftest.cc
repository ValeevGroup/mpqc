
#include <string.h>

#include <sys/stat.h>
#include <unistd.h>
#include <new.h>

#include <util/keyval/keyval.h>

#include <chemistry/qc/scf/clscf.h>
#include <chemistry/qc/scf/hsosscf.h>
#include <chemistry/qc/scf/ossscf.h>
#include <chemistry/qc/scf/tcscf.h>
#include <chemistry/qc/scf/xscf.h>
#include <chemistry/qc/scf/mcscf.h>
#include <chemistry/molecule/coor.h>
#include <math/optimize/qnewton.h>
#include <math/optimize/gdiis.h>
#include <math/optimize/efc.h>
#include <math/optimize/update.h>
#include <util/group/message.h>

// Force linkages:
#ifndef __PIC__
const ClassDesc &fl0  = CLSCF::class_desc_;
const ClassDesc &fl0a = HSOSSCF::class_desc_;
const ClassDesc &fl0b = OSSSCF::class_desc_;
const ClassDesc &fl0c = TCSCF::class_desc_;
const ClassDesc &fl0d = XSCF::class_desc_;
const ClassDesc &fl0e = MCSCF::class_desc_;
const ClassDesc &fl1a = RedundMolecularCoor::class_desc_;
const ClassDesc &fl1b = CartMolecularCoor::class_desc_;
const ClassDesc &fl1c = SymmMolecularCoor::class_desc_;
const ClassDesc &fl2 = QNewtonOpt::class_desc_;
const ClassDesc &fl3 = GDIISOpt::class_desc_;
const ClassDesc &fl4 = EFCOpt::class_desc_;
const ClassDesc &fl5 = BFGSUpdate::class_desc_;
# ifndef PARAGON
#   include <util/group/messshm.h>
    const ClassDesc &fl8 = ShmMessageGrp::class_desc_;
    const ClassDesc &fl9 = ProcMessageGrp::class_desc_;
# endif
# ifdef PARAGON
#  include <util/group/messpgon.h>
    const ClassDesc &fl10 = ParagonMessageGrp::class_desc_;
# endif
#endif

static RefMessageGrp
init_mp(const char *inputfile)
{
  RefMessageGrp grp;

  // if we are on a paragon then use a ParagonMessageGrp
  // otherwise read the message group from the input file
#if defined(PARAGON)
  grp = new ParagonMessageGrp;
#else
  RefKeyVal keyval = new ParsedKeyVal(inputfile);
  grp = keyval->describedclassvalue("message");
  keyval = 0;
#endif

  if (grp.nonnull()) MessageGrp::set_default_messagegrp(grp);
  else grp = MessageGrp::get_default_messagegrp();

  return grp;
}

main(int argc, char**argv)
{
  // the output stream is standard out
  ostream& o = cout;

  char *input =      (argc > 1)? argv[1] : SRCDIR "/mpqc.in";
  char *keyword =    (argc > 2)? argv[2] : "mole";
  char *optkeyword = (argc > 3)? argv[3] : "opt";

  init_mp(input);

  struct stat sb;
  RefMolecularEnergy mole;
  RefOptimize opt;

  if (stat("mpqctest.ckpt",&sb)==0 && sb.st_size) {
    StateInBinXDR si("mpqctest.ckpt");
    opt.restore_state(si);
    mole = opt->nlp();
  } else {
    // open keyval input
    RefKeyVal rpkv(new ParsedKeyVal(input));

    mole = rpkv->describedclassvalue(keyword);
    opt = rpkv->describedclassvalue(optkeyword);
    // opt->set_checkpoint();
    // opt->set_checkpoint_file("mpqctest.ckpt");
  }

  if (mole.nonnull()) {
    if (mole->gradient_implemented()) {
      if (opt.nonnull()) {
        //opt->print(o);
        opt->optimize();
      } else {
        mole->gradient().print("gradient");
      }
    } else if (mole->value_implemented()) {
      o << " value of mole is ";
      o << mole->energy() << endl;
      printf("%20.15f\n",mole->energy());
    }

    mole->print(o);
  }

  return 0;
}


#include <string.h>

#include <sys/stat.h>
#include <unistd.h>

#include <util/keyval/keyval.h>

#include "g92.h"

#include <chemistry/molecule/coor.h>
#include <math/optimize/qnewton.h>
#include <math/optimize/gdiis.h>
#include <math/optimize/efc.h>
#include <math/optimize/update.h>

// Force linkages:
#ifndef __PIC__
const ClassDesc &fl0  = Gaussian92SCF::class_desc_;
const ClassDesc &fl0a = Gaussian92UHF::class_desc_;

const ClassDesc &fl1 = IntMolecularCoor::class_desc_;
const ClassDesc &fl2 = QNewtonOpt::class_desc_;
const ClassDesc &fl3 = GDIISOpt::class_desc_;
const ClassDesc &fl4 = EFCOpt::class_desc_;
const ClassDesc &fl5 = BFGSUpdate::class_desc_;
#endif

main(int argc, char**argv)
{
  // the output stream is standard out
  SCostream& o = SCostream::cout;

  char *input =      (argc > 1)? argv[1] : SRCDIR "/mpqc.in";
  char *keyword =    (argc > 2)? argv[2] : "mole";
  char *optkeyword = (argc > 3)? argv[3] : "opt";

  struct stat sb;
  RefMolecularEnergy mole;
  RefOptimize opt;

  if (stat("g92test.ckpt",&sb)==0 && sb.st_size) {
    //StateInText si("g92test.ckpt");
    StateInBinXDR si("g92test.ckpt");
    opt.restore_state(si);
    mole = opt->nlp();
  } else {
    // open keyval input
    RefKeyVal rpkv(new ParsedKeyVal(input));

    mole = rpkv->describedclassvalue(keyword);
    opt = rpkv->describedclassvalue(optkeyword);
    opt->set_checkpoint();
    opt->set_checkpoint_file("g92test.ckpt");
  }

#if 0
  if (mole->gradient_implemented()) {
    if (opt.nonnull()) {
      //opt->print(o);
      opt->optimize();
    } else {
      o << "opt is null\n";
    }
  }
#else
  mole->do_hessian(1);
  mole->hessian();
  Gaussian92::castdown(mole)->normal_modes()->print("normal modes");
  Gaussian92::castdown(mole)->frequencies()->print("frequencies");
#endif

  return 0;
}

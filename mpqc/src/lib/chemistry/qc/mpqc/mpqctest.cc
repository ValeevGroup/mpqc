
#include <string.h>

#include <util/keyval/keyval.h>
#include <new.h>
#include "mpqc.h"
#include <chemistry/molecule/coor.h>
#include <math/optimize/qnewton.h>
#include <math/optimize/gdiis.h>
#include <math/optimize/efc.h>

// Force linkages:
#ifndef __PIC__
const ClassDesc &fl0 = MPSCF::class_desc_;
const ClassDesc &fl1 = IntMolecularCoor::class_desc_;
const ClassDesc &fl2 = QNewtonOpt::class_desc_;
const ClassDesc &fl3 = GDIISOpt::class_desc_;
const ClassDesc &fl4 = EFCOpt::class_desc_;
#endif

void die()
{
  fprintf(stderr,"die\n");
  abort();
}

main(int argc, char**argv)
{
  set_new_handler(die);
  
/* prepare for potential debugging */
#if (defined(SGI) || defined(SUN4)) && defined(USE_DEBUG)
  debug_init(argv[0]);
#endif

#if ((defined(SGI) || defined(SUN4)) && \
     (!defined(SABER))) && defined(USE_DEBUG)
  malloc_debug_on_error();
#endif

  // the output stream is standard out
  SCostream& o = SCostream::cout;

  char *input = (argv[1]) ? argv[1] : strdup(SRCDIR "/mpqc.in");

  // open keyval input
  RefKeyVal rpkv(new ParsedKeyVal(input));

  RefMolecularEnergy mole = rpkv->describedclassvalue("mole");
     
  if (mole.nonnull()) {
    mole->print(o);

    o << "energy = " << mole->energy() << endl;
    o << "gradient:\n";
    o++; mole->gradient().print(o); o--;
  } else {
    o << "mole is null\n";
  }

  RefOptimize opt = rpkv->describedclassvalue("opt");

  if (opt.nonnull()) {
    //opt->print(o);
    opt->optimize();
  } else {
    o << "opt is null\n";
  }

  return 0;
}

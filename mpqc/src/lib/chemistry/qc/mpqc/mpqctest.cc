
#include <string.h>

#include <math/optimize/opt.h>
#include <util/keyval/keyval.h>
#include <new.h>
#include "mpqc.h"
#include <chemistry/molecule/coor.h>
#include <math/optimize/opt.h>
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

#if ((defined(SGI) || defined(SUN4)) && (!defined(SABER))) && defined(USE_DEBUG)
  malloc_debug_on_error();
#endif

  // the output stream is standard out
  SCostream& o = SCostream::cout;

  char *input = (argv[1]) ? argv[1] : strdup(SRCDIR "/mpqc.in");

  // open keyval input
  RefKeyVal rpkv(new ParsedKeyVal(input));
  int nmole = rpkv->count("mole");
  int nopt = rpkv->count("opt");

  for (int i=0; i < nmole; i++) {
      RefMolecularEnergy mole = rpkv->describedclassvalue("mole",i);
     
      if (mole.nonnull()) {
          mole->print(o);

          o << "energy = " << mole->energy() << endl;
          o << "gradient:\n";
          o++; mole->gradient().print(o); o--;
        }
      else {
          o << "mole[" << i << "] is null\n";
        }
    }

  for (i=0; i < nopt; i++) {
      RefOptimize opt = rpkv->describedclassvalue("opt",i);

      if (opt.nonnull()) {
          //opt->print(o);

          opt->optimize();
        }
      else {
          o << "opt[" << i << "] is null\n";
        }
    }
  return 0;
}

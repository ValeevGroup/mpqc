
#include <math/optimize/opt.h>
#include <util/keyval/keyval.h>
#include <new.h>
#include "mpqc.h"

// Force MPSCF linkage:
const ClassDesc &tmp = MPSCF::class_desc_;

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

  // open keyval input
  RefKeyVal rpkv(new ParsedKeyVal(SRCDIR "/mpqc.in"));

  for (int i=0; rpkv->exists("mole",i); i++) {
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

  for (i=0; rpkv->exists("opt",i); i++) {
      RefOptimize opt = rpkv->describedclassvalue("opt",i);

      if (opt.nonnull()) {
          //opt->print(o);

          opt->optimize();
        }
      else {
          o << "opt[" << i << "] is null\n";
        }
    }
}

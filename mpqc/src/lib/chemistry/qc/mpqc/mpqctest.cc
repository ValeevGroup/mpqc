
#include <string.h>

#include <util/keyval/keyval.h>
#include <new.h>
#include "mpqc.h"
#include <chemistry/molecule/coor.h>
#include <math/optimize/qnewton.h>
#include <math/optimize/gdiis.h>
#include <math/optimize/efc.h>
#include <math/optimize/update.h>
#include <chemistry/molecule/molshape.h>
#include <chemistry/solvent/bem.h>

// Force linkages:
#ifndef __PIC__
const ClassDesc &fl0 = MPSCF::class_desc_;
const ClassDesc &fl1 = IntMolecularCoor::class_desc_;
const ClassDesc &fl2 = QNewtonOpt::class_desc_;
const ClassDesc &fl3 = GDIISOpt::class_desc_;
const ClassDesc &fl4 = EFCOpt::class_desc_;
const ClassDesc &fl5 = ConnollyShape2::class_desc_;
const ClassDesc &fl6 = BEMSolvent::class_desc_;
const ClassDesc &fl7 = BFGSUpdate::class_desc_;
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

  char *input =      (argc > 1)? argv[1] : SRCDIR "/mpqc.in";
  char *keyword =    (argc > 2)? argv[2] : "mole";
  char *optkeyword = (argc > 3)? argv[3] : "opt";

  // open keyval input
  RefKeyVal rpkv(new ParsedKeyVal(input));

  RefMolecularEnergy mole = rpkv->describedclassvalue(keyword);

  if (mole->gradient_implemented()) {
      RefOptimize opt = rpkv->describedclassvalue(optkeyword);

      if (opt.nonnull()) {
          //opt->print(o);
          opt->optimize();
        } else {
            o << "opt is null\n";
          }
    }

  return 0;
}


#include <util/misc/formio.h>
#include <chemistry/qc/psi/psi.h>
#include <math/optimize/opt.h>
#include <util/keyval/keyval.h>
#include <new>

using namespace std;

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
  ostream& o = cout;

  ParsedKeyVal* pkv;
  Ref<KeyVal> rpkv(pkv = new ParsedKeyVal());
  pkv->read( SRCDIR "/psi.in");
  pkv = 0; // should only use rpkv

  int i;
  for (i=0; rpkv->exists("mole",i); i++) {
      Ref<MolecularEnergy> mole = rpkv->describedclassvalue("mole",i);

      if (mole.nonnull()) {
          mole->print(o);

          o << "energy = " << mole->energy() << endl;
          o << "gradient:\n";
          o << incindent; mole->gradient().print(o); o << decindent;
        }
      else {
          o << "mole[" << i << "] is null\n";
        }
    }

  for (i=0; rpkv->exists("opt",i); i++) {
      Ref<Optimize> opt = rpkv->describedclassvalue("opt",i);

      if (opt.nonnull()) {
          //opt->print(o);

          opt->optimize();
        }
      else {
          o << "opt[" << i << "] is null\n";
        }
      opt->function()->print();
    }
}

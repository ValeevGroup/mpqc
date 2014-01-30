
#include <util/misc/formio.h>
#include <chemistry/qc/psi/psiwfn.h>
#include <math/optimize/opt.h>
#include <util/keyval/keyval.h>

#include <chemistry/qc/psi/linkage.h>
#include <math/optimize/linkage.h>

using namespace std;
using namespace sc;

void die()
{
  fprintf(stderr,"die\n");
  abort();
}

int
main(int argc, char**argv)
{

  set_new_handler(die);
  
  // the output stream is standard out
  ostream& o = cout;

  ParsedKeyVal* pkv;
  Ref<KeyVal> rpkv(pkv = new ParsedKeyVal());
  pkv->read( SRCDIR "/psi.in");
  pkv = 0; // should only use rpkv

  int i, do_grad = 1;

  for (i=0; rpkv->exists("mole",i); i++) {
      Ref<MolecularEnergy> mole;
      mole << rpkv->describedclassvalue("mole",i);
      if (do_grad)
	mole->do_gradient(1);
      else
	mole->do_gradient(0);

      if (mole) {
	  mole->print(o);

          o << "energy = " << mole->energy() << endl;
	  if (do_grad) {
	      o << "gradient:\n";
	      o << incindent; mole->gradient().print(o); o << decindent;
	  }
        }
      else {
          o << "mole[" << i << "] is null\n";
        }
    }

  for (i=0; rpkv->exists("opt",i); i++) {
      Ref<Optimize> opt;
      opt << rpkv->describedclassvalue("opt",i);

      if (opt) {
          //opt->print(o);

          opt->optimize();
        }
      else {
          o << "opt[" << i << "] is null\n";
        }
      opt->function()->print();
    }
}

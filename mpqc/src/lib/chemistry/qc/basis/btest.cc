
#include <iostream.h>
#include <util/keyval/keyval.h>
#include <chemistry/molecule/molecule.h>
#include <chemistry/qc/basis/basis.h>

int
main(int, char**)
{
  RefKeyVal keyval = new ParsedKeyVal(SRCDIR "/btest.kv");

  for (int i=0; i<keyval->count("test"); i++) {
      RefGaussianBasisSet gbs = keyval->describedclassvalue("test", i);
      gbs->print();

      fflush(stdout);
      cout.flush();

      StateOutText out("btest.out");
      gbs.save_state(out);
      StateInText in("btest.out");
      gbs.restore_state(in);
      gbs->print();
    }

  return 0;
}

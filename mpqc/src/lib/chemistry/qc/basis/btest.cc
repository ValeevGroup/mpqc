
#include <util/keyval/keyval.h>
#include <chemistry/qc/basis/basis.h>

int
main(int, char *argv[])
{
  char *filename = (argv[1]) ? argv[1] : SRCDIR "/btest.kv";
  
  RefKeyVal keyval = new ParsedKeyVal(filename);

  for (int i=0; i<keyval->count("test"); i++) {
      RefGaussianBasisSet gbs = keyval->describedclassvalue("test", i);

      //gbs->print();

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

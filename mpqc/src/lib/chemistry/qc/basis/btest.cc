
#include <iostream.h>
#include <util/keyval/keyval.h>
#include <chemistry/molecule/molecule.h>
#include <chemistry/qc/basis/basis.h>

main()
{
  RefKeyVal keyval = new ParsedKeyVal(SRCDIR "/btest.in");

  RefMolecule molecule = keyval->describedclassvalue("molecule");

  RefGaussianBasisSet gbs = new GaussianBasisSet(molecule,"sto-3g");

  gbs->print();

  StateOutText out("btest.out");

  gbs.save_state(out);

  StateInText in("btest.out");

  gbs.restore_state(in);

  gbs->print();
}

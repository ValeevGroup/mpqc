
#include "molecule.h"
#include "simple.h"
#include "simpleQCList.h"

main()
{
  ParsedKeyVal kv("moltest.in");

  Molecule mol(PrefixKeyVal("molecule",kv));
  mol.print();

  SimpleCoList simp(PrefixKeyVal("simp",kv));
  simp.print();
}


#include "molecule.h"
#include "symm.h"
#include "symmQCList.h"
#include "simple.h"
#include "simpleQCList.h"

main()
{
  RefKeyVal kv(new ParsedKeyVal("moltest.in"));

  RefDescribedClass val = kv->describedclassvalue("fixed");
  RefSymmCoList fixed = val;
  if (val.nonnull() && fixed.null()) {
      fprintf(stderr,"could not convert type %s to SymmCoList\n",
              val->class_name());
      abort();
    }
  val = 0;
  printf("the fixed list\n");
  fixed->print();

  Molecule mol(PrefixKeyVal("molecule",*kv.pointer()));
  mol.print();

  SimpleCoList simp(PrefixKeyVal("simp",*kv.pointer()));
  simp.print();
}

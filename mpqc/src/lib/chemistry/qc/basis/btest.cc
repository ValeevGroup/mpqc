
#include <iostream.h>
#include <util/keyval/keyval.h>
#include <chemistry/molecule/molecule.h>
#include <chemistry/qc/basis/basis.h>
#include <chemistry/qc/basis/petite.h>
#include <chemistry/qc/integral/integralv2.h>

int
main(int, char**)
{
  RefKeyVal keyval = new ParsedKeyVal(SRCDIR "/btest.kv");

  for (int i=0; i<keyval->count("test"); i++) {
      RefGaussianBasisSet gbs = keyval->describedclassvalue("test", i);

      PetiteList pl(gbs);
      pl.print();
      
      RefSCMatrix aotoso = pl.aotoso();
      aotoso.print("aotoso");

      RefSCMatrix sotoao = aotoso.i();
      sotoao.print("sotoao");

      RefSymmSCMatrix s(gbs->basisdim());
      s.assign(0.0);

      RefSCElementOp op = new GaussianOverlapIntv2(gbs);
      s.element_op(op);
      op=0;

      s.print("overlap");

      RefSymmSCMatrix ss = s.clone();
      ss.assign(0.0);
      ss.accumulate_transform(aotoso.t(),s);
      ss.print("blocked ss");
      
      s.assign(0.0);
      s.accumulate_transform(sotoao.t(),ss);
      s.print("unblocked s");
      
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

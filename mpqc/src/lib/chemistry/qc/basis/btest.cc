
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

      gbs->molecule()->print();
      printf("\n");

      PetiteList pl(gbs);
      pl.print(stdout,0);
      
      RefSCMatrix aotoso(pl.AO_basisdim(), pl.SO_basisdim());
      RefSCElementOp sotrans = new AOSO_Transformation(gbs);
      aotoso.element_op(sotrans);
      sotrans=0;
      aotoso.print("aotoso");

      RefSymmSCMatrix s(pl.AO_basisdim());
      s.assign(0.0);
#if 1
      RefSCElementOp op = new GaussianOverlapIntv2(gbs);
      s.element_op(op);
#else
      RefSCElementOp op = new GaussianKineticIntv2(gbs);
      s.element_op(op);
      op=0;
      op = new GaussianNuclearIntv2(gbs);
      s.element_op(op);
#endif      
      op=0;

      RefSymmSCMatrix ss(pl.SO_basisdim());
      ss.assign(0.0);
      ss.accumulate_transform(aotoso.t(),s);

      ss.print("blocked ss");

      RefSCMatrix sotoao = aotoso.i();
      RefSymmSCMatrix us = s.clone();
      us.assign(0.0);
      us.accumulate_transform(sotoao.t(),ss);
      
      us = us-s;
      us.print("us - s");
      
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


#include <util/keyval/keyval.h>
#include <chemistry/qc/basis/basis.h>
#include <chemistry/qc/intv2/integralv2.h>

int
main(int, char *argv[])
{
  char *filename = (argv[1]) ? argv[1] : SRCDIR "/btest.kv";
  
  RefKeyVal keyval = new ParsedKeyVal(filename);
  
  IntegralV2 intv2;

  for (int i=0; i<keyval->count("test"); i++) {
      RefGaussianBasisSet gbs = keyval->describedclassvalue("test", i);

      RefSymmSCMatrix s(gbs->basisdim());
      RefSCElementOp ov =
        new OneBodyIntOp(new OneBodyIntIter(intv2.overlap_int(gbs)));
      s.assign(0.0);
      s.element_op(ov);
      ov=0;
      s.print("overlap");
      
      RefGaussianBasisSet gbs2 = keyval->describedclassvalue("test2", i);
      RefSCMatrix ssq(gbs->basisdim(),gbs2->basisdim());
      ov = new OneBodyIntOp(new OneBodyIntIter(intv2.overlap_int(gbs,gbs2)));
      ssq.assign(0.0);
      ssq.element_op(ov);
      ssq.print("overlap sq");
      ov=0;
      
      //gbs->print();

      fflush(stdout);
      cout.flush();

      StateOutText out("btest.out");
      gbs.save_state(out);
      StateInText in("btest.out");
      gbs.restore_state(in);
      //gbs->print();
    }

  return 0;
}


#include <util/keyval/keyval.h>
#include <chemistry/qc/basis/basis.h>
#include <chemistry/qc/intv3/intv3.h>

int
main(int, char *argv[])
{
  char *filename = (argv[1]) ? argv[1] : SRCDIR "/btest.kv";
  
  RefKeyVal keyval = new ParsedKeyVal(filename);
  
  RefIntegral intgrl = new IntegralV3;

  for (int i=0; i<keyval->count("test"); i++) {
      RefGaussianBasisSet gbs = keyval->describedclassvalue("test", i);
      intgrl->set_basis(gbs);

      RefSymmSCMatrix s(gbs->basisdim(), gbs->matrixkit());
      RefSCElementOp ov =
        new OneBodyIntOp(new OneBodyIntIter(intgrl->overlap()));
      s.assign(0.0);
      s.element_op(ov);
      ov=0;
      s.print("overlap");
      
      RefGaussianBasisSet gbs2 = keyval->describedclassvalue("test2", i);
      RefSCMatrix ssq(gbs->basisdim(),gbs2->basisdim(),gbs->matrixkit());
      intgrl->set_basis(gbs,gbs2);
      ov = new OneBodyIntOp(new OneBodyIntIter(intgrl->overlap()));
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

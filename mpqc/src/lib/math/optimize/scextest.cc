
#include <util/keyval/keyval.h>
#include <math/scmat/local.h>
#include <math/optimize/diis.h>
#include <math/optimize/scextrap.h>
#include <math/optimize/scextrapmat.h>

// Force linkages:
#ifndef __PIC__
const ClassDesc &fl0 = DIIS::class_desc_;
#endif

int
main()
{
  int i;
  
  RefKeyVal keyval = new ParsedKeyVal( SRCDIR "/scextest.in");

  RefSelfConsistentExtrapolation extrap
      = keyval->describedclassvalue("scextrap");

  RefSCDimension dim = new LocalSCDimension(3, "test_dim");

  RefSymmSCMatrix datamat(dim);
  datamat.assign(0.0);
  datamat->shift_diagonal(2.0);

  RefDiagSCMatrix val(dim);
  RefSCMatrix vec(dim,dim);

  // solve f(x) = x

  i = 0;
  while (i < 100 && !extrap->converged()) {
      datamat.diagonalize(val,vec);
      for (int j=0; j<datamat.dim().n(); j++) {
          double v = val.get_element(j);
          val.set_element(j, sqrt(v));
        }
      RefSymmSCMatrix newdatamat(dim);
      newdatamat.assign(0.0);
      newdatamat.accumulate_transform(vec, val);
      RefSymmSCMatrix errormat = newdatamat - datamat;

      datamat.assign(newdatamat);
      RefSCExtrapData data = new SymmSCMatrixSCExtrapData(datamat);
      RefSCExtrapError error = new SymmSCMatrixSCExtrapError(errormat);

      cout << "Iteration " << i << ":" << endl;

      datamat.print("Datamat:");
      errormat.print("Errormat:");

      extrap->extrapolate(data, error);

      i++;
    }

  return 0;
}

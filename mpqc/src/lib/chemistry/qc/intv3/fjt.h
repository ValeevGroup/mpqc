#ifdef __GNUG__
#pragma interface
#endif

#ifndef _chemistry_qc_intv3_fjt_h
#define _chemistry_qc_intv3_fjt_h

#include <util/ref/ref.h>

class FJT: public VRefCount {
  private:
    double **gtable;

    int maxj;
    double *denomarray;
    double wval_infinity;
    int itable_infinity;

    double *int_fjttable;

    int ngtable() const { return maxj + 7; }
  public:
    FJT(int n);
    ~FJT();
    // Returns J-1 doubles.  The use may read/write these values.
    // They will be overwritten with the next call to values.
    // They will be deleted during the call to ~FJT.
    double *values(int J, double T);
};
REF_dec(FJT);

#endif

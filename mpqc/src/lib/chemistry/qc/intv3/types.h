
#ifndef _chemistry_qc_intv3_types_h
#define _chemistry_qc_intv3_types_h

#include <chemistry/qc/basis/gaussbas.h>

/* Types that are used for integrals, but for which we don't need all
 * of the sgen utilities, are defined here. */

class der_centersv3_t {
  public:
    int n;
    RefGaussianBasisSet cs[4];
    int num[4];
    RefGaussianBasisSet ocs; /* The omitted center's centers_t. */
    int onum;        /* The omitted center's number. */
};

#endif


#ifndef _chemistry_qc_intv2_int_cplus_h
#define _chemistry_qc_intv2_int_cplus_h

#include <chemistry/qc/basis/basis.h>
#include <chemistry/qc/intv2/atoms.h>

class RefKeyVal;
int int_read_basis(const RefKeyVal&, char*, const char*, basis_t&);
int int_read_centers(const RefKeyVal&, centers_t&);

centers_t * int_centers_from_gbs(const RefGaussianBasisSet&);

#endif

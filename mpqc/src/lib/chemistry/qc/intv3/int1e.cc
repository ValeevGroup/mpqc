#ifdef __GNUG__
#pragma implementation
#endif

#include <chemistry/qc/intv3/int1e.h>
#include <chemistry/qc/intv3/utils.h>

Int1eV3::Int1eV3(const RefGaussianBasisSet&b1,
                 const RefGaussianBasisSet&b2,
                 int order)
{
  exponent_weighted = -1;
  scale_shell_result = 0;
  result_scale_factor = 1.0;
  three_center = 0;
  init_order = -1;
  buff = 0;
  cartesianbuffer = 0;

  bs1_ = b1;
  bs2_ = b2;

  int_initialize_offsets1();
  int_initialize_1e(0,order);
}

Int1eV3::~Int1eV3()
{
  int_done_1e();
  int_done_offsets1();
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// eval: (c-set-style "CLJ")

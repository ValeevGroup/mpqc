#ifdef __GNUG__
#pragma implementation
#endif

#include <chemistry/qc/intv3/int1e.h>
#include <chemistry/qc/intv3/utils.h>

extern "C" {
#include <chemistry/qc/intv2/atoms.h>
#include <chemistry/qc/intv2/atomsfree.h>
}

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

  cs1 = int_centers_from_gbs(b1);
  if (b1 == b2) {
      cs2 = cs1;
    }
  else {
      cs2 = int_centers_from_gbs(b2);
    }

  int_initialize_offsets1(cs1,cs2);
  int_initialize_1e(0,order);
}

Int1eV3::~Int1eV3()
{
  int_done_1e();
  int_done_offsets1(cs1,cs2);

  free_centers(cs1);
  free(cs1);
  if (bs1_ != bs2_) {
      free_centers(cs2);
      free(cs2);
    }
}

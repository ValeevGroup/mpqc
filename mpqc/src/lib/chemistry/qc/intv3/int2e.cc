#ifdef __GNUG__
#pragma implementation
#endif

#include <chemistry/qc/intv3/int2e.h>
#include <chemistry/qc/intv3/utils.h>

Int2eV3::Int2eV3(const RefGaussianBasisSet& b1,
                 const RefGaussianBasisSet& b2,
                 const RefGaussianBasisSet& b3,
                 const RefGaussianBasisSet& b4,
                 int order, int storage)
{

  bs1_ = b1;
  bs2_ = b2;
  bs3_ = b3;
  bs4_ = b4;

  int_cs1 = int_centers_from_gbs(b1);

  if (b2 == b1) {
      int_cs2 = int_cs1;
      bs2_ = bs1_;
    }
  else {
      int_cs2 = int_centers_from_gbs(b2);
      bs2_ = b2;
    }

  if (b3 == b1) {
      int_cs3 = int_cs1;
    }
  else if (b3 == b2) {
      int_cs3 = int_cs2;
    }
  else {
      int_cs3 = int_centers_from_gbs(b3);
    }

  if (b4 == b1) {
      int_cs4 = int_cs1;
    }
  else if (b4 == b2) {
      int_cs4 = int_cs2;
    }
  else if (b4 == b3) {
      int_cs4 = int_cs3;
    }
  else {
      int_cs4 = int_centers_from_gbs(b4);
    }

  int_initialize_offsets2(int_cs1,int_cs2,int_cs3,int_cs4);
  int_initialize_erep(storage,order,int_cs1,int_cs2,int_cs3,int_cs4);
}

Int2eV3::~Int2eV3()
{
  int_done_offsets2(int_cs1,int_cs2,int_cs3,int_cs4);
  int_done_erep();
}

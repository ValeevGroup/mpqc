#ifdef __GNUG__
#pragma implementation
#endif

#include <chemistry/qc/intv3/int2e.h>
#include <chemistry/qc/intv3/utils.h>

Int2eV3::Int2eV3(const RefGaussianBasisSet& b1,
                 const RefGaussianBasisSet& b2,
                 const RefGaussianBasisSet& b3,
                 const RefGaussianBasisSet& b4,
                 int order, int storage) :
  store(0), int_Qvec(0), int_Rvec(0)
{
  bs1_ = b1;
  bs2_ = b2;
  bs3_ = b3;
  bs4_ = b4;

  if (bs2_.null()) bs2_ = bs1_;
  if (bs3_.null()) bs3_ = bs2_;
  if (bs4_.null()) bs4_ = bs3_;

  int_initialize_offsets2();
  int_initialize_erep(storage,order,bs1_,bs2_,bs3_,bs4_);
  if (order==0 && ((storage-used_storage_) > 0)) {
    init_bounds();
    init_storage(storage-used_storage_);
  } else if (order==1) {
    init_bounds_1der();
  }
}

Int2eV3::~Int2eV3()
{
  int_done_offsets2();
  int_done_erep();
  if (int_integral_storage) {
    done_storage();
    done_bounds();
    done_bounds_1der();
  }
}

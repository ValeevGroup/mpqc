
#include <stdio.h>

#include <chemistry/qc/intv3/tbintv3.h>

TwoBodyIntV3::TwoBodyIntV3(const RefGaussianBasisSet& b1,
                           const RefGaussianBasisSet& b2,
                           const RefGaussianBasisSet& b3,
                           const RefGaussianBasisSet& b4,
                           int storage):
  TwoBodyInt(b1,b2,b3,b4)
{
  int2ev3_ = new Int2eV3(b1,b2,b3,b4,0,storage);
  buffer_ = int2ev3_->buffer();
}

TwoBodyIntV3::~TwoBodyIntV3()
{
}

void
TwoBodyIntV3::compute_shell(int is, int js, int ks, int ls)
{
  int2ev3_->set_redundant(redundant());
  int2ev3_->erep(is,js,ks,ls);
}

int
TwoBodyIntV3::log2_shell_bound(int is, int js, int ks, int ls)
{
  return int2ev3_->erep_4bound(is,js,ks,ls);
}

//////////////////////////////////////////////////////////////////////////

TwoBodyDerivIntV3::TwoBodyDerivIntV3(const RefGaussianBasisSet& b1,
                                     const RefGaussianBasisSet& b2,
                                     const RefGaussianBasisSet& b3,
                                     const RefGaussianBasisSet& b4,
                                     int storage):
  TwoBodyDerivInt(b1,b2,b3,b4)
{
  int2ev3_ = new Int2eV3(b1,b2,b3,b4,1,storage);
  buffer_ = int2ev3_->buffer();
}

TwoBodyDerivIntV3::~TwoBodyDerivIntV3()
{
}

void
TwoBodyDerivIntV3::compute_shell(int is, int js, int ks, int ls,
                                 DerivCenters&c)
{
  int center;
  der_centersv3_t dercenters;
  int sh[4], sz[4];

  sh[0]=is; sh[1]=js; sh[2]=ks; sh[3]=ls;

  int2ev3_->erep_all1der(sh,sz,&dercenters);

  c.clear();
  for (int i=0; i<dercenters.n; i++) {
      if (dercenters.cs[i] == int2ev3_->cs1()) center = 0;
      else if (dercenters.cs[i] == int2ev3_->cs2()) center = 1;
      else if (dercenters.cs[i] == int2ev3_->cs3()) center = 2;
      else center = 3;
      c.add_center(center,dercenters.num[i]);
    }
  if (dercenters.ocs == int2ev3_->cs1()) center = 0;
  else if (dercenters.ocs == int2ev3_->cs2()) center = 1;
  else if (dercenters.ocs == int2ev3_->cs3()) center = 2;
  else center = 3;
  c.add_omitted(center,dercenters.onum);
}


int
TwoBodyDerivIntV3::log2_shell_bound(int is, int js, int ks, int ls)
{
  return int2ev3_->erep_4bound_1der(is,js,ks,ls);
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// eval: (c-set-style "CLJ")


#include <chemistry/qc/intv3/obintv3.h>

////////////////////////////////////////////////////////////////////////////
// OneBodyIntV3

OneBodyIntV3::OneBodyIntV3(const RefGaussianBasisSet&bs1,
                           const RefGaussianBasisSet&bs2,
                           IntegralFunction ifunc):
  OneBodyInt(bs1,bs2)
{
  int1ev3_ = new Int1eV3(bs1,bs2,0);
  intfunc_ = ifunc;
  buffer_ = int1ev3_->buffer();
}

OneBodyIntV3::~OneBodyIntV3()
{
}

void
OneBodyIntV3::compute_shell(int i, int j)
{
  (int1ev3_->*intfunc_)(i, j);
}

////////////////////////////////////////////////////////////////////////////
// PointChargeIntV3

PointChargeIntV3::PointChargeIntV3(
    const RefGaussianBasisSet&bs1,
    const RefGaussianBasisSet&bs2,
    const RefPointChargeData&dat):
  OneBodyInt(bs1,bs2),
  data_(dat),
  ncharge(0),
  charge(0),
  position(0)
{
  int1ev3_ = new Int1eV3(bs1,bs2,0);
  buffer_ = int1ev3_->buffer();
}

PointChargeIntV3::~PointChargeIntV3()
{
  for (int i=0; i<ncharge; i++) {
      delete position[i];
    }
  if (ncharge) {
      delete[] charge;
      delete[] position;
    }
}

void
PointChargeIntV3::reinitialize()
{
  PointBag_double *charges = data_->charges;

  int nchargenew = charges->length();

  int realloc_charges;
  if (charges->length() != ncharge) {
      ncharge = charges->length();
      realloc_charges = 1;
    }
  else {
      realloc_charges = 0;
    }

  if (ncharge && realloc_charges) {
    position = new double*[ncharge];
    charge = new double[ncharge];
  }
  
  int i = 0;
  for (Pix pix= charges->first(); pix!=0; charges->next(pix)) {
    if (realloc_charges) position[i] = new double[3];
    charge[i] = charges->get(pix);
    for (int j=0; j<3; j++) {
      position[i][j] = charges->point(pix)[j];
    }
    i++;
  }
}

void
PointChargeIntV3::compute_shell(int i,int j)
{
  int1ev3_->point_charge(i,j,ncharge,charge,position);
}

////////////////////////////////////////////////////////////////////////////
// EfieldDotVectorIntV3

EfieldDotVectorIntV3::EfieldDotVectorIntV3(
    const RefGaussianBasisSet&bs1,
    const RefGaussianBasisSet&bs2,
    const RefEfieldDotVectorData&dat) :
  OneBodyInt(bs1,bs2),
  data_(dat)
{
  int1ev3_ = new Int1eV3(bs1,bs2,0);
  buffer_ = int1ev3_->buffer();
}

EfieldDotVectorIntV3::~EfieldDotVectorIntV3()
{
}

void
EfieldDotVectorIntV3::compute_shell(int i,int j)
{
  int nbfi = basis1()->shell(i).nfunction();
  int nbfj = basis2()->shell(j).nfunction();
  int nint = nbfi*nbfj;
  int nint3 = nint*3;
  double *tmp;
  int ii,jj;

  int1ev3_->efield(i,j,data_->position);

  tmp = int1ev3_->buffer();
  for (ii=0; ii<nint; ii++) {
      double tmpval = 0.0;
      for (jj=0; jj<3; jj++) {
          tmpval += *tmp++ * data_->vector[jj];
        }
      buffer_[ii] = tmpval;
    }
}

////////////////////////////////////////////////////////////////////////////
// DipoleIntV3

DipoleIntV3::DipoleIntV3(const RefGaussianBasisSet&bs1,
                         const RefGaussianBasisSet&bs2,
                         const RefDipoleData&dat) :
  OneBodyInt(bs1,bs2),
  data_(dat)
{
  if (data_.null()) {
      data_ = new DipoleData;
    }
}

DipoleIntV3::~DipoleIntV3()
{
}

void
DipoleIntV3::compute_shell(int i,int j)
{
  int1ev3_->dipole(i,j,data_->origin);
}

////////////////////////////////////////////////////////////////////////////
// OneBodyDerivIntV3

OneBodyDerivIntV3::OneBodyDerivIntV3(const RefGaussianBasisSet&bs1,
                                     const RefGaussianBasisSet&bs2,
                                     IntegralFunction ifunc):
  OneBodyDerivInt(bs1,bs2)
{
  int1ev3_ = new Int1eV3(bs1,bs2,1);
  intfunc_ = ifunc;
  buffer_ = int1ev3_->buffer();
}

OneBodyDerivIntV3::~OneBodyDerivIntV3()
{
}

void
OneBodyDerivIntV3::compute_shell(int i, int j, DerivCenters& c)
{
  (int1ev3_->*intfunc_)(i,j,0,basis1()->shell_to_center(i));
  c.clear();
  c.add_center(0,basis1(),i);
  c.add_omitted(1,basis2(),j);
}

void
OneBodyDerivIntV3::compute_shell(int i, int j, int c)
{
  (int1ev3_->*intfunc_)(i,j,0,c);
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// eval: (c-set-style "CLJ")

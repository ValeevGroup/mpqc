
#include <stdio.h>

extern "C" {
#include <chemistry/qc/intv2/int_macros.h>
#include <chemistry/qc/intv2/atoms.h>
#include <chemistry/qc/intv2/atomsfree.h>
#include <chemistry/qc/intv2/comp_1e.gbl>
#include <chemistry/qc/intv2/offsets.gbl>
}

#include <chemistry/qc/intv2/int_cplus.h>
#include <chemistry/qc/intv2/obintv2.h>

////////////////////////////////////////////////////////////////////////////
// OneBodyIntv2

OneBodyIntv2::OneBodyIntv2(const RefGaussianBasisSet&bs1,
                           const RefGaussianBasisSet&bs2) :
  OneBodyInt(bs1,bs2)
{
  c1 = int_centers_from_gbs(bs1);
  if (bs1 == bs2 || bs2.null()) {
      c2 = c1;
      same_center=1;
    }
  else {
      c2 = int_centers_from_gbs(bs2);
      same_center=0;
    }
  int_initialize_1e(0,0,c1,c2);
  int_initialize_offsets1(c1,c2);
}

OneBodyIntv2::~OneBodyIntv2()
{
  int_done_offsets1(c1,c2);
  int_done_1e();
  free_centers(c1);
  free(c1);
  if (!same_center) {
    free_centers(c2);
    free(c2);
  }
}

////////////////////////////////////////////////////////////////////////////
// GaussianOverlapIntv2

GaussianOverlapIntv2::GaussianOverlapIntv2(const RefGaussianBasisSet&bs1,
                                           const RefGaussianBasisSet&bs2) :
  OneBodyIntv2(bs1,bs2)
{
}

GaussianOverlapIntv2::~GaussianOverlapIntv2()
{
}

void
GaussianOverlapIntv2::compute_shell(int i, int j)
{
  int_shell_overlap(c1,c2,buffer_,i,j);
}

////////////////////////////////////////////////////////////////////////////
// GaussianKineticIntv2

GaussianKineticIntv2::GaussianKineticIntv2(const RefGaussianBasisSet&bs1,
                                           const RefGaussianBasisSet&bs2) :
  OneBodyIntv2(bs1,bs2)
{
}

GaussianKineticIntv2::~GaussianKineticIntv2()
{
}

void
GaussianKineticIntv2::compute_shell(int i, int j)
{
  int_shell_kinetic(c1,c2,buffer_,i,j);
}

////////////////////////////////////////////////////////////////////////////
// GaussianPointChargeIntv2

void
GaussianPointChargeIntv2::reinitialize()
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

GaussianPointChargeIntv2::GaussianPointChargeIntv2(
    const RefGaussianBasisSet&bs1,
    const RefGaussianBasisSet&bs2,
    const RefPointChargeData&dat):
  OneBodyIntv2(bs1,bs2),
  data_(dat),
  ncharge(0),
  charge(0),
  position(0)
{
}

GaussianPointChargeIntv2::~GaussianPointChargeIntv2()
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
GaussianPointChargeIntv2::compute_shell(int i,int j)
{
  int_shell_point_charge(c1,c2,buffer_,i,j,ncharge,charge,position);
}

////////////////////////////////////////////////////////////////////////////
// GaussianNuclearIntv2

GaussianNuclearIntv2::GaussianNuclearIntv2(const RefGaussianBasisSet&bs1,
                                           const RefGaussianBasisSet&bs2) :
  GaussianPointChargeIntv2(bs1, bs2,
                           new PointChargeData(bs1->molecule()->charges()))
{
}

GaussianNuclearIntv2::~GaussianNuclearIntv2()
{
}

////////////////////////////////////////////////////////////////////////////
// GaussianEfieldDotVectorIntv2

GaussianEfieldDotVectorIntv2::GaussianEfieldDotVectorIntv2(
    const RefGaussianBasisSet&bs1,
    const RefGaussianBasisSet&bs2,
    const RefEfieldDotVectorData&dat) :
  OneBodyIntv2(bs1,bs2),
  data_(dat)
{
  int biggest_shell = bs1->max_nfunction_in_shell() *
                      bs2->max_nfunction_in_shell();
  if (biggest_shell) {
      buffer3_ = new double[biggest_shell * 3];
    }
  else {
      buffer3_ = 0;
    }
}

GaussianEfieldDotVectorIntv2::~GaussianEfieldDotVectorIntv2()
{
  if (buffer3_) delete[] buffer3_;
}

void
GaussianEfieldDotVectorIntv2::compute_shell(int i,int j)
{
  int nbfi = basis1()->operator()(i).nfunction();
  int nbfj = basis2()->operator()(j).nfunction();
  int nint = nbfi*nbfj;
  int nint3 = nint*3;
  double *tmp;
  int ii,jj;

  tmp = buffer3_;
  for (ii=0;ii<nint3;ii++) *tmp++ = 0.0;

  int_accum_shell_efield(c1,c2,buffer3_,i,j,data_->position);

  tmp = buffer3_;
  for (ii=0; ii<nint; ii++) {
      buffer_[ii] = 0.0;
      for (jj=0; jj<3; jj++) {
          buffer_[ii] += *tmp++ * data_->vector[jj];
        }
    }
}

////////////////////////////////////////////////////////////////////////////
// GaussianDipoleIntv2

GaussianDipoleIntv2::GaussianDipoleIntv2(const RefGaussianBasisSet&bs1,
                                         const RefGaussianBasisSet&bs2,
                                         const RefDipoleData&dat) :
  OneBodyIntv2(bs1,bs2),
  data_(dat)
{
  if (data_.null()) {
      data_ = new DipoleData;
    }
}

GaussianDipoleIntv2::~GaussianDipoleIntv2()
{
}

void
GaussianDipoleIntv2::compute_shell(int i,int j)
{
  int nbfi = basis1()->operator()(i).nfunction();
  int nbfj = basis2()->operator()(j).nfunction();
  int nint3 = nbfi*nbfj*3;
  double *tmp;
  int ii;

  tmp = buffer_;
  for (ii=0;ii<nint3;ii++) *tmp++ = 0.0;

  int_accum_shell_dipole(c1,c2,buffer_,i,j,data_->origin);
}

////////////////////////////////////////////////////////////////////////////
// OneBodyDerivIntv2

OneBodyDerivIntv2::OneBodyDerivIntv2(const RefGaussianBasisSet&bs1,
                                     const RefGaussianBasisSet&bs2) :
  OneBodyDerivInt(bs1,bs2)
{
  c1 = int_centers_from_gbs(bs1);
  if (bs1 == bs2) {
      c2 = c1;
      same_center = 1;
    }
  else {
      c2 = int_centers_from_gbs(bs2);
      same_center = 0;
    }
  int_initialize_offsets1(c1,c2);
  intbuf = int_initialize_1e(0,1,c1,c2);
}

OneBodyDerivIntv2::~OneBodyDerivIntv2()
{
  int_done_offsets1(c1,c2);
  int_done_1e();
  free_centers(c1);
  free(c1);
  if (!same_center) {
    free_centers(c2);
    free(c2);
  }
}

void
OneBodyDerivIntv2::compute_hcore_shell(int center, int ish, int jsh)
{
  int_shell_hcore_1der(c1,c2,intbuf,ish,jsh,c1,center);

  int ni = INT_SH_NFUNC((c1),ish);
  int nj = INT_SH_NFUNC((c2),jsh);

  memcpy(buffer_, intbuf, sizeof(double)*ni*nj*3);
}

void
OneBodyDerivIntv2::compute_overlap_shell(int center, int ish, int jsh)
{
  int_shell_overlap_1der(c1,c2,intbuf,ish,jsh,c1,center);

  int ni = INT_SH_NFUNC((c1),ish);
  int nj = INT_SH_NFUNC((c2),jsh);

  memcpy(buffer_, intbuf, sizeof(double)*ni*nj*3);
}

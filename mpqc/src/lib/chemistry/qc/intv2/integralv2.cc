
#include <stdio.h>

#include <math/array/math_lib.h>
#include <math/topology/pointbag.h>

#include <chemistry/molecule/molecule.h>
#include <chemistry/qc/intv2/int_libv2.h>
#include <chemistry/qc/intv2/integralv2.h>

// initialize the transformations
static SphericalTransformV2 trans2(2);
static SphericalTransformV2 trans3(3);
static SphericalTransformV2 trans4(4);

static ISphericalTransformV2 itrans2(2);
static ISphericalTransformV2 itrans3(3);
static ISphericalTransformV2 itrans4(4);

////////////////////////////////////////////////////////////////////////////

SphericalTransformIterV2::SphericalTransformIterV2(int l, int inverse)
{
  if (l==2) {
    if (inverse) transform_ = &itrans2;
    else transform_ = &trans2;
  } else if (l==3) {
    if (inverse) transform_ = &itrans3;
    else transform_ = &trans3;
  } else if (l==4) {
    if (inverse) transform_ = &itrans4;
    else transform_ = &trans4;
  } else {
    fprintf(stderr, "SphericalTransformIterV2: cannot handle l = %d\n", l);
    abort();
  }
}

////////////////////////////////////////////////////////////////////////////
// OneBodyIntv2

OneBodyIntv2::OneBodyIntv2(const RefGaussianBasisSet&bs) :
  OneBodyInt(bs)
{
  c1 = c2 = int_centers_from_gbs(bs);
  int_initialize_1e(0,0,c1,c1);
  int_initialize_offsets1(c1,c1);
  same_center=1;
}

OneBodyIntv2::OneBodyIntv2(const RefGaussianBasisSet&bs1,
                           const RefGaussianBasisSet&bs2) :
  OneBodyInt(bs1,bs2)
{
  c1 = int_centers_from_gbs(bs1);
  c2 = int_centers_from_gbs(bs2);
  int_initialize_1e(0,0,c1,c2);
  int_initialize_offsets1(c1,c2);
  same_center=0;
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

GaussianOverlapIntv2::GaussianOverlapIntv2(const RefGaussianBasisSet&bs_) :
  OneBodyIntv2(bs_)
{
}

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

GaussianKineticIntv2::GaussianKineticIntv2(const RefGaussianBasisSet&bs_) :
  OneBodyIntv2(bs_)
{
}

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
GaussianPointChargeIntv2::init(PointBag_double*charges)
{
  ncharge = charges->length();
  
  if (ncharge) {
    position = new double*[ncharge];
    charge = new double[ncharge];
  }
  
  int i = 0;
  for (Pix pix= charges->first(); pix!=0; charges->next(pix)) {
    position[i] = new double[3];
    charge[i] = charges->get(pix);
    for (int j=0; j<3; j++) {
      position[i][j] = charges->point(pix)[j];
    }
    i++;
  }
}

GaussianPointChargeIntv2::GaussianPointChargeIntv2(PointBag_double*charges,
                                           const RefGaussianBasisSet&bs_) :
  OneBodyIntv2(bs_)
{
  init(charges);
  delete charges;
}

GaussianPointChargeIntv2::GaussianPointChargeIntv2(PointBag_double*charges,
					   const RefGaussianBasisSet&bs1,
                                           const RefGaussianBasisSet&bs2):
  OneBodyIntv2(bs1,bs2)
{
  init(charges);
  delete charges;
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

GaussianNuclearIntv2::GaussianNuclearIntv2(const RefGaussianBasisSet&bs_) :
  GaussianPointChargeIntv2(bs_->molecule()->charges(),bs_)
{
}

GaussianNuclearIntv2::GaussianNuclearIntv2(PointBag_double *charges,
                                           const RefGaussianBasisSet&bs1,
                                           const RefGaussianBasisSet&bs2) :
  GaussianPointChargeIntv2(charges,bs1,bs2)
{
}

GaussianNuclearIntv2::~GaussianNuclearIntv2()
{
}

////////////////////////////////////////////////////////////////////////////
// GaussianEfieldDotVectorIntv2

GaussianEfieldDotVectorIntv2::GaussianEfieldDotVectorIntv2(
    const RefGaussianBasisSet&bs_,
    double *p,
    double *v) :
  OneBodyIntv2(bs_)
{
  int biggest_shell = bs_->max_nfunction_in_shell();
  if (biggest_shell) {
      buffer3_ = new double[biggest_shell * biggest_shell * 3];
    }
  else {
      buffer3_ = 0;
    }
  if (p) position(p);
  if (v) vector(v);
}

GaussianEfieldDotVectorIntv2::GaussianEfieldDotVectorIntv2(
    const RefGaussianBasisSet&bs1,
    const RefGaussianBasisSet&bs2,
    double *p,
    double *v) :
  OneBodyIntv2(bs1,bs2)
{
  int biggest_shell = bs1->max_nfunction_in_shell() *
                      bs2->max_nfunction_in_shell();
  if (biggest_shell) {
      buffer3_ = new double[biggest_shell * 3];
    }
  else {
      buffer3_ = 0;
    }
  if (p) position(p);
  if (v) vector(v);
}

GaussianEfieldDotVectorIntv2::~GaussianEfieldDotVectorIntv2()
{
  if (buffer3_) delete[] buffer3_;
}

void
GaussianEfieldDotVectorIntv2::position(const double *a)
{
  for (int i=0; i<3; i++) position_[i] = a[i];
}

void
GaussianEfieldDotVectorIntv2::vector(const double *a)
{
  for (int i=0; i<3; i++) vector_[i] = a[i];
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

  int_accum_shell_efield(c1,c2,buffer3_,i,j,position_);

  tmp = buffer3_;
  for (ii=0; ii<nint; ii++) {
      buffer_[ii] = 0.0;
      for (jj=0; jj<3; jj++) {
          buffer_[ii] += *tmp++ * vector_[jj];
        }
    }
}

////////////////////////////////////////////////////////////////////////////
// GaussianDipoleIntv2

GaussianDipoleIntv2::GaussianDipoleIntv2(const RefGaussianBasisSet&bs_,
                                         const double* o) :
  OneBodyIntv2(bs_)
{
  if (o) origin(o);
}

GaussianDipoleIntv2::GaussianDipoleIntv2(const RefGaussianBasisSet&bs1,
                                         const RefGaussianBasisSet&bs2,
                                         const double* o) :
  OneBodyIntv2(bs1,bs2)
{
  if (o) origin(o);
}

GaussianDipoleIntv2::~GaussianDipoleIntv2()
{
}

void
GaussianDipoleIntv2::origin(const double *a)
{
  for (int i=0; i<3; i++) origin_[i] = a[i];
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

  int_accum_shell_dipole(c1,c2,buffer_,i,j,origin_);
}

/////////////////////////////////////////////////////////////////////////////

#define CLASSNAME IntegralV2
#define PARENTS public Integral
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
IntegralV2::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = Integral::_castdown(cd);
  return do_castdowns(casts,cd);
}

IntegralV2::IntegralV2()
{
}

IntegralV2::IntegralV2(StateIn& s) :
  Integral(s)
{
}

IntegralV2::IntegralV2(const RefKeyVal& k) :
  Integral(k)
{
}

void
IntegralV2::save_data_state(StateOut& s)
{
  Integral::save_data_state(s);
}

CartesianIter *
IntegralV2::new_cartesian_iter(int l)
{
  return new CartesianIterV2(l);
}

RedundantCartesianIter *
IntegralV2::new_redundant_cartesian_iter(int l)
{
  return new RedundantCartesianIterV2(l);
}

RedundantCartesianSubIter *
IntegralV2::new_redundant_cartesian_sub_iter(int l)
{
  return new RedundantCartesianSubIterV2(l);
}

SphericalTransformIter *
IntegralV2::new_spherical_transform_iter(int l, int inv)
{
  return new SphericalTransformIterV2(l, inv);
}

RefOneBodyInt
IntegralV2::overlap_int(const RefGaussianBasisSet& gbs)
{
  return new GaussianOverlapIntv2(gbs);
}

RefOneBodyInt
IntegralV2::overlap_int(const RefGaussianBasisSet& gbs1,
                        const RefGaussianBasisSet& gbs2)
{
  return new GaussianOverlapIntv2(gbs1, gbs2);
}

RefOneBodyInt
IntegralV2::kinetic_int(const RefGaussianBasisSet& gbs)
{
  return new GaussianKineticIntv2(gbs);
}

RefOneBodyInt
IntegralV2::kinetic_int(const RefGaussianBasisSet& gbs1,
                        const RefGaussianBasisSet& gbs2)
{
  return new GaussianKineticIntv2(gbs1, gbs2);
}

RefOneBodyInt
IntegralV2::point_charge_int(PointBag_double* pb,
                             const RefGaussianBasisSet& gbs)
{
  return new GaussianPointChargeIntv2(pb, gbs);
}

RefOneBodyInt
IntegralV2::point_charge_int(PointBag_double *pb,
                             const RefGaussianBasisSet& gbs1,
                             const RefGaussianBasisSet& gbs2)
{
  return new GaussianPointChargeIntv2(pb, gbs1, gbs2);
}

RefOneBodyInt
IntegralV2::nuclear_int(const RefGaussianBasisSet& gbs)
{
  return new GaussianNuclearIntv2(gbs);
}

RefOneBodyInt
IntegralV2::nuclear_int(PointBag_double *pb,
                        const RefGaussianBasisSet& gbs1,
                        const RefGaussianBasisSet& gbs2)
{
  return new GaussianNuclearIntv2(pb, gbs1, gbs2);
}

RefOneBodyInt
IntegralV2::efield_dot_vector_int(const RefGaussianBasisSet& gbs,
                                 double *position,
                                 double *vector)
{
  return new GaussianEfieldDotVectorIntv2(gbs, position, vector);
}

RefOneBodyInt
IntegralV2::efield_dot_vector_int(const RefGaussianBasisSet& gbs1,
                                 const RefGaussianBasisSet& gbs2,
                                 double *position,
                                 double *vector)
{
  return new GaussianEfieldDotVectorIntv2(gbs1, gbs2, position, vector);
}

RefOneBodyInt
IntegralV2::dipole_int(const RefGaussianBasisSet& gbs,
                      double *origin)
{
  return new GaussianDipoleIntv2(gbs, origin);
}

RefOneBodyInt
IntegralV2::dipole_int(const RefGaussianBasisSet& gbs1,
                       const RefGaussianBasisSet& gbs2,
                       double *origin)
{
  return new GaussianDipoleIntv2(gbs1, gbs2, origin);
}

RefTwoBodyInt
IntegralV2::two_body_int(const RefGaussianBasisSet& gbs)
{
  return new TwoBodyIntV2(gbs);
}

RefTwoBodyInt
IntegralV2::two_body_int(const RefGaussianBasisSet& gbs1,
                         const RefGaussianBasisSet& gbs2,
                         const RefGaussianBasisSet& gbs3,
                         const RefGaussianBasisSet& gbs4)
{
  return new TwoBodyIntV2(gbs1, gbs2, gbs3, gbs4);
}

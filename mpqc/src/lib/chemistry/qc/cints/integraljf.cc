
#include <stdio.h>
#include <chemistry/qc/cints/cints.h>
#include <chemistry/qc/cints/integraljf.h>

static double
int_pow(double a, int p)
{
  switch (p) {
  case 0:
    return 1.0;
    break;

  case 1:
    return a;
    break;

  case 2:
    return a*a;
    break;

  case 3:
    return a*a*a;
    break;

  case 4:
    return a*a*a*a;
    break;

  case 5:
    return a*a*a*a*a;
    break;

  default:
    {
      register int i;
      double b = 1.0;
    
      for(i=0; i<p; i++) b = b*a;
      return b;
    }
  }
}

////////////////////////////////////////////////////////////////////////////
// OneBodyIntJF

#define MAXAM 5

int OneBodyIntJF::maxfact = 100;

void
OneBodyIntJF::init()
{
  int i,j;
  
  // find the max angular momentum and set up constant arrays
  maxam=MAXAM;

  // lci is binomial coefficient, l choose i, in one-d array
  if (lci)
    delete[] lci;

  lci = new int[ioff(maxam+1)];
  lci[0]=0;

  for (i=1; i <= maxam; i++)
    lci[ioff(i)] = lci[ioff(i)-1] = 1;

  for (i=2; i <= maxam; i++)
    for (j=1; j < i; j++)
      lci[ioff(i)+j] = lci[ioff(i-1)+j-1]+lci[ioff(i-1)+j];

  lci[ioff(maxam+1)] = 1;
  
  // df[i] gives (i-1)!!, so that (-1)!! is defined...
  // we shouldn't need both this and lci with the range needed on df[]
  if (df)
    delete[] df;
  
  df = new double[maxfact*2];
  
  df[0] = 1.0;
  df[1] = 1.0;
  df[2] = 1.0;
  for (i=3; i < maxfact*2; i++)
    df[i] = (i-1)*df[i-2];

  // (2i-1)!! a useful thing, num_ser in the integral expansion
  if (num_ser)
    delete[] num_ser;

  num_ser = new unsigned int[maxam+1];
  
  num_ser[0] = 1;
  for (i=1; i <= maxam; i++)
    num_ser[i] = (2*i-1)*num_ser[i-1];
}

OneBodyIntJF::OneBodyIntJF(const RefGaussianBasisSet&bs, OneBodyIntIter *it) :
  lci(0), num_ser(0), df(0),
  OneBodyInt(bs,it)
{
  init();
}

OneBodyIntJF::OneBodyIntJF(const RefGaussianBasisSet&bs1,
                           const RefGaussianBasisSet&bs2,
                           OneBodyIntIter *it) :
  lci(0), num_ser(0), df(0),
  OneBodyInt(bs1,bs2,it)
{
  init();
}

OneBodyIntJF::~OneBodyIntJF()
{
  if (lci) {
    delete[] lci;
    lci=0;
  }
  if (df) {
    delete[] df;
    df=0;
  }
  if (num_ser) {
    delete[] num_ser;
    num_ser=0;
  }
}

double
OneBodyIntJF::overlap_int(
  double a1, int l1, int m1, int n1, double norm1,
  double a2, int l2, int m2, int n2, double norm2,
  double ab2, double gam,
  Point& PA, Point& PB, int am)
{
  double Ix, Iy, Iz;
  double I;
  int i, j, k, l;
  int imax, jmax, kmax;
  double tval, tval1, tval2 ;
  double norm_fact ;

  norm_fact = norm1*norm2; 

  tval1 = 2*gam;
  imax = (l1+l2)/2;
  Ix = 0.0;
  for (i=0; i <= imax; i++) {
    tval = f_n(i*2, l1, l2, PA[0], PB[0]);
    tval2 = int_pow(tval1, i);
    //tval2 = pow(tval1, (double) i);
    Ix += tval*(num_ser[i])/(tval2);
  }

  jmax = (m1+m2)/2;
  Iy = 0.0;
  for (j=0; j <= jmax; j++) {
    tval = f_n(j*2, m1, m2, PA[1], PB[1]);
    tval2 = int_pow(tval1, j);
    //tval2 = pow(tval1, (double) j);
    Iy += tval*num_ser[j]/(tval2);
  }

  kmax = (n1+n2)/2;
  Iz = 0.0;
  for (k=0; k <= kmax; k++) {
    tval = f_n(k*2, n1, n2, PA[2], PB[2]);
    tval2 = int_pow(tval1, k);
    //tval2 = pow(tval1, (double) k);
    Iz += tval*num_ser[k]/(tval2);
  }
 
  I = exp(-1*a1*a2*ab2/gam)*Ix*Iy*Iz*sqrt(M_PI/gam)*(M_PI/gam);

  return I*norm_fact;
}

double
OneBodyIntJF::f_n(int k, int l1, int l2, double A, double B)
{
  double sum = 0.0;
  int i, j, itmp;

  int* iol1 = lci + ioff(l1);
  int* iol2 = lci + ioff(l2);
  
  for (i=0; i <= l1; i++) {
    j = k-i;
    if (j > l2)
      continue;
    
    itmp = iol1[i] * iol2[j];
    sum += itmp * int_pow(A, (l1-i)) * int_pow(B, (l2-j));
  }

  return sum;
}

#if 0
////////////////////////////////////////////////////////////////////////////
// GaussianKineticIntv2

GaussianKineticIntJF::GaussianKineticIntJF(const RefGaussianBasisSet&bs_,
                                           OneBodyIntIter *it):
  OneBodyIntJF(bs_,it)
{
}

GaussianKineticIntJF::GaussianKineticIntJF(const RefGaussianBasisSet&bs1,
                                           const RefGaussianBasisSet&bs2,
                                           OneBodyIntIter *it):
  OneBodyIntJF(bs1,bs2,it)
{
}

GaussianKineticIntJF::~GaussianKineticIntJF()
{
}

void GaussianKineticIntJF::compute_shell(int i, int j, double * buf)
{
  //int_shell_kinetic(c1,c2,buf,i,j);
}

////////////////////////////////////////////////////////////////////////////
// GaussianPointChargeIntv2

void
GaussianPointChargeIntJF::init(PointBag_double*charges)
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

GaussianPointChargeIntJF::GaussianPointChargeIntJF(PointBag_double*charges,
					   const RefGaussianBasisSet&bs_,
                                           OneBodyIntIter *it) :
  OneBodyIntJF(bs_,it)
{
  init(charges);
  delete charges;
}

GaussianPointChargeIntJF::GaussianPointChargeIntJF(PointBag_double*charges,
					   const RefGaussianBasisSet&bs1,
					   const RefGaussianBasisSet&bs2,
                                           OneBodyIntIter *it) :
  OneBodyIntJF(bs1,bs2,it)
{
  init(charges);
  delete charges;
}

GaussianPointChargeIntJF::~GaussianPointChargeIntJF()
{
  for (int i=0; i<ncharge; i++) {
      delete position[i];
    }
  if (ncharge) {
      delete[] charge;
      delete[] position;
    }
}

void GaussianPointChargeIntJF::compute_shell(int i,int j,double*buf)
{
  //int_shell_point_charge(c1,c2,buf,i,j,ncharge,charge,position);
}

////////////////////////////////////////////////////////////////////////////
// GaussianNuclearIntJF

GaussianNuclearIntJF::GaussianNuclearIntJF(const RefGaussianBasisSet&bs_,
                                           OneBodyIntIter *it) :
  GaussianPointChargeIntJF(bs_->molecule()->charges(),bs_,it)
{
}

GaussianNuclearIntJF::GaussianNuclearIntJF(PointBag_double *charges,
                                           const RefGaussianBasisSet&bs1,
                                           const RefGaussianBasisSet&bs2,
                                           OneBodyIntIter *it):
  GaussianPointChargeIntJF(charges,bs1,bs2,it)
{
}

GaussianNuclearIntJF::~GaussianNuclearIntJF()
{
}
#endif


#include <stdio.h>
#include <chemistry/qc/cints/cints.h>
#include <chemistry/qc/cints/integraljf.h>

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
OneBodyIntJF::overlap_int(int l1, int m1, int n1, int l2, int m2, int n2,
                          double PA[3], double PB[3],
                          double gam)
{
  int i, i2, ii, jj, imax, itmp;
  double Ix, I=1.0;
  double oo2gam = 0.5/gam;
  double gtmp;

  int al1[] = { l1, m1, n1 };
  int al2[] = { l2, m2, n2 };

  for (int x=0; x < 3; x++) {
    int l1x = al1[x];
    int l2x = al2[x];

    int *lcip1 = lci + (l1x*(l1x+1)>>1);
    int *lcip2 = lci + (l2x*(l2x+1)>>1);

    double da = PA[x];
    double db = PB[x];
  
    Ix = 0.0;
    gtmp=1.0;
    imax = (l1x+l2x)>>1;

    for (i=0,i2=0; i <= imax; i++, gtmp *= oo2gam, i2 += 2) {
      double sum = 0.0;

      for (ii=0; ii <= l1x; ii++) {
        jj = i2-ii;
        if (jj < 0)
          break;
        else if (jj > l2x)
          continue;
    
        itmp = lcip1[ii] * lcip2[jj];
        if (!itmp)
          continue;
    
        sum += itmp * pow(da, (long int)(l1x-ii)) *
                      pow(db, (long int)(l2x-jj));
      }

      Ix += sum * num_ser[i] * gtmp;
    }

    I *= Ix;
  }

  return I;
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

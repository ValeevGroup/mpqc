
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
OneBodyIntJF::overlap_intx(int l1, int l2, double da, double db, double oo2gam)
{
  int i, i2, ii, jj, imax, itmp;
  int *lcip1, *lcip2;
  double Ix = 0.0;
  double gtmp = 1.0;

  lcip1 = lci + (l1*(l1+1)>>1);
  lcip2 = lci + (l2*(l2+1)>>1);

  imax = (l1+l2)>>1;

  for (i=0,i2=0; i <= imax; i++, i2++, i2++, gtmp *= oo2gam) {
    double sum = 0.0;

    ii=0; jj=i2;
    while (jj > l2) {
      ii++;
      jj--;
    }

    for (; (jj > -1) && (ii <= l1); ii++, jj--) {
    
      itmp = lcip1[ii] * lcip2[jj] * num_ser[i];
    
      sum += itmp * pow(da, (long int)(l1-ii)) *
                    pow(db, (long int)(l2-jj));
    }

    Ix += sum * gtmp;
  }

  return Ix;
}

double
OneBodyIntJF::overlap_int(int l1, int m1, int n1, int l2, int m2, int n2,
                          double PA[3], double PB[3],
                          double oo2gam)
{
  return overlap_intx(l1,l2,PA[0],PB[0],oo2gam) *
         overlap_intx(m1,m2,PA[1],PB[1],oo2gam) *
         overlap_intx(n1,n2,PA[2],PB[2],oo2gam);
}

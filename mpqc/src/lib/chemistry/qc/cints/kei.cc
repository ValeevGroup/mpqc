
#include <chemistry/molecule/localdef.h>
#include <chemistry/qc/cints/cints.h>
#include <chemistry/qc/cints/integraljf.h>

////////////////////////////////////////////////////////////////////////////
// GaussianOverlapIntJF

GaussianKineticIntJF::GaussianKineticIntJF(const RefGaussianBasisSet&bs_,
                                           OneBodyIntIter *it) :
  OneBodyIntJF(bs_,it)
{
}

GaussianKineticIntJF::GaussianKineticIntJF(const RefGaussianBasisSet&bs1,
                                           const RefGaussianBasisSet&bs2,
                                           OneBodyIntIter *it) :
  OneBodyIntJF(bs1,bs2,it)
{
}

GaussianKineticIntJF::~GaussianKineticIntJF()
{
}

void
GaussianKineticIntJF::compute_shell(int i, int j, double * buf)
{
  int isiz = bs1->operator()(i).nfunction();
  int jsiz = bs2->operator()(j).nfunction();

  memset(buf,0,isiz*jsiz*sizeof(double));

  AtomicCenter& ai = bs1->molecule()->atom(bs1->shell_to_center(i));
  AtomicCenter& aj = bs2->molecule()->atom(bs2->shell_to_center(j));

  GaussianShell& gsi = bs1->operator()(i);
  GaussianShell& gsj = bs2->operator()(j);

  double ab2 = dist(ai.point(), aj.point());
  ab2 *= ab2;

  // loop over primitives
  for (int pi=0; pi < gsi.nprimitive(); pi++) {
    double a1 = gsi.exponent(pi);

    for (int pj=0; pj < gsj.nprimitive(); pj++) {
      double a2 = gsj.exponent(pj);
      double gam = a1+a2;

      Point P, PA, PB;

      P[0] = (ai[0]*a1 + aj[0]*a2)/gam;
      P[1] = (ai[1]*a1 + aj[1]*a2)/gam;
      P[2] = (ai[2]*a1 + aj[2]*a2)/gam;

      PA[0] = P[0] - ai[0];
      PA[1] = P[1] - ai[1];
      PA[2] = P[2] - ai[2];

      PB[0] = P[0] - aj[0];
      PB[1] = P[1] - aj[1];
      PB[2] = P[2] - aj[2];

      // loop over general contractions
      int ioffset=0;
      for (int ci=0; ci < gsi.ncontraction(); ci++) {
        double inorm = gsi.coefficient_unnorm(ci,pi);
        int ami = gsi.am(ci);
        
        int joffset=0;
        for (int cj=0; cj < gsj.ncontraction(); cj++) {
          double jnorm = gsj.coefficient_unnorm(cj,pj);
          int amj = gsj.am(cj);

          int ij=ioffset;
          for (int ii=0; ii <= ami; ii++) {
            int l1 = ami-ii;

            for (int jj=0; jj <= ii; jj++,ij++) {
              int m1 = ii - jj;
              int n1 = jj;

              int ijkl=ij*jsiz+joffset;
              for (int kk=0; kk <= amj; kk++) {
                int l2 = amj-kk;

                for (int ll=0; ll <= kk; ll++,ijkl++) {
                  int m2 = kk-ll;
                  int n2 = ll;

                  int am = l1+m1+n1+l2+m2+n2;

                  buf[ijkl] += ke_int(a1, l1, m1, n1, inorm,
                                      a2, l2, m2, n2, jnorm,
                                      ab2, gam, PA, PB, am);
                }
              }
            }
          }
          joffset += ioff(amj+1);
        }
        ioffset += ioff(ami+1);
      }
    }
  }
}

double
GaussianKineticIntJF::ke_int(double a1, int l1, int m1, int n1, double norm1,
                             double a2, int l2, int m2, int n2, double norm2,
                             double ab2, double gam,
                             Point& PA, Point& PB, int am)
{
  double Ix, Iy, Iz;
  double I1 = 0.0, I2 = 0.0, I3 = 0.0, I4 = 0.0 ;

  I2 = overlap_int(a1, l1+1, m1, n1, norm1, a2, l2+1, m2, n2, norm2, 
                   ab2, gam, PA, PB, am+2);
  I1 = overlap_int(a1, l1-1, m1, n1, norm1, a2, l2-1, m2, n2, norm2, 
                   ab2, gam, PA, PB, am-2);
  I3 = overlap_int(a1, l1+1, m1, n1, norm1, a2, l2-1, m2, n2, norm2, 
                   ab2, gam, PA, PB, am);
  I4 = overlap_int(a1, l1-1, m1, n1, norm1, a2, l2+1, m2, n2, norm2, 
                   ab2, gam, PA, PB, am);
 
  Ix = (l1*l2*I1/2.0 + 2*a1*a2*I2 - a1*l2*I3 - a2*l1*I4);

  I1 = 0.0;
  I2 = 0.0;
  I3 = 0.0;
  I4 = 0.0;

  I2 = overlap_int(a1, l1, m1+1, n1, norm1, a2, l2, m2+1, n2, norm2, 
                   ab2, gam, PA, PB, am+2);
  I1 = overlap_int(a1, l1, m1-1, n1, norm1, a2, l2, m2-1, n2, norm2, 
                   ab2, gam, PA, PB, am-2);
  I3 = overlap_int(a1, l1, m1+1, n1, norm1, a2, l2, m2-1, n2, norm2, 
                   ab2, gam, PA, PB, am);
  I4 = overlap_int(a1, l1, m1-1, n1, norm1, a2, l2, m2+1, n2, norm2, 
                   ab2, gam, PA, PB, am);
  
  Iy = (m1*m2*I1/2.0 + 2*a1*a2*I2 - a1*m2*I3 - a2*m1*I4);

  I1 = 0.0;
  I2 = 0.0;
  I3 = 0.0;
  I4 = 0.0;

  I2 = overlap_int(a1, l1, m1, n1+1, norm1, a2, l2, m2, n2+1, norm2, 
                   ab2, gam, PA, PB, am+2);
  I1 = overlap_int(a1, l1, m1, n1-1, norm1, a2, l2, m2, n2-1, norm2, 
                   ab2, gam, PA, PB, am-2);
  I3 = overlap_int(a1, l1, m1, n1+1, norm1, a2, l2, m2, n2-1, norm2, 
                   ab2, gam, PA, PB, am);
  I4 = overlap_int(a1, l1, m1, n1-1, norm1, a2, l2, m2, n2+1, norm2, 
                   ab2, gam, PA, PB, am);
  
  Iz = (n1*n2*I1/2.0 + 2*a1*a2*I2 - a1*n2*I3 - a2*n1*I4);

  return((Ix+Iy+Iz));
}


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
      double oogam = 1.0/gam;

      double prefact = exp(-a1*a2*ab2*oogam)*pow(M_PI*oogam,1.5);

      double Px, Py, Pz;
      Px = (ai[0]*a1 + aj[0]*a2)*oogam;
      Py = (ai[1]*a1 + aj[1]*a2)*oogam;
      Pz = (ai[2]*a1 + aj[2]*a2)*oogam;

      double PA[3], PB[3];

      PA[0] = Px - ai[0];
      PA[1] = Py - ai[1];
      PA[2] = Pz - ai[2];

      PB[0] = Px - aj[0];
      PB[1] = Py - aj[1];
      PB[2] = Pz - aj[2];

      oogam *= 0.5;
      
      // loop over general contractions
      int ioffset=0;
      for (int ci=0; ci < gsi.ncontraction(); ci++) {
        double inorm = gsi.coefficient_unnorm(ci,pi);
        int ami = gsi.am(ci);
        
        int joffset=0;
        for (int cj=0; cj < gsj.ncontraction(); cj++) {
          double jnorm = gsj.coefficient_unnorm(cj,pj);
          int amj = gsj.am(cj);

          double norm_prefact = inorm*jnorm*prefact;
          
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

                  buf[ijkl] += norm_prefact*ke_int(l1, m1, n1, l2, m2, n2,
                                                   PA, PB, oogam, a1, a2);
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
GaussianKineticIntJF::ke_int(int l1, int m1, int n1, int l2, int m2, int n2,
                             double PA[3], double PB[3],
                             double oo2gam, double a1, double a2)
{
  double I = 0;
  double I1, I2, I3, I4, Ilmn;

  double pax = PA[0];
  double pbx = PB[0];
  
  Ilmn = overlap_intx(m1,m2,PA[1],PB[1],oo2gam) *
         overlap_intx(n1,n2,PA[2],PB[2],oo2gam);
  I2 = Ilmn * overlap_intx(l1+1, l2+1, pax, pbx, oo2gam);
  I3 = Ilmn * overlap_intx(l1+1, l2-1, pax, pbx, oo2gam);
  I1 = Ilmn * overlap_intx(l1-1, l2-1, pax, pbx, oo2gam);
  I4 = Ilmn * overlap_intx(l1-1, l2+1, pax, pbx, oo2gam);
 
  I += (0.5*l1*l2*I1 + 2*a1*a2*I2 - a1*l2*I3 - a2*l1*I4);

  pax = PA[1];
  pbx = PB[1];
  
  Ilmn = overlap_intx(l1,l2,PA[0],PB[0],oo2gam) *
         overlap_intx(n1,n2,PA[2],PB[2],oo2gam);
  I2 = Ilmn * overlap_intx(m1+1, m2+1, pax, pbx, oo2gam);
  I3 = Ilmn * overlap_intx(m1+1, m2-1, pax, pbx, oo2gam);
  I1 = Ilmn * overlap_intx(m1-1, m2-1, pax, pbx, oo2gam);
  I4 = Ilmn * overlap_intx(m1-1, m2+1, pax, pbx, oo2gam);
  
  I += (0.5*m1*m2*I1 + 2*a1*a2*I2 - a1*m2*I3 - a2*m1*I4);

  pax = PA[2];
  pbx = PB[2];
  
  Ilmn = overlap_intx(l1,l2,PA[0],PB[0],oo2gam) *
         overlap_intx(m1,m2,PA[1],PB[1],oo2gam);
  I2 = Ilmn * overlap_intx(n1+1, n2+1, pax, pbx, oo2gam);
  I3 = Ilmn * overlap_intx(n1+1, n2-1, pax, pbx, oo2gam);
  I1 = Ilmn * overlap_intx(n1-1, n2-1, pax, pbx, oo2gam);
  I4 = Ilmn * overlap_intx(n1-1, n2+1, pax, pbx, oo2gam);
  
  I += (0.5*n1*n2*I1 + 2*a1*a2*I2 - a1*n2*I3 - a2*n1*I4);

  return I;
}

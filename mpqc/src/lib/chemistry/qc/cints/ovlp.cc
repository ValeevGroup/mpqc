
#include <chemistry/molecule/localdef.h>
#include <chemistry/qc/cints/cints.h>
#include <chemistry/qc/cints/integraljf.h>

////////////////////////////////////////////////////////////////////////////
// GaussianOverlapIntJF

GaussianOverlapIntJF::GaussianOverlapIntJF(const RefGaussianBasisSet&bs_,
                                           OneBodyIntIter *it) :
  OneBodyIntJF(bs_,it)
{
}

GaussianOverlapIntJF::GaussianOverlapIntJF(const RefGaussianBasisSet&bs1,
                                           const RefGaussianBasisSet&bs2,
                                           OneBodyIntIter *it) :
  OneBodyIntJF(bs1,bs2,it)
{
}

GaussianOverlapIntJF::~GaussianOverlapIntJF()
{
}

void
GaussianOverlapIntJF::compute_shell(int i, int j, double * buf)
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
      
      double PA[3], PB[3];

      double px, py, pz;
      px = (ai[0]*a1 + aj[0]*a2)*oogam;
      py = (ai[1]*a1 + aj[1]*a2)*oogam;
      pz = (ai[2]*a1 + aj[2]*a2)*oogam;

      PA[0] = px - ai[0];
      PA[1] = py - ai[1];
      PA[2] = pz - ai[2];

      PB[0] = px - aj[0];
      PB[1] = py - aj[1];
      PB[2] = pz - aj[2];

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

                  buf[ijkl] += norm_prefact*overlap_int(l1, m1, n1, l2, m2, n2,
                                                        PA, PB, gam);
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


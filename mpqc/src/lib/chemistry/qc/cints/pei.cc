
#include <math.h>

#include <chemistry/molecule/localdef.h>
#include <chemistry/qc/cints/cints.h>
#include <chemistry/qc/cints/integraljf.h>

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

  vint = new double***[maxam];
  vint_done = new char**[maxam];
  
  for (i=0; i < maxam; i++) {
    vint[i] = new double**[maxam];
    vint_done[i] = new char*[maxam];

    for (int j=0; j < maxam; j++) {
      vint_done[i][j] = new char[2*maxam];

      int vijlen = 2*maxam-i-j;
      vint[i][j] = new double*[vijlen];

      for (int k=0; k < vijlen; k++) {
        vint[i][j][k] = new double[ioffset(i+1)*ioffset(j+1)];
      }
    }
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
  int i;
  for (i=0; i < ncharge; i++)
    delete position[i];

  if (ncharge) {
    delete[] charge;
    delete[] position;
    charge=0;
    position=0;
  }
  
  if (vint && vint_done) {
    for (i=0; i < maxam; i++) {
      for (int j=0; j < maxam; j++) {
        int vijlen = 2*maxam-i-j;
        for (int k=0; k < vijlen; k++)
          delete[] vint[i][j][k];

        delete[] vint_done[i][j];
        delete[] vint[i][j];
      }
      delete[] vint[i];
      delete[] vint_done[i];
    }

    delete[] vint;
    delete[] vint_done;

    vint=0;
    vint_done=0;
  }
}

void
GaussianPointChargeIntJF::compute_shell(int i,int j,double*buf)
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
      double oo2gam = 0.5*oogam;

      double K_peint = 2*M_PI*oogam*exp(-a1*a2*ab2*oogam);

      double P[3], PA[3], PB[3];

      P[0] = (ai[0]*a1 + aj[0]*a2)*oogam;
      P[1] = (ai[1]*a1 + aj[1]*a2)*oogam;
      P[2] = (ai[2]*a1 + aj[2]*a2)*oogam;

      PA[0] = P[0] - ai[0];
      PA[1] = P[1] - ai[1];
      PA[2] = P[2] - ai[2];

      PB[0] = P[0] - aj[0];
      PB[1] = P[1] - aj[1];
      PB[2] = P[2] - aj[2];

      // loop over general contractions
      int ioffst=0;
      for (int ci=0; ci < gsi.ncontraction(); ci++) {
        double inorm = gsi.coefficient_unnorm(ci,pi);
        int ami = gsi.am(ci);
        
        int joffset=0;
        for (int cj=0; cj < gsj.ncontraction(); cj++) {
          double jnorm = gsj.coefficient_unnorm(cj,pj);
          int amj = gsj.am(cj);

          double K_norm = inorm*jnorm*K_peint;
          
          init_peint(ami, amj, pi, pj, gam, oo2gam, K_norm, P, PA, PB);

          int ij=ioffst;
          int index=0;
          for (int ii=0; ii <= ami; ii++) {
            int l1 = ami-ii;

            for (int jj=0; jj <= ii; jj++,ij++) {
              int m1 = ii - jj;
              int n1 = jj;

              int ijkl=ij*jsiz+joffset;
              for (int kk=0; kk <= amj; kk++) {
                int l2 = amj-kk;

                for (int ll=0; ll <= kk; ll++,ijkl++,index++) {
                  int m2 = kk-ll;
                  int n2 = ll;

                  buf[ijkl] += -vint[ami][amj][0][index];
                }
              }
            }
          }
          joffset += ioffset(amj+1);
        }
        ioffst += ioffset(ami+1);
      }
    }
  }
}

void
GaussianPointChargeIntJF::init_peint(int am1, int am2, int pi, int pj,
                                     double gam, double oo2gam, double K_peint,
                                     double P[3], double PA[3], double PB[3])
{
  int pt;
  int am;
  double *I0, *I1, *I2, *I3, *I4, *I5;
  int ij;
  int L0[3];
  int L1[3];
  int ii, kk, jj, ll;
  int xyz, c12, c13, c14;
  int t20, t30, t40;
  int a, c;
  int i, j, k;
  int T2max;
  int T1, T2, T3, T4;
  double t;
  double pc2;
  double sum;
  double term;
  double F[20];
  double abc[3][3]; // PA, PB, PC
  double *via1a2m = vint[am1][am2][0];
  double *via1a2mi;
  
  am = am1+am2;
  abc[0][0] = PA[0];
  abc[0][1] = PA[1];
  abc[0][2] = PA[2];
  abc[1][0] = PB[0];
  abc[1][1] = PB[1];
  abc[1][2] = PB[2];

  // zero the storage for targets
  memset((char *)via1a2m, 0, (ioffset(am1+1)*ioffset(am2+1)*sizeof(double)));

  // for each center
  for (pt=0; pt < ncharge; pt++) {
    via1a2mi = via1a2m;

    abc[2][0] = P[0] - position[pt][0];
    abc[2][1] = P[1] - position[pt][1];
    abc[2][2] = P[2] - position[pt][2];

    pc2 = abc[2][0]*abc[2][0]+abc[2][1]*abc[2][1]+abc[2][2]*abc[2][2];
    t = gam*pc2;

    // zero out the flags for all intermediates below this target class
    for (i=0; i <= am1; i++) {
      for (j=0; j <= am2; j++) {
        for (k=0; k <= (am1+am2-i-j); k++) {
          vint_done[i][j][k] = 0;
        }
      }
    }

    // convert to aux nuc attr int
    calc_f(F, am, t);
    for (j=0; j <= am; j++)
      F[j] *= K_peint;

    if (am1==0 && am2==0) {
      *via1a2m += F[0]*charge[pt];
    } else {
      if (am1) {
        I0 = do_pe_doublet(abc, F, am1-1, am2, 0, pt, oo2gam);
        I1 = do_pe_doublet(abc, F, am1-1, am2, 1, pt, oo2gam);

        if (am1>1) {
          I2 = do_pe_doublet(abc, F, am1-2, am2, 0, pt, oo2gam);
          I3 = do_pe_doublet(abc, F, am1-2, am2, 1, pt, oo2gam);
        }
        if (am2) {
          I4 = do_pe_doublet(abc, F, am1-1, am2-1, 0, pt, oo2gam);
          I5 = do_pe_doublet(abc, F, am1-1, am2-1, 1, pt, oo2gam);
        }

        // create all am components of si
        ij = 0;

        T2max = (am1*(am1+1)*(am2+1)*(am2+2))>>2;
        T1 = T2 = T3 = T4 = 0;

        for (ii = 0; ii <= am1; ii++) {
          L0[0] = am1 - ii;

          for (jj = 0; jj <= ii; jj++) {
            L0[1] = ii - jj;
            L0[2] = jj ;

            if (L0[0])
              xyz=0;
            else if(L0[1])
              xyz=1;
            else if(L0[2])
              xyz=2;

            // create all am components of sj
            for (kk = 0; kk <= am2; kk++) {
              L1[0] = am2 - kk;

              for (ll = 0; ll <= kk; ll++) {
                L1[1] = kk - ll;
                L1[2] = ll ;

                if (T2==T2max) {
                  L0[xyz] = L0[xyz] - 1;
                  am1 = am1 - 1;
                  a = c = 0;

                  if (am1) {
                    i = am1 - L0[0];
                    a = i + ioffset(i) - L0[1];
                  }

                  if (am2) {
                    i = am2-L1[0];
                    c = i + ioffset(i) - L1[1];
                  }

                  T2 = a*ioffset(am2+1)+c;

                  L0[xyz] = L0[xyz] - 1;
                  am1 = am1 - 1;
                  a = c = 0;

                  if (am1) {
                    i = am1-L0[0];
                    a = i + ioffset(i) - L0[1];
                  }

                  if (am2) {
                    i = am2-L1[0];
                    c = i + ioffset(i) - L1[1];
                  }

                  T3 = a*ioffset(am2+1)+c;

                  L0[xyz] = L0[xyz] + 1;
                  am1 = am1 + 1;
                  L1[0] = L1[0] - 1;
                  am2 = am2 - 1;
                  a = c = 0;

                  if (am1) {
                    i = am1-L0[0];
                    a = i + ioffset(i) - L0[1];
                  }

                  if (am2) {
                    i = am2-L1[0];
                    c = i + ioffset(i) - L1[1];
                  }

                  T4 = a*ioffset(am2+1)+c;
                  L1[0] = L1[0] + 1;
                  am2 = am2 + 1;
                  L0[xyz] = L0[xyz] + 1;
                  am1 = am1 + 1;
                }

                term = abc[0][xyz]*I0[T2] - abc[2][xyz]*I1[T2];
                T2++;
                if (L0[xyz] > 1) {
                  term += oo2gam*(L0[xyz]-1)*(I2[T3] - I3[T3]);
                  T3++;
                }
                if (L1[xyz]) {
                  term += oo2gam*(L1[xyz])*(I4[T4] - I5[T4]);
                  T4++;
                }

                term *= charge[pt];
                *via1a2mi += term;
                ij++;
                via1a2mi++;
              }
            }
          }
        }
      } else if (am2) {
        I0 = do_pe_doublet(abc, F, am1, am2-1, 0, pt, oo2gam);
        I1 = do_pe_doublet(abc, F, am1, am2-1, 1, pt, oo2gam);

        if (am2>1) {
          I2 = do_pe_doublet(abc, F, am1, am2-2, 0, pt, oo2gam);
          I3 = do_pe_doublet(abc, F, am1, am2-2, 1, pt, oo2gam);
        }
        ij = 0;

        T2max = am2*(am2+1)/2;
        T1 = T2 = T3 = T4 = 0;
        // create all am components of si
        L0[0] = 0;
        L0[1] = 0;
        L0[2] = 0;
        
        // create all am components of sj
        for (kk = 0; kk <= am2; kk++) {
          L1[0] = am2 - kk;

          for (ll = 0; ll <= kk; ll++) {
            L1[1] = kk - ll;
            L1[2] = ll ;

            if (L1[0])
              xyz = 0;
            else if (L1[1])
              xyz = 1;
            else if (L1[2])
              xyz = 2;

            if (T2==T2max) {
              L1[xyz] = L1[xyz] - 1;
              am2 = am2 - 1;
              a = c = 0;
              if (am1) {
                i = am1-L0[0];
                a = i + ioffset(i) - L0[1];
              }
              if (am2) {
                i = am2-L1[0];
                c = i + ioffset(i) - L1[1];
              }
              T2 = a*ioffset(am2+1)+c;

              L1[xyz] = L1[xyz] - 1;
              am2 = am2 - 1;
              a = c = 0;
              if (am1) {
                i = am1-L0[0];
                a = i + ioffset(i) - L0[1];
              }
              if (am2) {
                i = am2-L1[0];
                c = i + ioffset(i) - L1[1];
              }
              T3 = a*ioffset(am2+1)+c;

              L1[xyz] = L1[xyz] + 2;
              am2 = am2 + 2;
            }

            term = abc[1][xyz]*I0[T2] - abc[2][xyz]*I1[T2];
            T2++;
            if (L1[xyz] > 1) {
              term += oo2gam*(L1[xyz]-1)*(I2[T3] - I3[T3]);
              T3++;
            }

            term *= charge[pt];
            *via1a2mi += term;
            via1a2mi++;
          }
        }
      }
    }
  }
}
  
void
GaussianPointChargeIntJF::calc_f(double F[], int n, double t)
{
  int i, m, k;
  int m2;
  double t2;
  double num;
  double sum;
  double term1, term2;
  const double K = 1.0/M_2_SQRTPI;
  double et;


  if (t > 20.0) {
    t2 = 2*t;
    et = exp(-t);
    t = sqrt(t);
    F[0] = K*erf(t)/t;
    for (m=0; m <= n-1; m++) {
      F[m+1] = ((2*m + 1)*F[m] - et)/t2;
    }
  } else {
    et = exp(-t);
    t2 = 2*t;
    m2 = 2*n;
    num = df[m2];
    i=0;
    sum = 1.0/(m2+1);
    do {
      i++;
      num = num*t2;
      term1 = num/df[m2+2*i+2];
      sum += term1;
    } while (fabs(term1) > 1.0e-17 && i < maxfact);

    F[n] = sum*et;
    for (m=n-1;m>=0;m--) {
      F[m] = (t2*F[m+1] + et)/(2*m+1);
    }
  }
}

double *
GaussianPointChargeIntJF::do_pe_doublet(double abc[3][3], double F[20],
                                        int am1, int am2, int m, int si,
                                        double oo2gam)
{
  int am;
  double *I0, *I1, *I2, *I3, *I4, *I5;
  int L0[3], L1[3];
  int ii, kk, jj, ll;
  int xyz, c12, c13, c14;
  int t20, t30, t40;
  int i, j;
  int a, c;
  int T2max;
  int T1, T2, T3, T4;
  double t;
  double *via1a2m = vint[am1][am2][m];
  double *via1a2mi = via1a2m;
  
  am = am1+am2;
  if (vint_done[am1][am2][m])
    return vint[am1][am2][m];

  if (am1==0 && am2==0) {
    *via1a2m = F[m];
    return via1a2m;

  } else {
    if (am1) {
      I0 = do_pe_doublet(abc, F, am1-1, am2, m, si, oo2gam);
      I1 = do_pe_doublet(abc, F, am1-1, am2, m+1, si, oo2gam);

      if (am1 > 1) {
        I2 = do_pe_doublet(abc, F, am1-2, am2, m, si, oo2gam);
        I3 = do_pe_doublet(abc, F, am1-2, am2, m+1, si, oo2gam);
      }

      if (am2) {
        I4 = do_pe_doublet(abc, F, am1-1, am2-1, m, si, oo2gam);
        I5 = do_pe_doublet(abc, F, am1-1, am2-1, m+1, si, oo2gam);
      }

      // create all am components of si
      T2max = (am1*(am1+1)*(am2+1)*(am2+2))/4;
      T1 = T2 = T3 = T4 = 0;
      for (ii = 0; ii <= am1; ii++) {
        L0[0] = am1 - ii;

        for (jj = 0; jj <= ii; jj++) {
          L0[1] = ii - jj;
          L0[2] = jj ;

          if(L0[0])
            xyz = 0;
          else if(L0[1])
            xyz = 1;
          else if(L0[2])
            xyz = 2;

          // create all am components of sj
          for (kk = 0; kk <= am2; kk++) {
            L1[0] = am2 - kk;

            for (ll = 0; ll <= kk; ll++) {
              L1[1] = kk - ll;
              L1[2] = ll ;

              if (T2==T2max) {
                L0[xyz] = L0[xyz] - 1;
                am1 = am1 - 1;
                a = c = 0;
                if (am1) {
                  i = am1-L0[0];
                  a = i + ioffset(i) - L0[1];
                }
                if (am2) {
                  i = am2-L1[0];
                  c = i + ioffset(i) - L1[1];
                }
                T2 = a*ioffset(am2+1)+c;

                L0[xyz] = L0[xyz] - 1;
                am1 = am1 - 1;
                a = c = 0;
                if (am1){
                  i = am1-L0[0];
                  a = i + ioffset(i) - L0[1];
                }
                if (am2) {
                  i = am2-L1[0];
                  c = i + ioffset(i) - L1[1];
                }
                T3 = a*ioffset(am2+1)+c;

                L0[xyz] = L0[xyz] + 1;
                am1 = am1 + 1;
                L1[0] = L1[0] - 1;
                am2 = am2 - 1;
                a = c = 0;
                if (am1) {
                  i = am1-L0[0];
                  a = i + ioffset(i) - L0[1];
                }
                if (am2) {
                  i = am2-L1[0];
                  c = i + ioffset(i) - L1[1];
                }
                T4 = a*ioffset(am2+1)+c;

                L1[0] = L1[0] + 1;
                am2 = am2 + 1;
                L0[xyz] = L0[xyz] + 1;
                am1 = am1 + 1;
              }

              *via1a2mi = abc[0][xyz]*I0[T2] - abc[2][xyz]*I1[T2];
              T2++;
              if (L0[xyz]>1) {
                *via1a2mi += oo2gam*(L0[xyz]-1)*(I2[T3] - I3[T3]);
                T3++;
              }
              if (L1[xyz]) {
                *via1a2mi += oo2gam*(L1[xyz])*(I4[T4] - I5[T4]);
                T4++;
              }
              via1a2mi++;
            }
          }
        }
      }
      return via1a2m;

    } else if (am2) {
      I0 = do_pe_doublet(abc, F, am1, am2-1, m, si, oo2gam);
      I1 = do_pe_doublet(abc, F, am1, am2-1, m+1, si, oo2gam);

      if (am2 > 1) {
        I2 = do_pe_doublet(abc, F, am1, am2-2, m, si, oo2gam);
        I3 = do_pe_doublet(abc, F, am1, am2-2, m+1, si, oo2gam);
      }

      T2max = am2*(am2+1)/2;
      T1 = T2 = T3 = T4 = 0;

      // create all am components of si - all s since am1==0
      L0[0] = 0;
      L0[1] = 0;
      L0[2] = 0;

      // create all am components of sj
      for (kk = 0; kk <= am2; kk++) {
        L1[0] = am2 - kk;

        for (ll = 0; ll <= kk; ll++) {
          L1[1] = kk - ll;
          L1[2] = ll ;

          if (L1[0])
            xyz = 0;
          else if (L1[1])
            xyz = 1;
          else if (L1[2])
            xyz = 2;

          if (T2==T2max) {
            L1[xyz] = L1[xyz] - 1;
            am2 = am2 - 1;
            a = c = 0;
            if (am1) {
              i = am1-L0[0];
              a = i + ioffset(i) - L0[1];
            }
            if (am2) {
              i = am2-L1[0];
              c = i + ioffset(i) - L1[1];
            }
            T2 = a*ioffset(am2+1)+c;

            L1[xyz] = L1[xyz] - 1;
            am2 = am2 - 1;
            a = c = 0;
            if (am1) {
              i = am1-L0[0];
              a = i + ioffset(i) - L0[1];
            }
            if (am2) {
              i = am2-L1[0];
              c = i + ioffset(i) - L1[1];
            }
            T3 = a*ioffset(am2+1)+c;
            L1[xyz] = L1[xyz] + 2;
            am2 = am2 + 2;
          }

          *via1a2mi = abc[1][xyz]*I0[T2] - abc[2][xyz]*I1[T2];
          T2++;
          if (L1[xyz]>1) {
            *via1a2mi += oo2gam*(L1[xyz]-1)*(I2[T3] - I3[T3]);
            T3++;
          }
          via1a2mi++;
        }
      }

      vint_done[am1][am2][m] = 1;
      return via1a2m;
    }
  }
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

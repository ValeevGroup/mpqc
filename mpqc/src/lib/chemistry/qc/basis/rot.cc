
#ifdef __GNUC__
#pragma implementation
#endif

#include <chemistry/qc/intv2/transform.h>
#include <chemistry/qc/basis/rot.h>

Rotation::Rotation(int a, SymmetryOperation& so, int pure = 0):
  _am(0),
  _n(0),
  r(0)
{
  if (pure) init_pure(a,so);
  else init(a,so);
}

// Compute the transformation matrices for general pure am
// by summing contributions from the cartesian components
// using the P (xyz) transformation matrix.  This is done as a
// matrix outer product, keeping only the unique terms.
void
Rotation::init_pure(int a, SymmetryOperation&so)
{
  done();

  _am=a;
  
  SphericalTransformIter I(_am);
  SphericalTransformIter J(_am, 1);
  RedundantCartesianSubIter K(_am);
  int lI[3];
  int m, iI;
  
  _n = I.n();

  r = new double*[_n];
  for (m=0; m<_n; m++) {
      r[m] = new double[_n];
      memset(r[m],0,sizeof(double)*_n);
    }

  for (I.start(); I; I.next()) {
      for (J.start(); J; J.next()) {
          double coef = I.coef()*J.coef();
          double tmp = 0.0;
#if 0
          printf("T(%d,%d) += %6.4f", I.bfn(), J.bfn(), coef);
          double tmp2 = coef;
          int iJ, lJ[3];
          for (m=0; m < 3; m++) {
              lI[m] = I.l(m);
              lJ[m] = J.l(m);
            }
          
          for (m=0; m < _am; m++) {
              for (iI=0; lI[iI]==0; iI++);
              lI[iI]--;
              for (iJ=0; lJ[iJ]==0; iJ++);
              lJ[iJ]--;
              tmp2 *= so(iI,iJ);
              printf(" * so(%d,%d) [=%4.2f]",
                     iI,iJ,so(iI,iJ));
            }
          printf(" = %8.6f\n", tmp2);
          tmp += tmp2;
#else          
          for (K.start(J.a(), J.b(), J.c()); K; K.next()) {
              //printf("T(%d,%d) += %6.4f", I.bfn(), J.bfn(), coef);
              double tmp2 = coef;
              for (m=0; m < 3; m++) {
                  lI[m] = I.l(m);
                }
      
              for (m=0; m < _am; m++) {
                  for (iI=0; lI[iI]==0; iI++);
                  lI[iI]--;
                  tmp2 *= so(iI,K.axis(m));
                  //printf(" * so(%d,%d) [=%4.2f]",
                  //       iI,K.axis(m),so(iI,K.axis(m)));
                }
              //printf(" = %8.6f\n", tmp2);
              tmp += tmp2;
            }
#endif
          r[I.bfn()][J.bfn()] += tmp;
        }
    }
}

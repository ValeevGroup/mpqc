
#ifdef __GNUC__
#pragma implementation
#endif

#include <chemistry/qc/intv2/transform.h>
#include <chemistry/qc/basis/rot.h>

void
Rotation::done() {
  if (r) {
    for (int i=0; i < _n; i++) {
      if (r[i]) delete[] r[i];
    }
    delete[] r;
    r=0;
  }
  _n=0;
}

Rotation::Rotation(int n) :
  _am(0),
  _n(n),
  r(0)
{
  if (_n) {
    r = new double*[_n];
    for (int i=0; i < _n; i++)
      r[i] = new double[_n];
  }
}

Rotation::Rotation(const Rotation& rot) :
  _am(0),
  _n(0),
  r(0)
{
  *this = rot;
}

Rotation::Rotation(int a, SymmetryOperation& so, int pure):
  _am(0),
  _n(0),
  r(0)
{
  if (pure)
    init_pure(a,so);
  else
    init(a,so);
}

Rotation::~Rotation()
{
  done();
}

Rotation&
Rotation::operator=(const Rotation& rot)
{
  done();

  _n = rot._n;
  _am = rot._am;

  if (_n && rot.r) {
    r = new double*[_n];
    for (int i=0; i < _n; i++) {
      r[i] = new double[_n];
      memcpy(r[i],rot.r[i],sizeof(double)*_n);
    }
  }

  return *this;
}

// Compute the transformation matrices for general cartesian shells
// using the P (xyz) transformation matrix.  This is done as a
// matrix outer product, keeping only the unique terms.
// Written by clj...blame him
void
Rotation::init(int a, SymmetryOperation&so)
{
  done();

  _am=a;
  
  CartesianIter I(_am);
  RedundantCartesianIter J(_am);
  int lI[3];
  int k, iI;
  
  _n = I.n();
  r = new double*[_n];

  for (I.start(); I; I.next()) {
    r[I.bfn()] = new double[_n];
    memset(r[I.bfn()],0,sizeof(double)*_n);

    for (J.start(); J; J.next()) {
      double tmp = 1.0;

      for (k=0; k < 3; k++) {
        lI[k] = I.l(k);
      }
      
      for (k=0; k < _am; k++) {
        for (iI=0; lI[iI]==0; iI++);
        lI[iI]--;
        //tmp *= so(iI,J.axis(k));
        tmp *= so(J.axis(k),iI);
      }

      r[I.bfn()][J.bfn()] += tmp;
    }
  }
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
                  //tmp2 *= so(iI,K.axis(m));
                  tmp2 *= so(K.axis(m),iI);
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

// returns the result of rot*this
Rotation
Rotation::operate(const Rotation& rot) const
{
  if (_n != rot._n) {
    fprintf(stderr,"Rotation::operate(): dimensions don't match\n");
    fprintf(stderr,"  %d != %d\n",rot._n,_n);
    abort();
  }
  
  Rotation ret(_n);
  ret._am = _am;
  
  for (int i=0; i < _n; i++) {
    for (int j=0; j < _n; j++) {
      double t=0;
      for (int k=0; k < _n; k++)
        t += rot.r[i][k] * r[k][j];
      ret.r[i][j] = t;
    }
  }

  return ret;
}

Rotation
Rotation::sim_transform(const Rotation& rot) const
{
  int i,j,k;

  if (rot._n != _n) {
    fprintf(stderr,"Rotation::sim_transform(): dimensions don't match\n");
    fprintf(stderr,"  %d != %d\n",rot._n,_n);
    abort();
  }
  
  Rotation ret(_n), foo(_n);
  ret._am = foo._am = _am;

  // foo = r * d
  for (i=0; i < _n; i++) {
    for (j=0; j < _n; j++) {
      double t=0;
      for (k=0; k < _n; k++)
        t += rot.r[i][k] * r[k][j];
      foo.r[i][j] = t;
    }
  }

  // ret = (r*d)*r~ = foo*r~
  for (i=0; i < _n; i++) {
    for (j=0; j < _n; j++) {
      double t=0;
      for (k=0; k < _n; k++)
        t += foo.r[i][k]*rot.r[j][k];
      ret.r[i][j]=t;
    }
  }

  return ret;
}
    
double
Rotation::trace() const {
  double t=0;
  for (int i=0; i < _n; i++)
    t += r[i][i];
  return t;
}

void
Rotation::print() const
{
  for (int i=0; i < _n; i++) {
    printf("%5d ",i+1);
    for (int j=0; j < _n; j++) {
      printf(" %10.7f",r[i][j]);
    }
    printf("\n");
  }
}

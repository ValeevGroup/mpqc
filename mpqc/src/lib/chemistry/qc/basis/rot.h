
#ifndef _chemistry_qc_basis_rot_h
#define _chemistry_qc_basis_rot_h

//#ifdef __GNUC__
//#pragma interface
//#endif

#include <math/symmetry/pointgrp.h>
#include <chemistry/qc/basis/gaussshell.h>

class Rotation {
  private:
    int _n;
    int _am;
    double **r;
    
    void done() {
      if (r) {
        for (int i=0; i < _n; i++) {
          if (r[i]) delete[] r[i];
        }
        delete[] r;
        r=0;
      }
      _n=0;
    }

  public:
    inline void init(int a, SymmetryOperation&so);
    
    Rotation(int a, SymmetryOperation& so) : _am(0), _n(0), r(0) { init(a,so);}
    ~Rotation() { done(); }

    int am() const { return _am; }
    int dim() const { return _n; }
    
    double& operator()(int i, int j) { return r[i][j]; }
    double* operator[](int i) { return r[i]; }
    
    void print() const;
    
    inline double trace() const {
      double t=0;
      for (int i=0; i < _n; i++)
        t += r[i][i];
      return t;
    }
};


// Compute the transformation matrices for general cartesian shells
// using the P (xyz) transformation matrix.  This is done as a
// matrix outer product, keeping only the unique terms.
// Written by clj...blame him
inline void
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
        tmp *= so(iI,J.axis(k));
      }

      r[I.bfn()][J.bfn()] += tmp;
    }
  }
}

inline void
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

#endif


#ifndef _chemistry_qc_basis_rot_h
#define _chemistry_qc_basis_rot_h

#ifdef __GNUC__
#pragma interface
#endif

#include <math/symmetry/pointgrp.h>
#include <chemistry/qc/basis/gaussshell.h>

class Rotation {
  private:
    int _n;
    int _am;
    double **r;
    
    void done();

  public:
    void init(int a, SymmetryOperation&so);
    void init_pure(int a, SymmetryOperation&so);
    
    Rotation(int a, SymmetryOperation& so, int pure = 0);
    ~Rotation();

    int am() const { return _am; }
    int dim() const { return _n; }
    
    double& operator()(int i, int j) { return r[i][j]; }
    double* operator[](int i) { return r[i]; }
    
    void print() const;
    
    double trace() const;
};


#endif

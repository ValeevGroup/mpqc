
#ifndef _chemistry_qc_basis_shellrot_h
#define _chemistry_qc_basis_shellrot_h

#ifdef __GNUC__
#pragma interface
#endif

#include <math/symmetry/pointgrp.h>

class Integral;
class ShellRotation {
  private:
    int n_;
    int am_;
    double **r;
    
    void done();

  public:
    void init(int a, SymmetryOperation&, const RefIntegral&);
    void init_pure(int a, SymmetryOperation&, const RefIntegral&);
    
    ShellRotation(int n);
    ShellRotation(const ShellRotation&);
    ShellRotation(int a, SymmetryOperation&, const RefIntegral&, int pure =0);
    virtual ~ShellRotation();

    ShellRotation& operator=(const ShellRotation&);
    
    int am() const { return am_; }
    int dim() const { return n_; }
    
    double& operator()(int i, int j) { return r[i][j]; }
    double* operator[](int i) { return r[i]; }
    
    ShellRotation operate(const ShellRotation&) const;
    ShellRotation sim_transform(const ShellRotation&) const;
    
    double trace() const;
    
    void print() const;
};

#endif

// Local Variables:
// mode: c++
// eval: (c-set-style "ETS")


#ifndef _chemistry_qc_basis_gaussshell_h
#define _chemistry_qc_basis_gaussshell_h

#include <stdio.h>
#include <math/topology/point.h>

class CartesianIter
{
 private:
  int a_;
  int c_;
  int l_;
  int bfn_;
 public:
  inline CartesianIter(int l):l_(l) {}
  inline ~CartesianIter() {}
  inline void start() { bfn_=a_=c_=0; }
  inline void next() { if (c_<l_-a_) c_++; else {c_=0; a_++;} bfn_++; }
  inline operator int() { return a_<=l_; }
  inline int a() { return a_; }
  inline int b() { return l_-a_-c_; }
  inline int c() { return c_; }
  inline int l() { return l_; }
  inline int bfn() { return bfn_; }
};

class KeyVal;
class GaussianShell
{
 public:
  enum PrimitiveType { Normalized, Unnormalized };
  enum GaussianType { Cartesian, Pure };
 private:
  int nprim;
  int ncon;
  int nfunc;
  int* l;
  int* puream;
  double* exp;
  double** coef;  // contraction coefficients for unnormalized primitives
  double** norm;
  double shell_normalization(int);
  void convert_coef();
  void normalize_shell();
  void compute_nfunc();
  PrimitiveType keyval_init(KeyVal&,int,int);
  static const char* amtypes;
  static const char* AMTYPES;
 public:
  // Users of GaussianShell must pass pointers to newed memory that is kept
  // by GaussianShell and deleted by the destructor.
  // The arguments for the following ctor are:
  //   ncn is the number of contracted functions (1 except for SP and gen. con.)
  //   nprm is the number of primitives
  //   e gives the exponents (length nprm)
  //   am gives the angular momentum (length ncn)
  //   pure is 1 for pure am and 0 for cartesian (length ncn)
  //   c are the contraction coefficients (length ncn by nprm)
  //   pt describes whether the primitive functions are to be considered
  //     normalized or unnormalized.  this effects whether or not c is
  //     manipulated to give the correct normalization.
  GaussianShell(
    int ncn,
    int nprm,
    double* e,
    int* am,
    int* pure,
    double** c,
    PrimitiveType pt = GaussianShell::Normalized);
  // In this ctor pure is either GaussianShell::Cartesian or
  // Gaussian::Pure and all of the contracted functions are treated in that way.
  // (The user doesn\'t need to compute generate a int*pure vector in this case.)
  GaussianShell(
    int ncn,
    int nprm,
    double* e,
    int* am,
    GaussianType pure,
    double** c,
    PrimitiveType pt = GaussianShell::Normalized);
  GaussianShell(KeyVal&);
  GaussianShell(KeyVal&,int pure);
  ~GaussianShell();
  inline int nprimitive() { return nprim; }
  inline int ncontraction() { return ncon; }
  inline int nfunction() { return nfunc; }
  inline int am(int con) { return l[con]; }
  int max_am();
  inline char amchar(int con) { return amtypes[l[con]]; }
  int nfunction(int con);
  inline int is_cartesian(int con) { return !puream[con]; }
  inline int is_pure(int con) { return puream[con]; }
  // returns the con coef for unnormalized primitives
  inline double coefficient_unnorm(int con,int prim) { return coef[con][prim]; }
  // returns the con coef for normalized primitives
  double coefficient_norm(int con,int prim);
  inline double normalization(int con,int bfn) { return norm[con][bfn]; }
  inline double exponent(int iprim) { return exp[iprim]; }

  // compute the value of this shell at offset r
  int values(cart_point& r, double* basis_values);
  int grad_values(cart_point& r, double* g_values, double* basis_values=0);

  void print(FILE*fp=stdout);
};
  
#endif

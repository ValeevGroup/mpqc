
/* nihmatrix.h -- declarations for the vector and matrix classes
 *
 *      THIS SOFTWARE FITS THE DESCRIPTION IN THE U.S. COPYRIGHT ACT OF A
 *      "UNITED STATES GOVERNMENT WORK".  IT WAS WRITTEN AS A PART OF THE
 *      AUTHOR'S OFFICIAL DUTIES AS A GOVERNMENT EMPLOYEE.  THIS MEANS IT
 *      CANNOT BE COPYRIGHTED.  THIS SOFTWARE IS FREELY AVAILABLE TO THE
 *      PUBLIC FOR USE WITHOUT A COPYRIGHT NOTICE, AND THERE ARE NO
 *      RESTRICTIONS ON ITS USE, NOW OR SUBSEQUENTLY.
 *
 *  Author:
 *      E. T. Seidl
 *      Bldg. 12A, Rm. 2033
 *      Computer Systems Laboratory
 *      Division of Computer Research and Technology
 *      National Institutes of Health
 *      Bethesda, Maryland 20892
 *      Internet: seidl@alw.nih.gov
 *      March, 1993
 */

#ifndef _mathQC_h
#define _mathQC_h

#include <math.h>
#include <iostream.h>

#include <util/unix/cct_cprot.h>
#include <util/class/class.h>
#include <util/state/state.h>
#include <util/container/array.h>

class RefPoint;

#define V_BASE virtual public SavableState

class DMatrix;
class KeyVal;

class DVector : V_BASE {
#   define CLASSNAME DVector
#   define HAVE_CTOR
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  protected:
    int n;
    double *d;

  public:
    DVector();
    DVector(int sz);
    DVector(int sz, double *dp);
    DVector(const DVector& da);
    DVector(RefPoint&);
    DVector(StateIn&);
    DVector(KeyVal&);
    ~DVector();

    void resize(int);
    inline void init() { if(d) delete[] d; d=0; n=0; }
    inline void zero() { for (int i=0; i < n; i++) d[i]=0; }

    DVector& operator=(const DVector&);
    DVector& operator=(RefPoint&);

    DVector& operator+=(const double);
    DVector& operator*=(const double);

    DVector& operator+=(const DVector&);
    DVector& operator-=(const DVector&);
    DVector& operator*=(const DVector&);

    void save_data_state(StateOut& so);
    void restore_data_state(int, StateIn& si);

    inline int dim() const { return n; }
    inline double& operator[](int i) { return d[i]; }
    const double& operator[](int i) const { return d[i]; }
    double operator()(int i) const { return d[i]; }
    double * pointer() { return d; }
    const double * pointer() const { return d; }

    const double maxval() const;

    double dot(const DVector&) const; // dot or inner product
    double norm() const; // the 2 norm
    void normalize();
    DMatrix ccross(const DVector&) const; // direct or outer product
    DVector cross(const DVector&) const;

    void print(const char* =0,ostream& =cout, int =10) const;
    void print(ostream&) const;
  };
ARRAY_dec(DVector);
DescribedClass_REF_dec(DVector);

DVector operator+(const DVector&, const double);
DVector operator+(const double, const DVector&);

DVector operator*(const DVector&, const double);
DVector operator*(const double, const DVector&);

DVector operator+(const DVector&, const DVector&);
DVector operator-(const DVector&, const DVector&);
DVector operator*(const DVector&, const DVector&);

DVector operator*(DVector&, DMatrix&);
DVector operator*(DMatrix&, DVector&);

DMatrix ccross(const DVector&, const DVector&);
double dot(const DVector&, const DVector&);

//////////////////////////////////////////////////////////////////////////

class DMatrix : V_BASE {
#   define CLASSNAME DMatrix
#   define HAVE_CTOR
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  protected:
    int n1;
    int n2;
    double *d;
    double **dp;

  public:
    DMatrix();
    DMatrix(int, int);
    DMatrix(int, int, double*);
    DMatrix(const DMatrix& da);
    DMatrix(KeyVal&);
    DMatrix(StateIn&);
    ~DMatrix();

    void init();
    void resize(int,int);
    inline void zero() { for (int i=0; i < n1*n2; i++) d[i]=0; }

    DMatrix& operator=(const DMatrix&);

    DMatrix& operator+=(const double);
    DMatrix& operator*=(const double);

    DMatrix& operator+=(const DMatrix&);
    DMatrix& operator-=(const DMatrix&);
    DMatrix& operator*=(const DMatrix&);

  // functions which return a modified this
    DMatrix transpose() const;
    DMatrix inverse() const;

  // functions which alter this
    void transpose_IP();
    double invert();
    
  // other functions
    inline double trace() const {
      if(n1!=n2) return 0; double t=0;
      for(int i=0; i < n1; i++) t += dp[i][i];
      return t;
      }

    double solve_lin(DVector&);
    DMatrix eigenvectors();
    DVector eigenvalues();
    void diagonalize(DVector&,DMatrix&,double =1.0e-15);

    void save_data_state(StateOut&);
    void restore_data_state(int,StateIn&);

    inline double* operator[](int i) { return dp[i]; }
    inline double& operator()(int i,int j) { return dp[i][j]; }

    inline const int nrow() const { return n1; }
    inline const int ncol() const { return n2; }

    const double maxval() const;

    void print(const char* =0,ostream& =cout, int =10) const;
    void print(ostream& os) const;
  };
DescribedClass_REF_dec(DMatrix);

DMatrix operator+(DMatrix&,DMatrix&);
DMatrix operator-(DMatrix&,DMatrix&);
DMatrix operator*(DMatrix&,DMatrix&);

DMatrix operator*(DMatrix&,double);
DMatrix operator*(double,DMatrix&);

//////////////////////////////////////////////////////////////////////////

#undef V_BASE
#endif /* _mathQC_h */

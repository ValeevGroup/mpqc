
#ifndef _chemistry_qc_wfn_obwfn_h
#define _chemistry_qc_wfn_obwfn_h

#ifdef __GNUC__
#pragma interface
#endif

#include <chemistry/qc/wfn/wfn.h>

class OneBodyWavefunction: public Wavefunction
{
#   define CLASSNAME OneBodyWavefunction
#   include <util/state/stated.h>
#   include <util/class/classda.h>
 protected:
    ResultRefSymmSCMatrix _density;
 public:
    OneBodyWavefunction(const RefKeyVal&);
    OneBodyWavefunction(StateIn&);
    ~OneBodyWavefunction();
    void save_data_state(StateOut&);

    virtual RefSCMatrix eigenvectors() = 0;
    virtual double occupation(int vectornum) = 0;
    double orbital(const SCVector3& r, int iorb);
    double orbital_density(const SCVector3& r, int iorb, double* orbval = 0);

    double density(const SCVector3&);
    RefSymmSCMatrix density();

    void print(SCostream&o=SCostream::cout);
};

// This is useful as an initial guess for other one body wavefunctions
class HCoreWfn: public OneBodyWavefunction {
#   define CLASSNAME HCoreWfn
#   define HAVE_KEYVAL_CTOR
//#   include <util/state/stated.h>
#   include <util/class/classd.h>
  public:
    HCoreWfn(const RefKeyVal&);
    ~HCoreWfn();
    RefSCMatrix eigenvectors();
    double occupation(int vectornum);
    void compute();
};

#endif

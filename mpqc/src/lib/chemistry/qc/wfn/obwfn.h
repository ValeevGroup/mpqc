
#ifndef _chemistry_qc_wfn_obwfn_h
#define _chemistry_qc_wfn_obwfn_h

#ifdef __GNUC__
#pragma interface
#endif

#include <chemistry/qc/wfn/wfn.h>
#include <chemistry/qc/wfn/hcore.h>

SavableState_REF_fwddec(OneBodyWavefunction);
class OneBodyWavefunction: public Wavefunction {
#   define CLASSNAME OneBodyWavefunction
#   include <util/state/stated.h>
#   include <util/class/classda.h>
 protected:
    ResultRefSymmSCMatrix density_;
    AccResultRefSCMatrix eigenvectors_;

 public:
    OneBodyWavefunction(StateIn&);
    OneBodyWavefunction(const RefKeyVal&);
    ~OneBodyWavefunction();

    void save_data_state(StateOut&);

    virtual RefSCMatrix eigenvectors() = 0;
    virtual double occupation(int irrep, int vectornum) = 0;

    virtual RefSCMatrix projected_eigenvectors(const RefOneBodyWavefunction&);
    virtual RefSCMatrix hcore_guess();

    double orbital(const SCVector3& r, int iorb);
    double orbital_density(const SCVector3& r, int iorb, double* orbval = 0);

    double density(const SCVector3&);
    RefSymmSCMatrix density();

    void print(ostream&o=cout);
};
SavableState_REF_dec(OneBodyWavefunction);

// This is useful as an initial guess for other one body wavefunctions
class HCoreWfn: public OneBodyWavefunction {
#   define CLASSNAME HCoreWfn
#   define HAVE_STATEIN_CTOR
#   define HAVE_KEYVAL_CTOR
#   include <util/class/classd.h>
  private:
    RefAccumHCore accumh;

    int nirrep_;
    int *docc;
    int *socc;
    
    void compute();

  public:
    HCoreWfn(StateIn&);
    HCoreWfn(const RefKeyVal&);
    ~HCoreWfn();

    void save_data_state(StateOut&);

    double occupation(int irrep, int vectornum);

    RefSCMatrix eigenvectors();
};
SavableState_REF_dec(HCoreWfn);

#endif

// Local Variables:
// mode: c++
// eval: (c-set-style "ETS")

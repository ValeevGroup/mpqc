
#ifndef _chemistry_qc_wfn_obwfn_h
#define _chemistry_qc_wfn_obwfn_h

#ifdef __GNUC__
#pragma interface
#endif

#include <chemistry/qc/wfn/wfn.h>
#include <chemistry/qc/wfn/hcore.h>

SavableState_REF_fwddec(OneBodyWavefunction);
class OneBodyWavefunction: public Wavefunction
{
#   define CLASSNAME OneBodyWavefunction
#   include <util/state/stated.h>
#   include <util/class/classda.h>
 protected:
    ResultRefSymmSCMatrix _density;
    AccResultRefSCMatrix _eigenvectors;

    virtual void form_density(const RefSCMatrix& vec,
                              const RefSymmSCMatrix& density,
                              const RefSymmSCMatrix& density_diff,
                              const RefSymmSCMatrix& open_density,
                              const RefSymmSCMatrix& open_density_diff);

 public:
    OneBodyWavefunction(const OneBodyWavefunction&);
    OneBodyWavefunction(const RefKeyVal&);
    OneBodyWavefunction(StateIn&);
    ~OneBodyWavefunction();

    OneBodyWavefunction & operator=(const OneBodyWavefunction&);
    
    void save_data_state(StateOut&);

    virtual RefSCMatrix eigenvectors() = 0;
    virtual double occupation(int vectornum) = 0;

    virtual RefSCMatrix projected_eigenvectors(const RefOneBodyWavefunction&);

    // returns a matrix which transforms AO's to orthogonal AO's
    // can be overridden, but defaults to S^-1/2
    virtual RefSymmSCMatrix ao_to_orthog_ao();

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
    int ndocc;
    int nsocc;
    
    RefAccumHCore _accumh;
    
    void compute();

  public:
    HCoreWfn(const OneBodyWavefunction&);
    HCoreWfn(const RefKeyVal&);
    HCoreWfn(const HCoreWfn&);
    HCoreWfn(StateIn&);
    ~HCoreWfn();

    HCoreWfn & operator=(const HCoreWfn&);
    
    void save_data_state(StateOut&);

    double occupation(int vectornum);

    RefSCMatrix eigenvectors();
};
SavableState_REF_dec(HCoreWfn);

#endif

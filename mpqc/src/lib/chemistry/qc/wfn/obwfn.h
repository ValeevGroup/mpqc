
#ifndef _chemistry_qc_wfn_obwfn_h
#define _chemistry_qc_wfn_obwfn_h

#include <chemistry/qc/wfn/wfn.h>

class OneBodyWavefunction: public Wavefunction
{
#   define CLASSNAME OneBodyWavefunction
#   include <util/state/stated.h>
#   include <util/class/classda.h>
 private:
    ResultRefSymmSCMatrix _density;
 public:
    OneBodyWavefunction(KeyVal&);
    OneBodyWavefunction(StateIn&);
    ~OneBodyWavefunction();
    void save_data_state(StateOut&);

    virtual RefSCMatrix eigenvectors() = 0;
    virtual double occupation(int vectornum) = 0;
    double orbital(cart_point& r, int iorb);
    double orbital_density(cart_point& r, int iorb, double* orbval = 0);

    double density(cart_point&);
    RefSymmSCMatrix density();

    void print(SCostream&o=SCostream::cout);
};

#endif


#ifndef _chemistry_qc_psi_psici_h
#define _chemistry_qc_psi_psici_h

#include <chemistry/qc/wfn/wfn.h>
#include "file11.h"
#include "psiinput.h"

class PSI_CI: public Wavefunction
{
#define CLASSNAME PSI_CI
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#include <util/state/stated.h>
#include <util/class/classd.h>

 private:
    PSI_Input_CI psi_in;
 protected:
    void compute();
 public:
    PSI_CI(const RefKeyVal&);
    PSI_CI(StateIn&);
    virtual ~PSI_CI();
    void save_data_state(StateOut&);

    //double orbital(cart_point& r, int iorb);
    //double orbital_density(cart_point& r, int iorb, double* orbval = 0);

    double density(cart_point&);
    RefSymmSCMatrix density();

    void print(ostream&o=cout);

    int gradient_implemented();
    int value_implemented();
};

#endif

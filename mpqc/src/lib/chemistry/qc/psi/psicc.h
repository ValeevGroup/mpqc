
#ifndef _chemistry_qc_psi_psicc_h
#define _chemistry_qc_psi_psicc_h

#include <chemistry/qc/wfn/wfn.h>
#include "file11.h"
#include "psiinput.h"

class PSI_CCSD_T: public Wavefunction
{
#define CLASSNAME PSI_CCSD_T
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#include <util/state/stated.h>
#include <util/class/classd.h>

 private:
    PSI_Input_CC psi_in;
 protected:
    void compute();
 public:
    PSI_CCSD_T(const RefKeyVal&);
    PSI_CCSD_T(StateIn&);
    virtual ~PSI_CCSD_T();
    void save_data_state(StateOut&);

    //double orbital(cart_point& r, int iorb);
    //double orbital_density(cart_point& r, int iorb, double* orbval = 0);

    double density(cart_point&);
    RefSymmSCMatrix density();

    void print(ostream&o=cout);

    int gradient_implemented();
    int value_implemented();
};

class PSI_CCSD: public Wavefunction
{
#define CLASSNAME PSI_CCSD
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#include <util/state/stated.h>
#include <util/class/classd.h>

 private:
    PSI_Input_CC psi_in;
 protected:
    void compute();
 public:
    PSI_CCSD(const RefKeyVal&);
    PSI_CCSD(StateIn&);
    virtual ~PSI_CCSD();
    void save_data_state(StateOut&);

    //double orbital(cart_point& r, int iorb);
    //double orbital_density(cart_point& r, int iorb, double* orbval = 0);

    double density(cart_point&);
    RefSymmSCMatrix density();

    void print(ostream&o=cout);

    int gradient_implemented();
    int value_implemented();
};

class PSI_CCSDT: public Wavefunction
{
#define CLASSNAME PSI_CCSDT
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#include <util/state/stated.h>
#include <util/class/classd.h>

 private:
    PSI_Input_CC psi_in;
 protected:
    void compute();
 public:
    PSI_CCSDT(const RefKeyVal&);
    PSI_CCSDT(StateIn&);
    virtual ~PSI_CCSDT();
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

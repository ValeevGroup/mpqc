
#ifdef __GNUC__
#pragma interface
#endif

#ifndef _chemistry_qc_psi_psi_h
#define _chemistry_qc_psi_psi_h

#include <chemistry/qc/wfn/wfn.h>
#include <chemistry/qc/psi/psiinput.h>

class PsiWfn: public Wavefunction {
#   define CLASSNAME PsiWfn
#   include <util/state/stated.h>
#   include <util/class/classda.h>
  protected:
    PSI_Input psi_in;
    void compute();

    virtual void write_input(int conv) = 0;
    virtual double read_energy() = 0;
    void write_basic_input(int conv, const char *wfn);
  public:
    PsiWfn(const RefKeyVal&);
    PsiWfn(StateIn&);
    virtual ~PsiWfn();
    void save_data_state(StateOut&);

    double density(const SCVector3&);
    RefSymmSCMatrix density();

    void print(ostream&o=cout);

    int spin_polarized();
    int nelectron();

    int gradient_implemented() const;
    int value_implemented() const;
};

class PsiHF: public PsiWfn {
#   define CLASSNAME PsiHF
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  protected:
    void write_input(int conv);
    double read_energy();
  public:
    PsiHF(const RefKeyVal&);
    PsiHF(StateIn&);
    ~PsiHF();
    void save_data_state(StateOut&);
};

class PsiCCSD: public PsiWfn {
#   define CLASSNAME PsiCCSD
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  protected:
    void write_input(int conv);
    double read_energy();
  public:
    PsiCCSD(const RefKeyVal&);
    PsiCCSD(StateIn&);
    ~PsiCCSD();
    void save_data_state(StateOut&);
};

class PsiCCSD_T: public PsiWfn {
#   define CLASSNAME PsiCCSD_T
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  protected:
    void write_input(int conv);
    double read_energy();
  public:
    PsiCCSD_T(const RefKeyVal&);
    PsiCCSD_T(StateIn&);
    ~PsiCCSD_T();
    void save_data_state(StateOut&);
};

class PsiCCSDT: public PsiWfn {
#   define CLASSNAME PsiCCSDT
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  protected:
    void write_input(int conv);
    double read_energy();
  public:
    PsiCCSDT(const RefKeyVal&);
    PsiCCSDT(StateIn&);
    ~PsiCCSDT();
    void save_data_state(StateOut&);
    int gradient_implemented() const;
};

class PsiCI: public PsiWfn {
#   define CLASSNAME PsiCI
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  protected:
    void write_input(int conv);
    double read_energy();
  public:
    PsiCI(const RefKeyVal&);
    PsiCI(StateIn&);
    ~PsiCI();
    void save_data_state(StateOut&);
};

#endif

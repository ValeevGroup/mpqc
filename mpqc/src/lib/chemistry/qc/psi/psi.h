
#ifdef __GNUC__
#pragma interface
#endif

#ifndef _chemistry_qc_psi_psi_h
#define _chemistry_qc_psi_psi_h

#include <chemistry/qc/wfn/wfn.h>
#include <chemistry/qc/psi/psiinput.h>

class PsiWfn: public Wavefunction {
  protected:
    PSI_Input psi_in;
    void compute();

    virtual void write_input(int conv) = 0;
    virtual double read_energy() = 0;
    void write_basic_input(int conv, const char *wfn);
  public:
    PsiWfn(const Ref<KeyVal>&);
    PsiWfn(StateIn&);
    virtual ~PsiWfn();
    void save_data_state(StateOut&);

    double density(const SCVector3&);
    RefSymmSCMatrix density();

    void print(std::ostream&o=ExEnv::out0()) const;

    int spin_polarized();
    int nelectron();

    int gradient_implemented() const;
    int value_implemented() const;
};

class PsiHF: public PsiWfn {
  protected:
    void write_input(int conv);
    double read_energy();
  public:
    PsiHF(const Ref<KeyVal>&);
    PsiHF(StateIn&);
    ~PsiHF();
    void save_data_state(StateOut&);
};

class PsiCCSD: public PsiWfn {
  protected:
    void write_input(int conv);
    double read_energy();
  public:
    PsiCCSD(const Ref<KeyVal>&);
    PsiCCSD(StateIn&);
    ~PsiCCSD();
    void save_data_state(StateOut&);
};

class PsiCCSD_T: public PsiWfn {
  protected:
    void write_input(int conv);
    double read_energy();
  public:
    PsiCCSD_T(const Ref<KeyVal>&);
    PsiCCSD_T(StateIn&);
    ~PsiCCSD_T();
    void save_data_state(StateOut&);
};

class PsiCCSDT: public PsiWfn {
  protected:
    void write_input(int conv);
    double read_energy();
  public:
    PsiCCSDT(const Ref<KeyVal>&);
    PsiCCSDT(StateIn&);
    ~PsiCCSDT();
    void save_data_state(StateOut&);
    int gradient_implemented() const;
};

class PsiCI: public PsiWfn {
  protected:
    void write_input(int conv);
    double read_energy();
  public:
    PsiCI(const Ref<KeyVal>&);
    PsiCI(StateIn&);
    ~PsiCI();
    void save_data_state(StateOut&);
};

#endif

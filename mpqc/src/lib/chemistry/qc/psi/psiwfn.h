
#ifdef __GNUC__
#pragma interface
#endif

#ifndef _chemistry_qc_psi_psiwfn_h
#define _chemistry_qc_psi_psiwfn_h

#include <chemistry/qc/wfn/wfn.h>
#include <chemistry/qc/psi/psiexenv.h>

class PsiWavefunction: public Wavefunction {

    Ref<PsiExEnv> exenv_;

  protected:
    virtual void write_input(int conv) = 0;
    virtual void write_basic_input(int conv, const char *wfn);

  public:
    PsiWavefunction(const Ref<KeyVal>&);
    PsiWavefunction(StateIn&);
    ~PsiWavefunction();

    void save_data_state(StateOut&);
    void compute();
    void print(std::ostream&o=ExEnv::out0()) const;
    RefSymmSCMatrix density();
    int nelectron();
    Ref<PsiExEnv> get_psi_exenv() const { return exenv_; };
    Ref<PsiInput> get_psi_input() const { return exenv_->get_psi_input(); };
};

class PsiSCF: public PsiWavefunction {
  protected:
    void write_basic_input(int conv);
  public:
    PsiSCF(const Ref<KeyVal>&);
    PsiSCF(StateIn&);
    ~PsiSCF();
    void save_data_state(StateOut&);
};

class PsiCLHF: public PsiSCF {
  protected:
    void write_input(int conv);
  public:
    PsiCLHF(const Ref<KeyVal>&);
    PsiCLHF(StateIn&);
    ~PsiCLHF();
    int spin_polarized() { return 0;};
    int gradient_implemented() const { return 1;};
};

class PsiROHF: public PsiSCF {
  protected:
    void write_input(int conv);
  public:
    PsiROHF(const Ref<KeyVal>&);
    PsiROHF(StateIn&);
    ~PsiROHF();
    int spin_polarized() { return 0;};
    int gradient_implemented() const { return 1;};
};

class PsiCCSD: public PsiWavefunction {
  protected:
    void write_input(int conv);
  public:
    PsiCCSD(const Ref<KeyVal>&);
    PsiCCSD(StateIn&);
    ~PsiCCSD();
    void save_data_state(StateOut&);
    int spin_polarized() { return 0;};
    int gradient_implemented() const { return 1;};
};

class PsiCCSD_T: public PsiWavefunction {
  protected:
    void write_input(int conv);
  public:
    PsiCCSD_T(const Ref<KeyVal>&);
    PsiCCSD_T(StateIn&);
    ~PsiCCSD_T();
    void save_data_state(StateOut&);
    int spin_polarized() { return 0;};
    int gradient_implemented() const { return 0;};
};

#endif

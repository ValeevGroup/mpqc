
#ifdef __GNUC__
#pragma interface
#endif

#ifndef _chemistry_qc_psi_psiwfn_h
#define _chemistry_qc_psi_psiwfn_h

#include <chemistry/qc/wfn/wfn.h>
#include <chemistry/qc/psi/psiexenv.h>

class PsiWavefunction: public Wavefunction {

    Ref<PsiExEnv> exenv_;

    int* read_occ(const Ref<KeyVal> &keyval, const char *name, int nirrep);

  protected:
    int nirrep_;
    int *docc_;
    int *socc_;
    int *frozen_docc_;
    int *frozen_uocc_;
    char *memory_;
    virtual void write_input(int conv) = 0;

  public:
    PsiWavefunction(const Ref<KeyVal>&);
    PsiWavefunction(StateIn&);
    ~PsiWavefunction();

    void save_data_state(StateOut&);
    virtual void write_basic_input(int conv);
    void compute();
    void print(std::ostream&o=ExEnv::out0()) const;
    RefSymmSCMatrix density();
    int nelectron();
    Ref<PsiExEnv> get_psi_exenv() const { return exenv_; };
    Ref<PsiInput> get_psi_input() const { return exenv_->get_psi_input(); };
};

class PsiSCF: public PsiWavefunction {
  public:
    PsiSCF(const Ref<KeyVal>&);
    PsiSCF(StateIn&);
    ~PsiSCF();
    void save_data_state(StateOut&);

    enum RefType {rhf, rohf, uhf};
    virtual PsiSCF::RefType reftype() const =0;
};

class PsiCLHF: public PsiSCF {
  protected:
    void write_input(int conv);
  public:
    PsiCLHF(const Ref<KeyVal>&);
    PsiCLHF(StateIn&);
    ~PsiCLHF();

    void write_basic_input(int conv);
    int spin_polarized() { return 0;};
    int gradient_implemented() const { return 1;};
    PsiSCF::RefType reftype() const { return rhf;};
};

class PsiROHF: public PsiSCF {
  protected:
    void write_input(int conv);
  public:
    PsiROHF(const Ref<KeyVal>&);
    PsiROHF(StateIn&);
    ~PsiROHF();

    void write_basic_input(int conv);
    int spin_polarized() { return 0;};
    int gradient_implemented() const { return 1;};
    PsiSCF::RefType reftype() const { return rohf;};
};

class PsiUHF: public PsiSCF {
  protected:
    void write_input(int conv);
  public:
    PsiUHF(const Ref<KeyVal>&);
    PsiUHF(StateIn&);
    ~PsiUHF();

    void write_basic_input(int conv);
    int spin_polarized() { return 1;};
    int gradient_implemented() const { return 1;};
    PsiSCF::RefType reftype() const { return uhf;};
};

class PsiCCSD: public PsiWavefunction {
    Ref<PsiSCF> reference_;
  protected:
    void write_input(int conv);
  public:
    PsiCCSD(const Ref<KeyVal>&);
    PsiCCSD(StateIn&);
    ~PsiCCSD();
    void save_data_state(StateOut&);
    int spin_polarized() { return reference_->spin_polarized();};
    int gradient_implemented() const;
};

class PsiCCSD_T: public PsiWavefunction {
    Ref<PsiSCF> reference_;
  protected:
    void write_input(int conv);
  public:
    PsiCCSD_T(const Ref<KeyVal>&);
    PsiCCSD_T(StateIn&);
    ~PsiCCSD_T();

    void save_data_state(StateOut&);
    int spin_polarized() { return reference_->spin_polarized();};
    int gradient_implemented() const;
};

#endif


#ifdef __GNUC__
#pragma implementation
#endif

#include <util/keyval/keyval.h>
#include <util/misc/libmisc.h>
#include <util/misc/formio.h>
#include <util/state/stateio.h>
#include <math/scmat/matrix.h>
#include <math/symmetry/pointgrp.h>
#include <chemistry/molecule/molecule.h>
#include <chemistry/qc/psi/psiwfn.h>

using namespace std;

//////////////////////////////////////////////////////////////////////////

static ClassDesc PsiWavefunction_cd(
  typeid(PsiWavefunction),"PsiWavefunction",2,"public Wavefunction",
  0, 0, 0);

PsiWavefunction::PsiWavefunction(const Ref<KeyVal>&keyval):
  Wavefunction(keyval)
{
  exenv_ << keyval->describedclassvalue("psienv");
}

PsiWavefunction::~PsiWavefunction()
{
}

PsiWavefunction::PsiWavefunction(StateIn&s):
  SavableState(s),
  Wavefunction(s)
{
  abort();
}

void
PsiWavefunction::save_data_state(StateOut&s)
{
  abort();
}

void
PsiWavefunction::print(ostream&o) const
{
  Wavefunction::print(o);
  exenv_->print(o);
}

void
PsiWavefunction::compute()
{  
  if (gradient_needed() && !gradient_implemented()) {
    ExEnv::out0() << scprintf("Gradient is not implemented for this Psi wavefunction") << endl;
    abort();
  }
  double energy_acc = desired_value_accuracy();
  double grad_acc = desired_gradient_accuracy();
  if (energy_acc > 1.0e-6) energy_acc = 1.0e-6;
  if (grad_acc > 1.0e-7) grad_acc = 1.0e-7;
  if (gradient_needed() && energy_acc > grad_acc/10.0)
      energy_acc = grad_acc/10.0;

  write_input((int)-log10(energy_acc));
  exenv_->run_psi();

  // read output
  if (gradient_needed()) {
    Ref<PsiFile11> file11 = exenv_->get_psi_file11();
    file11->open();

    set_energy(file11->get_energy(0));
    set_actual_value_accuracy(energy_acc);

    int natom_mpqc = molecule()->natom();
    int natom = file11->get_natom(0);
    if (natom != natom_mpqc) {
      ExEnv::out0() << scprintf("Number of atoms in MPQC and Psi3 do not match") << endl;
      abort();
    }
    RefSCVector gradientvec = basis()->matrixkit()->vector(moldim());
    for(int atom=0; atom<natom; atom++) {
      gradientvec[3*atom] = file11->get_grad(0,atom,0);
      gradientvec[3*atom+1] = file11->get_grad(0,atom,1);
      gradientvec[3*atom+2] = file11->get_grad(0,atom,2);
    }
    set_gradient(gradientvec);
    file11->close();
    system("/bin/rm -f /tmp/file11.dat");
  }
  else {
      double energy = 0.0;;
      set_energy(energy);
      set_actual_value_accuracy(energy_acc);
  }
}

RefSymmSCMatrix
PsiWavefunction::density()
{
  abort();
  return 0;
}

int
PsiWavefunction::nelectron()
{
  abort();
  return 0;
}

void
PsiWavefunction::write_basic_input(int conv, const char *wfn)
{
  const char *dertype = gradient_needed() ? "first" : "none";

  Ref<PsiInput> psiinput = exenv_->get_psi_input();
  psiinput->write_defaults(exenv_,wfn,dertype);
  psiinput->begin_section("input");
  psiinput->write_keyword("no_reorient","true");
  psiinput->write_keyword("basis","ccpvdz");
  psiinput->write_geom(molecule());
  psiinput->end_section();
}

//////////////////////////////////////////////////////////////////////////

static ClassDesc PsiSCF_cd(
  typeid(PsiSCF),"PsiSCF",1,"public PsiWavefunction",
  0, 0, 0);

PsiSCF::PsiSCF(const Ref<KeyVal>&keyval):
  PsiWavefunction(keyval)
{
}

PsiSCF::~PsiSCF()
{
}

PsiSCF::PsiSCF(StateIn&s):
  SavableState(s),
  PsiWavefunction(s)
{
  abort();
}

void
PsiSCF::save_data_state(StateOut&s)
{
  PsiWavefunction::save_data_state(s);
  abort();
}

void
PsiSCF::write_basic_input(int convergence)
{
  PsiWavefunction::write_basic_input(convergence, "SCF");
}

//////////////////////////////////////////////////////////////////////////

static ClassDesc PsiCLHF_cd(
  typeid(PsiCLHF),"PsiCLHF",1,"public PsiSCF",
  0, create<PsiCLHF>, create<PsiCLHF>);

PsiCLHF::PsiCLHF(const Ref<KeyVal>&keyval):
  PsiSCF(keyval)
{
}

PsiCLHF::~PsiCLHF()
{
}

PsiCLHF::PsiCLHF(StateIn&s):
  PsiSCF(s)
{
  abort();
}

void
PsiCLHF::write_input(int convergence)
{
  Ref<PsiInput> input = get_psi_input();
  input->open();
  write_basic_input(convergence);
  input->write_keyword("default:reftype","rhf");
  input->close();
}

//////////////////////////////////////////////////////////////////////////

static ClassDesc PsiROHF_cd(
  typeid(PsiROHF),"PsiROHF",1,"public PsiSCF",
  0, create<PsiROHF>, create<PsiROHF>);

PsiROHF::PsiROHF(const Ref<KeyVal>&keyval):
  PsiSCF(keyval)
{
}

PsiROHF::~PsiROHF()
{
}

PsiROHF::PsiROHF(StateIn&s):
  PsiSCF(s)
{
  abort();
}

void
PsiROHF::write_input(int convergence)
{
  Ref<PsiInput> input = get_psi_input();
  input->open();
  write_basic_input(convergence);
  input->write_keyword("default:reftype","rohf");
  input->close();
}

//////////////////////////////////////////////////////////////////////////

static ClassDesc PsiCCSD_cd(
  typeid(PsiCCSD),"PsiCCSD",1,"public PsiWavefunction",
  0, create<PsiCCSD>, create<PsiCCSD>);

PsiCCSD::PsiCCSD(const Ref<KeyVal>&keyval):
  PsiWavefunction(keyval)
{
}

PsiCCSD::~PsiCCSD()
{
}

PsiCCSD::PsiCCSD(StateIn&s):
  SavableState(s),
  PsiWavefunction(s)
{
  abort();
}

void
PsiCCSD::save_data_state(StateOut&s)
{
  PsiWavefunction::save_data_state(s);
  abort();
}

void
PsiCCSD::write_input(int convergence)
{
  Ref<PsiInput> input = get_psi_input();
  input->open();
  write_basic_input(convergence,"ccsd");
  input->write_keyword("default:reftype","rhf");
  input->write_keyword("default:memory","(20 MB)");
  input->close();
}

//////////////////////////////////////////////////////////////////////////

static ClassDesc PsiCCSD_T_cd(
  typeid(PsiCCSD_T),"PsiCCSD_T",1,"public PsiWavefunction",
  0, create<PsiCCSD_T>, create<PsiCCSD_T>);

PsiCCSD_T::PsiCCSD_T(const Ref<KeyVal>&keyval):
  PsiWavefunction(keyval)
{
}

PsiCCSD_T::~PsiCCSD_T()
{
}

PsiCCSD_T::PsiCCSD_T(StateIn&s):
  SavableState(s),
  PsiWavefunction(s)
{
  abort();
}

void
PsiCCSD_T::save_data_state(StateOut&s)
{
  PsiWavefunction::save_data_state(s);
  abort();
}

void
PsiCCSD_T::write_input(int convergence)
{
  Ref<PsiInput> input = get_psi_input();
  input->open();
  write_basic_input(convergence,"ccsd");
  input->write_keyword("default:reftype","rhf");
  input->write_keyword("default:memory","(20 MB)");
  input->close();
}

//////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:

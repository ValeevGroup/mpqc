
#ifdef __GNUC__
#pragma implementation
#endif

#include <cmath>

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
  if (exenv_.null()) {
    ExEnv::err0() << "PsiWavefunction::PsiWavefunction: no Psi execution environment object (psienv)" << endl;
    abort();
  }
  
  nirrep_ = molecule()->point_group()->char_table().order();
  docc_ = read_occ(keyval,"docc",nirrep_);
  socc_ = read_occ(keyval,"socc",nirrep_);
  frozen_docc_ = read_occ(keyval,"frozen_docc",nirrep_);
  frozen_uocc_ = read_occ(keyval,"frozen_uocc",nirrep_);

  int bytes = keyval->intvalue("memory");
  if (bytes <= 2000000)
    bytes = 2000000;
  int bytes_str_len = (int)ceil(log10((long double)bytes));
  memory_ = new char[bytes_str_len+5];
  sprintf(memory_,"(%ld B)",bytes);
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
    file11->remove();
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
PsiWavefunction::write_basic_input(int conv)
{
  const char *dertype = gradient_needed() ? "first" : "none";

  Ref<PsiInput> psiinput = get_psi_input();
  psiinput->write_defaults(exenv_,dertype);
  psiinput->write_keyword("default:memory",memory_);
  psiinput->begin_section("input");
  psiinput->write_keyword("no_reorient","true");
  psiinput->write_basis(basis());
  if (basis()->max_nfunction_in_shell() != basis()->max_ncartesian_in_shell())
    psiinput->write_keyword("puream","true");
  psiinput->write_geom(molecule());
  psiinput->end_section();
  psiinput->write_basis_sets(basis());
}

// Shamelessly borrowed from class SCF
int *
PsiWavefunction::read_occ(const Ref<KeyVal> &keyval, const char *name, int nirrep)
{
  int *occ = 0;
  if (keyval->exists(name)) {
    if (keyval->count(name) != nirrep) {
      ExEnv::err0() << indent
                   << "ERROR: PsiWavefunction: have " << nirrep << " irreps but "
                   << name << " vector is length " << keyval->count(name)
                   << endl;
      abort();
    }
    occ = new int[nirrep];
    for (int i=0; i<nirrep; i++) {
      occ[i] = keyval->intvalue(name,i);
    }
  }
  return occ;
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

//////////////////////////////////////////////////////////////////////////

static ClassDesc PsiCLHF_cd(
  typeid(PsiCLHF),"PsiCLHF",1,"public PsiSCF",
  0, create<PsiCLHF>, create<PsiCLHF>);

PsiCLHF::PsiCLHF(const Ref<KeyVal>&keyval):
  PsiSCF(keyval)
{
  if (!docc_) {
    ExEnv::err0() << indent
		  << "ERROR: PsiCLHF: did not find docc array" << endl;
    abort();
  }
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
PsiCLHF::write_basic_input(int convergence)
{
  Ref<PsiInput> input = get_psi_input();
  input->write_keyword("default:reference","rhf");
  if (docc_)
    input->write_keyword_array("default:docc",nirrep_,docc_);
}

void
PsiCLHF::write_input(int convergence)
{
  Ref<PsiInput> input = get_psi_input();
  input->open();
  PsiWavefunction::write_basic_input(convergence);
  write_basic_input(convergence);
  input->write_keyword("default:wfn","scf");
  input->close();
}

//////////////////////////////////////////////////////////////////////////

static ClassDesc PsiHSOSHF_cd(
  typeid(PsiHSOSHF),"PsiHSOSHF",1,"public PsiSCF",
  0, create<PsiHSOSHF>, create<PsiHSOSHF>);

PsiHSOSHF::PsiHSOSHF(const Ref<KeyVal>&keyval):
  PsiSCF(keyval)
{
  if (!docc_) {
    ExEnv::err0() << indent
		  << "ERROR: PsiHSOSHF: did not find docc array" << endl;
    abort();
  }
  if (!socc_) {
    ExEnv::err0() << indent
		  << "ERROR: PsiHSOSHF: did not find socc array" << endl;
    abort();
  }
}

PsiHSOSHF::~PsiHSOSHF()
{
}

PsiHSOSHF::PsiHSOSHF(StateIn&s):
  PsiSCF(s)
{
  abort();
}

void
PsiHSOSHF::write_basic_input(int convergence)
{
  Ref<PsiInput> input = get_psi_input();
  input->write_keyword("default:reference","rohf");
  if (docc_)
    input->write_keyword_array("default:docc",nirrep_,docc_);
  if (socc_)
    input->write_keyword_array("default:socc",nirrep_,socc_);
}

void
PsiHSOSHF::write_input(int convergence)
{
  Ref<PsiInput> input = get_psi_input();
  input->open();
  PsiWavefunction::write_basic_input(convergence);
  write_basic_input(convergence);
  input->write_keyword("default:wfn","scf");
  input->close();
}


//////////////////////////////////////////////////////////////////////////

static ClassDesc PsiUHF_cd(
  typeid(PsiUHF),"PsiUHF",1,"public PsiSCF",
  0, create<PsiUHF>, create<PsiUHF>);

PsiUHF::PsiUHF(const Ref<KeyVal>&keyval):
  PsiSCF(keyval)
{
  if (!docc_) {
    ExEnv::err0() << indent
		  << "ERROR: PsiUHF: did not find docc array" << endl;
    abort();
  }
  if (!socc_) {
    ExEnv::err0() << indent
		  << "ERROR: PsiUHF: did not find socc array" << endl;
    abort();
  }
}

PsiUHF::~PsiUHF()
{
}

PsiUHF::PsiUHF(StateIn&s):
  PsiSCF(s)
{
  abort();
}

void
PsiUHF::write_basic_input(int convergence)
{
  Ref<PsiInput> input = get_psi_input();
  input->write_keyword("default:reference","uhf");
  if (docc_)
    input->write_keyword_array("default:docc",nirrep_,docc_);
  if (socc_)
    input->write_keyword_array("default:socc",nirrep_,socc_);
}

void
PsiUHF::write_input(int convergence)
{
  Ref<PsiInput> input = get_psi_input();
  input->open();
  PsiWavefunction::write_basic_input(convergence);
  write_basic_input(convergence);
  input->write_keyword("default:wfn","scf");
  input->close();
}

//////////////////////////////////////////////////////////////////////////

static ClassDesc PsiCCSD_cd(
  typeid(PsiCCSD),"PsiCCSD",1,"public PsiWavefunction",
  0, create<PsiCCSD>, create<PsiCCSD>);

PsiCCSD::PsiCCSD(const Ref<KeyVal>&keyval):
  PsiWavefunction(keyval)
{
  reference_ << keyval->describedclassvalue("reference");
  if (reference_.null()) {
    ExEnv::err0() << "PsiCCSD::PsiCCSD: no reference wavefunction" << endl;
    abort();
  }
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

int
PsiCCSD::gradient_implemented() const
{
  int impl = 0;
  PsiSCF::RefType reftype = reference_->reftype();
  if (reftype == PsiSCF::rhf || reftype == PsiSCF::hsoshf)
    impl = 1;
  return impl;
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
  if (gradient_needed())
    reference_->do_gradient(1);
  else
    reference_->do_gradient(0);
    
  Ref<PsiInput> input = get_psi_input();
  input->open();
  PsiWavefunction::write_basic_input(convergence);
  reference_->write_basic_input(convergence);
  input->write_keyword("default:wfn","ccsd");
  if (frozen_docc_)
    input->write_keyword_array("default:frozen_docc",nirrep_,frozen_docc_);
  if (frozen_uocc_)
    input->write_keyword_array("default:frozen_uocc",nirrep_,frozen_uocc_);
  input->close();
}

//////////////////////////////////////////////////////////////////////////

static ClassDesc PsiCCSD_T_cd(
  typeid(PsiCCSD_T),"PsiCCSD_T",1,"public PsiWavefunction",
  0, create<PsiCCSD_T>, create<PsiCCSD_T>);

PsiCCSD_T::PsiCCSD_T(const Ref<KeyVal>&keyval):
  PsiWavefunction(keyval)
{
  reference_ << keyval->describedclassvalue("reference");
  if (reference_.null()) {
    ExEnv::err0() << "PsiCCSD_T::PsiCCSD_T: no reference wavefunction" << endl;
    abort();
  }

  PsiSCF::RefType reftype = reference_->reftype();
  if (reftype == PsiSCF::hsoshf) {
    ExEnv::err0() << "PsiCCSD_T::PsiCCSD_T: HSOSHF-based CCSD(T) has not been implemented yet" << endl;
    abort();
  }
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

int
PsiCCSD_T::gradient_implemented() const
{
  int impl = 0;
  PsiSCF::RefType reftype = reference_->reftype();
  return impl;
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
  if (gradient_needed())
    reference_->do_gradient(1);
  else
    reference_->do_gradient(0);
    
  Ref<PsiInput> input = get_psi_input();
  input->open();
  PsiWavefunction::write_basic_input(convergence);
  reference_->write_basic_input(convergence);
  input->write_keyword("default:wfn","ccsd");
  input->begin_section("psi");
  input->write_keyword("exec","(\"cints\" \"cscf\" \"transqt\" \"ccsort\" \"ccenergy\" \"cchbar\" \"cctriples\")");
  input->end_section();
  if (frozen_docc_)
    input->write_keyword_array("default:frozen_docc",nirrep_,frozen_docc_);
  if (frozen_uocc_)
    input->write_keyword_array("default:frozen_uocc",nirrep_,frozen_uocc_);
  input->close();
}

//////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:

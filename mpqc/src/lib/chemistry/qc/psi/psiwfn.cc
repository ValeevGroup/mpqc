
#ifdef __GNUC__
#pragma implementation
#endif

#include <stdexcept>
#include <math.h>

#include <util/keyval/keyval.h>
#include <util/misc/formio.h>
#include <util/state/stateio.h>
#include <math/scmat/matrix.h>
#include <math/symmetry/pointgrp.h>
#include <chemistry/molecule/molecule.h>
#include <chemistry/qc/psi/psiwfn.h>

using namespace std;

namespace sc {

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
  throw std::runtime_error("PsiWavefunction::PsiWavefunction(StateIn&) -- cannot restore state of Psi wave functions");
}

void
PsiWavefunction::save_data_state(StateOut&s)
{
  throw std::runtime_error("PsiWavefunction::save_data_state -- cannot save state of Psi wave functions, set savestate = no in your input file");
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
  psiinput->write_keyword("psi:memory",memory_);
  psiinput->begin_section("input");
  psiinput->write_keyword("no_reorient","true");
  psiinput->write_keyword("keep_ref_frame","true");
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
  if (!docc_ || !socc_) {
    if (keyval->exists("total_charge") && keyval->exists("multiplicity")) {
      charge_ = keyval->intvalue("total_charge");
      multp_ = keyval->intvalue("multiplicity");
      if (multp_ < 1) {
	ExEnv::err0() << indent
		      << "ERROR: PsiSCF: valid multiplicity has to be >= 1" << endl;
	abort();
      }
    }
    else {
      ExEnv::err0() << indent
		    << "ERROR: PsiSCF: multiplicity and total_charge need "
		    << "to be specified when docc (socc) are missing" << endl;
      abort();
    }
  }
}

PsiSCF::~PsiSCF()
{
}

PsiSCF::PsiSCF(StateIn&s):
  PsiWavefunction(s)
{
}

void
PsiSCF::save_data_state(StateOut&s)
{
  PsiWavefunction::save_data_state(s);
}

//////////////////////////////////////////////////////////////////////////

static ClassDesc PsiCLHF_cd(
  typeid(PsiCLHF),"PsiCLHF",1,"public PsiSCF",
  0, create<PsiCLHF>, create<PsiCLHF>);

PsiCLHF::PsiCLHF(const Ref<KeyVal>&keyval):
  PsiSCF(keyval)
{
  if (!docc_ && multp_ != 1) {
    ExEnv::err0() << indent
		  << "ERROR: PsiCLHF: multiplicity should be 1 for CLHF wave function" << endl;
    abort();
  }
}

PsiCLHF::~PsiCLHF()
{
}

PsiCLHF::PsiCLHF(StateIn&s):
  PsiSCF(s)
{
}

void
PsiCLHF::write_basic_input(int convergence)
{
  Ref<PsiInput> input = get_psi_input();
  input->write_keyword("psi:reference","rhf");
  if (docc_)
    input->write_keyword_array("psi:docc",nirrep_,docc_);
  else {
    input->write_keyword("psi:multp",multp_);
    input->write_keyword("psi:charge",charge_);
  }
}

void
PsiCLHF::write_input(int convergence)
{
  Ref<PsiInput> input = get_psi_input();
  input->open();
  PsiWavefunction::write_basic_input(convergence);
  write_basic_input(convergence);
  input->write_keyword("psi:wfn","scf");
  input->close();
}

//////////////////////////////////////////////////////////////////////////

static ClassDesc PsiHSOSHF_cd(
  typeid(PsiHSOSHF),"PsiHSOSHF",1,"public PsiSCF",
  0, create<PsiHSOSHF>, create<PsiHSOSHF>);

PsiHSOSHF::PsiHSOSHF(const Ref<KeyVal>&keyval):
  PsiSCF(keyval)
{
  if ((!docc_ || !socc_) && multp_ == 1) {
    ExEnv::err0() << indent
		  << "ERROR: PsiHSOSHF: multiplicity should be > 1 for HSOSHF wave function" << endl;
    abort();
  }
}

PsiHSOSHF::~PsiHSOSHF()
{
}

PsiHSOSHF::PsiHSOSHF(StateIn&s):
  PsiSCF(s)
{
}

void
PsiHSOSHF::write_basic_input(int convergence)
{
  Ref<PsiInput> input = get_psi_input();
  input->write_keyword("psi:reference","rohf");
  if (docc_)
    input->write_keyword_array("psi:docc",nirrep_,docc_);
  if (socc_)
    input->write_keyword_array("psi:socc",nirrep_,socc_);
}

void
PsiHSOSHF::write_input(int convergence)
{
  Ref<PsiInput> input = get_psi_input();
  input->open();
  PsiWavefunction::write_basic_input(convergence);
  write_basic_input(convergence);
  input->write_keyword("psi:wfn","scf");
  input->close();
}


//////////////////////////////////////////////////////////////////////////

static ClassDesc PsiUHF_cd(
  typeid(PsiUHF),"PsiUHF",1,"public PsiSCF",
  0, create<PsiUHF>, create<PsiUHF>);

PsiUHF::PsiUHF(const Ref<KeyVal>&keyval):
  PsiSCF(keyval)
{
}

PsiUHF::~PsiUHF()
{
}

PsiUHF::PsiUHF(StateIn&s):
  PsiSCF(s)
{
}

void
PsiUHF::write_basic_input(int convergence)
{
  Ref<PsiInput> input = get_psi_input();
  input->write_keyword("psi:reference","uhf");
  if (docc_)
    input->write_keyword_array("psi:docc",nirrep_,docc_);
  if (socc_)
    input->write_keyword_array("psi:socc",nirrep_,socc_);
}

void
PsiUHF::write_input(int convergence)
{
  Ref<PsiInput> input = get_psi_input();
  input->open();
  PsiWavefunction::write_basic_input(convergence);
  write_basic_input(convergence);
  input->write_keyword("psi:wfn","scf");
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
  PsiWavefunction(s)
{
  reference_ << SavableState::restore_state(s);
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
  SavableState::save_state(reference_.pointer(),s);
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
  input->write_keyword("psi:wfn","ccsd");
  if (frozen_docc_)
    input->write_keyword_array("psi:frozen_docc",nirrep_,frozen_docc_);
  if (frozen_uocc_)
    input->write_keyword_array("psi:frozen_uocc",nirrep_,frozen_uocc_);
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
  PsiWavefunction(s)
{
  reference_ << SavableState::restore_state(s);
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
  SavableState::save_state(reference_.pointer(),s);
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
  input->write_keyword("psi:wfn","ccsd");
  input->begin_section("psi");
  input->write_keyword("exec","(\"cints\" \"cscf\" \"transqt\" \"ccsort\" \"ccenergy\" \"cchbar\" \"cctriples\")");
  input->end_section();
  if (frozen_docc_)
    input->write_keyword_array("psi:frozen_docc",nirrep_,frozen_docc_);
  if (frozen_uocc_)
    input->write_keyword_array("psi:frozen_uocc",nirrep_,frozen_uocc_);
  input->close();
}

//////////////////////////////////////////////////////////////////////////

}

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:

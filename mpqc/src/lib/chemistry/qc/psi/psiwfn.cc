
#ifdef __GNUC__
#pragma implementation
#endif

#include <stdexcept>
#include <cmath>

#include <psifiles.h>
#include <ccfiles.h>

#include <util/keyval/keyval.h>
#include <util/misc/formio.h>
#include <util/class/scexception.h>
#include <util/state/stateio.h>
#include <math/scmat/matrix.h>
#include <math/symmetry/pointgrp.h>
#include <chemistry/molecule/molecule.h>
#include <chemistry/qc/mbptr12/mbptr12.h>
#include <chemistry/qc/mbptr12/spin.h>
#include <chemistry/qc/mbptr12/utils.h>
#include <chemistry/qc/psi/psiwfn.h>
#include <chemistry/qc/psi/psiinput.timpl.h>

#define TEST_V 0

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

  size_t bytes = keyval->sizevalue("memory");
  if (bytes <= 2000000)
    bytes = 2000000;
  int bytes_str_len = (int)ceil(log10((long double)bytes));
  memory_ = new char[bytes_str_len+5];
  sprintf(memory_,"(%ld B)",bytes);
}

PsiWavefunction::~PsiWavefunction()
{
    exenv_->run_psi_module("psiclean");
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
std::vector<int>
PsiWavefunction::read_occ(const Ref<KeyVal> &keyval, const char *name, int nirrep)
{
  std::vector<int> occ;
  if (keyval->exists(name)) {
    if (keyval->count(name) != nirrep) {
      ExEnv::err0() << indent
                   << "ERROR: PsiWavefunction: have " << nirrep << " irreps but "
                   << name << " vector is length " << keyval->count(name)
                   << endl;
      abort();
    }
    occ.resize(nirrep);
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
  if (docc_.empty() || socc_.empty()) {
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

unsigned int
PsiSCF::nmo()
{
    psi::PSIO& psio = exenv()->psio();
    int num_mo;
    psio.open(PSIF_CHKPT,PSIO_OPEN_OLD);
    psio.read_entry(PSIF_CHKPT,"::Num. MO's",reinterpret_cast<char*>(&num_mo),sizeof(int));
    psio.close(PSIF_CHKPT,1);
    return num_mo;
}

unsigned int
PsiSCF::nocc(SpinCase1 spin)
{
    int* doccpi = new int[nirrep_];
    psi::PSIO& psio = exenv()->psio();
    psio.open(PSIF_CHKPT,PSIO_OPEN_OLD);
    psio.read_entry(PSIF_CHKPT,"::Closed shells per irrep",reinterpret_cast<char*>(doccpi),nirrep_*sizeof(int));
    psio.close(PSIF_CHKPT,1);

    unsigned int nocc = 0;
    for(unsigned int h=0; h<nirrep_; ++h)
	nocc += doccpi[h];
    delete[] doccpi;

    return nocc;
}

const RefDiagSCMatrix&
PsiSCF::evals()
{
    if (evals_.nonnull())
	return evals_;

    // For now only handle closed-shell case
    PsiSCF::RefType ref = reftype();
    if (ref != PsiSCF::rhf)
	throw FeatureNotImplemented("PsiSCF::coefs() -- only closed-shell case is implemented",__FILE__,__LINE__);

    psi::PSIO& psio = exenv()->psio();
    // grab orbital info
    int num_mo;
    int* mopi = new int[nirrep_];
    psio.open(PSIF_CHKPT,PSIO_OPEN_OLD);
    psio.read_entry(PSIF_CHKPT,"::Num. MO's",reinterpret_cast<char*>(&num_mo),sizeof(int));
    psio.read_entry(PSIF_CHKPT,"::MO's per irrep",reinterpret_cast<char*>(mopi),nirrep_*sizeof(int));
    // get the eigenvalues
    double* E = new double[num_mo];
    psio.read_entry(PSIF_CHKPT,"::MO energies",reinterpret_cast<char*>(E),num_mo*sizeof(double));
    psio.close(PSIF_CHKPT,1);

    // convert raw matrices to SCMatrices
    RefSCDimension modim = new SCDimension(num_mo,nirrep_,mopi);
    for(unsigned int h=0; h<nirrep_; ++h)
	modim->blocks()->set_subdim(h,new SCDimension(mopi[h]));
    evals_ = basis_matrixkit()->diagmatrix(modim);  evals_.assign(E);
    evals_.print("Psi3 SCF eigenvalues");

    delete[] mopi;

    return evals_;
}

const RefSCMatrix&
PsiSCF::coefs()
{
    if (coefs_.nonnull())
	return coefs_;

    // For now only handle closed-shell case
    PsiSCF::RefType ref = reftype();
    if (ref != PsiSCF::rhf)
	throw FeatureNotImplemented("PsiSCF::coefs() -- only closed-shell case is implemented",__FILE__,__LINE__);

    psi::PSIO& psio = exenv()->psio();
    // grab orbital info
    int num_so, num_mo;
    int* mopi = new int[nirrep_];
    int* sopi = new int[nirrep_];
    psio.open(PSIF_CHKPT,PSIO_OPEN_OLD);
    psio.read_entry(PSIF_CHKPT,"::Num. SO",reinterpret_cast<char*>(&num_so),sizeof(int));
    psio.read_entry(PSIF_CHKPT,"::Num. MO's",reinterpret_cast<char*>(&num_mo),sizeof(int));
    psio.read_entry(PSIF_CHKPT,"::SO's per irrep",reinterpret_cast<char*>(sopi),nirrep_*sizeof(int));
    psio.read_entry(PSIF_CHKPT,"::MO's per irrep",reinterpret_cast<char*>(mopi),nirrep_*sizeof(int));
    // get MO coefficients in SO basis
    const unsigned int nsm = num_so * num_mo;
    double* C = new double[nsm];
    psio.read_entry(PSIF_CHKPT,"::MO coefficients",reinterpret_cast<char*>(C),nsm*sizeof(double));
    // get AO->SO matrix (MPQC AO equiv PSI3 BF)
    const unsigned int nss = num_so * num_so;
    double* ao2so = new double[nss];
    psio.read_entry(PSIF_CHKPT,"::SO->BF transmat",reinterpret_cast<char*>(ao2so),nsm*sizeof(double));
    psio.close(PSIF_CHKPT,1);

    // convert raw matrices to SCMatrices
    RefSCDimension sodim_nb = new SCDimension(num_so,1);
    sodim_nb->blocks()->set_subdim(0,new SCDimension(num_so));
    RefSCDimension sodim = new SCDimension(num_so,nirrep_,sopi);
    for(unsigned int h=0; h<nirrep_; ++h)
	sodim->blocks()->set_subdim(h,new SCDimension(sopi[h]));
    RefSCDimension modim = new SCDimension(num_mo,nirrep_,mopi);
    for(unsigned int h=0; h<nirrep_; ++h)
	modim->blocks()->set_subdim(h,new SCDimension(mopi[h]));
    RefSCMatrix C_so = basis_matrixkit()->matrix(sodim,modim);  C_so.assign(C);
    C_so.print("Psi3 eigenvector in SO basis");
    RefSCMatrix aotoso = basis_matrixkit()->matrix(sodim,sodim_nb);  aotoso.assign(ao2so);
    aotoso.print("Psi3 SO->AO matrix");
    coefs_ = aotoso.t() * C_so;
    coefs_.print("Psi3 eigenvector in AO basis");

    delete[] mopi;
    delete[] sopi;

    return coefs_;
}

//////////////////////////////////////////////////////////////////////////

static ClassDesc PsiCLHF_cd(
  typeid(PsiCLHF),"PsiCLHF",1,"public PsiSCF",
  0, create<PsiCLHF>, create<PsiCLHF>);

PsiCLHF::PsiCLHF(const Ref<KeyVal>&keyval):
  PsiSCF(keyval)
{
  if (docc_.empty() && multp_ != 1) {
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
  if (!docc_.empty())
    input->write_keyword_array("psi:docc",nirrep_,docc_);
  else {
    input->write_keyword("psi:multp",multp_);
    input->write_keyword("psi:charge",charge_);
    input->write_keyword("psi:reset_occupations",true);
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
  if ((docc_.empty() || socc_.empty()) && multp_ == 1) {
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
  if (!docc_.empty())
    input->write_keyword_array("psi:docc",nirrep_,docc_);
  if (!socc_.empty())
    input->write_keyword_array("psi:socc",nirrep_,socc_);
  if (docc_.empty() && socc_.empty()) {
    input->write_keyword("psi:multp",multp_);
    input->write_keyword("psi:charge",charge_);
    input->write_keyword("psi:reset_occupations",true);
  }
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
  if (!docc_.empty())
    input->write_keyword_array("psi:docc",nirrep_,docc_);
  if (!socc_.empty())
    input->write_keyword_array("psi:socc",nirrep_,socc_);
  if (docc_.empty() && socc_.empty()) {
    input->write_keyword("psi:multp",multp_);
    input->write_keyword("psi:charge",charge_);
    input->write_keyword("psi:reset_occupations",true);
  }
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

static ClassDesc PsiCorrWavefunction_cd(
  typeid(PsiCorrWavefunction),"PsiCorrWavefunction",1,"public PsiWavefunction",
  0, create<PsiCorrWavefunction>, create<PsiCorrWavefunction>);

PsiCorrWavefunction::PsiCorrWavefunction(const Ref<KeyVal>&keyval):
  PsiWavefunction(keyval)
{
  frozen_docc_ = read_occ(keyval,"frozen_docc",nirrep_);
  frozen_uocc_ = read_occ(keyval,"frozen_uocc",nirrep_);
  if (frozen_docc_.empty()) {
      frozen_docc_.resize(nirrep_);
      for(unsigned int h=0; h<nirrep_; ++h)
	  frozen_docc_[h] = 0;
  }
  if (frozen_uocc_.empty()) {
      frozen_uocc_.resize(nirrep_);
      for(unsigned int h=0; h<nirrep_; ++h)
	  frozen_uocc_[h] = 0;
  }

  reference_ << keyval->describedclassvalue("reference");
  if (reference_.null()) {
    ExEnv::err0() << "PsiCorrWavefunction::PsiCorrWavefunction: no reference wavefunction" << endl;
    abort();
  }
}

PsiCorrWavefunction::~PsiCorrWavefunction()
{
}

PsiCorrWavefunction::PsiCorrWavefunction(StateIn&s):
  PsiWavefunction(s)
{
  reference_ << SavableState::restore_state(s);
  s.get(frozen_docc_);
  s.get(frozen_uocc_);
}

void
PsiCorrWavefunction::save_data_state(StateOut&s)
{
  PsiWavefunction::save_data_state(s);
  SavableState::save_state(reference_.pointer(),s);
  s.put(frozen_docc_);
  s.put(frozen_uocc_);
}

void
PsiCorrWavefunction::write_input(int convergence)
{
  if (gradient_needed())
    reference_->do_gradient(1);
  else
    reference_->do_gradient(0);
    
  Ref<PsiInput> input = get_psi_input();
  PsiWavefunction::write_basic_input(convergence);
  reference_->write_basic_input(convergence);
  if (!frozen_docc_.empty())
    input->write_keyword_array("psi:frozen_docc",nirrep_,frozen_docc_);
  if (!frozen_uocc_.empty())
    input->write_keyword_array("psi:frozen_uocc",nirrep_,frozen_uocc_);
}

const Ref<MOIndexSpace>&
PsiCorrWavefunction::occ_act_sb(SpinCase1 spin)
{
    if (occ_act_sb_.nonnull())
	return occ_act_sb_;

    const int nmo = reference_->nmo();
    const int nocc = reference_->nocc(spin);
    int nfzc = 0;
    for(unsigned int h=0; h<nirrep_; ++h)
	nfzc += frozen_docc_[h];

    occ_act_sb_ = new MOIndexSpace("i","active occupied MOs (Psi3)",
				   reference_->coefs(),basis(),integral(),
				   reference_->evals(),nfzc,nmo-nocc,MOIndexSpace::symmetry);
    
    return occ_act_sb_;
}

const Ref<MOIndexSpace>&
PsiCorrWavefunction::vir_act_sb(SpinCase1 spin)
{
    if (vir_act_sb_.nonnull())
	return vir_act_sb_;

    const int nmo = reference_->nmo();
    const int nocc = reference_->nocc(spin);
    int nfzv = 0;
    for(unsigned int h=0; h<nirrep_; ++h)
	nfzv += frozen_uocc_[h];

    vir_act_sb_ = new MOIndexSpace("a","active virtual MOs (Psi3)",
				   reference_->coefs(),basis(),integral(),
				   reference_->evals(),nocc,nfzv,MOIndexSpace::symmetry);
    
    return vir_act_sb_;
}

//////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////

static ClassDesc PsiCCSD_cd(
  typeid(PsiCCSD),"PsiCCSD",1,"public PsiCC",
  0, create<PsiCCSD>, create<PsiCCSD>);

PsiCCSD::PsiCCSD(const Ref<KeyVal>&keyval):
  PsiCC(keyval)
{
}

PsiCCSD::~PsiCCSD()
{
}

PsiCCSD::PsiCCSD(StateIn&s):
  PsiCC(s)
{
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
  PsiCC::save_data_state(s);
}

void
PsiCCSD::write_input(int convergence)
{
  Ref<PsiInput> input = get_psi_input();
  input->open();
  PsiCorrWavefunction::write_input(convergence);
  input->write_keyword("psi:wfn","ccsd");
  input->close();
}

//////////////////////////////////////////////////////////////////////////

static ClassDesc PsiCCSD_T_cd(
  typeid(PsiCCSD_T),"PsiCCSD_T",1,"public PsiCC",
  0, create<PsiCCSD_T>, create<PsiCCSD_T>);

PsiCCSD_T::PsiCCSD_T(const Ref<KeyVal>&keyval):
  PsiCC(keyval)
{
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
  PsiCC(s)
{
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
  PsiCC::save_data_state(s);
  SavableState::save_state(reference_.pointer(),s);
}

void
PsiCCSD_T::write_input(int convergence)
{
  Ref<PsiInput> input = get_psi_input();
  input->open();
  PsiCorrWavefunction::write_input(convergence);
  input->write_keyword("psi:wfn","ccsd_t");
  input->close();
}

//////////////////////////////////////////////////////////////////////////

static ClassDesc PsiCCSD_PT2R12_cd(
  typeid(PsiCCSD_PT2R12),"PsiCCSD_PT2R12",1,"public PsiCC",
  0, create<PsiCCSD_PT2R12>, create<PsiCCSD_PT2R12>);

PsiCCSD_PT2R12::PsiCCSD_PT2R12(const Ref<KeyVal>&keyval):
  PsiCC(keyval)
{
    if (mp2_only_)
	PsiCC::do_test_t2_phases();

    mbptr12_ =  require_dynamic_cast<MBPT2_R12*>(
	keyval->describedclassvalue("mbpt2r12").pointer(),
	"PsiCCSD_PT2R12::PsiCCSD_PT2R12\n"
	);

    const Ref<R12IntEvalInfo> r12info = mbptr12_->r12evalinfo();
    const Ref<R12Technology> r12tech = r12info->r12tech();
    // cannot do gbc = false yet
    if (!r12tech->gbc())
	throw FeatureNotImplemented("PsiCCSD_PT2R12::PsiCCSD_PT2R12() -- gbc = false is not yet implemented",__FILE__,__LINE__);
}

PsiCCSD_PT2R12::~PsiCCSD_PT2R12()
{
}

PsiCCSD_PT2R12::PsiCCSD_PT2R12(StateIn&s):
  PsiCC(s)
{
    mbptr12_ << SavableState::restore_state(s);
}

int
PsiCCSD_PT2R12::gradient_implemented() const
{
  return 0;
}

void
PsiCCSD_PT2R12::save_data_state(StateOut&s)
{
  PsiCC::save_data_state(s);
  SavableState::save_state(mbptr12_.pointer(),s);
}

void
PsiCCSD_PT2R12::write_input(int convergence)
{
  Ref<PsiInput> input = get_psi_input();
  input->open();
  PsiCorrWavefunction::write_input(convergence);
  // if testing T2 transform, obtain MP1 amplitudes
  if (test_t2_phases_)
      input->write_keyword("psi:wfn","mp2");
  else
      input->write_keyword("psi:wfn","ccsd");
  input->close();
}

void
PsiCCSD_PT2R12::compute()
{
    const Ref<R12IntEvalInfo> r12info = mbptr12_->r12evalinfo();
    const Ref<R12Technology> r12tech = r12info->r12tech();
    // params
    const bool ebc = r12tech->ebc();

    // compute CCSD wave function
    PsiWavefunction::compute();

    // grab amplitudes
    RefSCMatrix T1_psi = T(1);
    RefSCMatrix T2_psi = T(2);
    RefSCMatrix Tau2_psi;
    if (mp2_only_) {
      Tau2_psi = T2_psi.clone();  Tau2_psi.assign(T2_psi);
    }
    else
      Tau2_psi = Tau2();

    // Compute intermediates
    const double mp2r12_energy = mbptr12_->value();
    Ref<R12IntEval> r12eval = mbptr12_->r12eval();
    RefSCMatrix Vpq[NSpinCases2];
    RefSCMatrix Vab[NSpinCases2];
    RefSCMatrix Via[NSpinCases2];
    RefSCMatrix Vij[NSpinCases2];
    RefSymmSCMatrix B[NSpinCases2];
    RefSymmSCMatrix X[NSpinCases2];
    RefSCMatrix T2_MP1[NSpinCases2];
    RefSCMatrix A[NSpinCases2];
    RefSCMatrix T1[NSpinCases2];
    RefSCMatrix T2[NSpinCases2];
    RefSCMatrix Tau2[NSpinCases2];

    const int nspincases2 = r12eval->nspincases2();
    for(int s=0; s<nspincases2; s++) {
	const SpinCase2 spincase2 = static_cast<SpinCase2>(s);
	const SpinCase1 spin1 = case1(spincase2);
	const SpinCase1 spin2 = case2(spincase2);

	Ref<R12IntEvalInfo> r12info = r12eval->r12info();
	const Ref<MOIndexSpace>& p1 = r12info->refinfo()->orbs(spin1);
	const Ref<MOIndexSpace>& p2 = r12info->refinfo()->orbs(spin2);
	const unsigned int np1 = p1->rank();
	const unsigned int np2 = p2->rank();

	const Ref<MOIndexSpace>& occ1_act = r12eval->occ_act(spin1);
	const Ref<MOIndexSpace>& occ2_act = r12eval->occ_act(spin2);
	const Ref<MOIndexSpace>& vir1_act = r12eval->vir_act(spin1);
	const Ref<MOIndexSpace>& vir2_act = r12eval->vir_act(spin2);

	Vpq[s] = r12eval->V(spincase2,p1,p2);
	Vij[s] = r12eval->V(spincase2);
	X[s] = r12eval->X(spincase2);
	B[s] = r12eval->B(spincase2);
	T2_MP1[s] = r12eval->T2(spincase2);
	if (!ebc)
	    A[s] = r12eval->A(spincase2);

	const Ref<MOIndexSpace>& v1 = r12eval->vir_act(spin1);
	const Ref<MOIndexSpace>& v2 = r12eval->vir_act(spin2);
	const unsigned int nv1 = v1->rank();
	const unsigned int nv2 = v2->rank();
	const Ref<MOIndexSpace>& o1 = r12eval->occ_act(spin1);
	const Ref<MOIndexSpace>& o2 = r12eval->occ_act(spin2);
	const unsigned int no1 = o1->rank();
	const unsigned int no2 = o2->rank();
	const unsigned int no1v2 = no1*nv2;
	const unsigned int nxy = Vpq[s].rowdim().n();

	// extract Vab and Via from Vpq
	Vab[s] = Vpq[s].kit()->matrix(Vpq[s].rowdim(),T2_MP1[s].coldim());
	Via[s] = Vpq[s].kit()->matrix(Vpq[s].rowdim(),new SCDimension(no1v2));

	typedef MOIndexMap OrbMap;
	OrbMap v1_to_p1 (*p1<<*v1);
	OrbMap v2_to_p2 (*p2<<*v2);
	OrbMap o1_to_p1 (*p1<<*o1);

	for(unsigned int xy=0; xy<nxy; ++xy) {
	    unsigned int ab = 0;
	    for(unsigned int a=0; a<nv1; ++a) {
		const unsigned int p = v1_to_p1[a];
		for(unsigned int b=0; b<nv1; ++b, ++ab) {
		    const unsigned int q = v2_to_p2[b];
		    const unsigned int pq = p*np2 + q;
		    const double elem = Vpq[s].get_element(xy,pq);
		    Vab[s].set_element(xy,ab,elem);
		}
	    }
	}

	for(unsigned int xy=0; xy<nxy; ++xy) {
	    unsigned int ia = 0;
	    for(unsigned int i=0; i<no1; ++i) {
		const unsigned int p = o1_to_p1[i];
		for(unsigned int a=0; a<nv2; ++a, ++ia) {
		    const unsigned int q = v1_to_p1[a];
		    const unsigned int pq = p*np2 + q;
		    const double elem = Vpq[s].get_element(xy,pq);
		    Via[s].set_element(xy,ia,elem);
		}
	    }
	}

	Vpq[s].print("Vpq matrix");
	Vab[s].print("Vab matrix");
#if TEST_V
	RefSCMatrix Vab_test = r12eval->V(spincase2,vir1_act,vir2_act);
	Vab_test.print("Vab matrix (test)");
	(Vab[s] - Vab_test).print("Vab - Vab (test): should be 0");
#endif
	Via[s].print("Via matrix");
#if TEST_V
	RefSCMatrix Via_test = r12eval->V(spincase2,occ1_act,vir2_act);
	Via_test.print("Via matrix (test)");
	(Via[s] - Via_test).print("Via - Via (test): should be 0");
#endif
	Vij[s].print("Vij matrix");
	T2_MP1[s].print("MP1 T2 amplitudes");
    }
    Ref<SCMatrixKit> localkit = Vpq[AlphaBeta].kit();

    // print out MPQC orbitals to compare to Psi orbitals below;
    Ref<MOIndexSpace> orbs_sb_mpqc = r12eval->r12info()->refinfo()->orbs_sb(Alpha);
    orbs_sb_mpqc->coefs().print("MPQC eigenvector");
    orbs_sb_mpqc->evals().print("MPQC eigenvalues");
    Ref<MOIndexSpace> occ_act = r12eval->occ_act(Alpha);
    Ref<MOIndexSpace> vir_act = r12eval->vir_act(Alpha);

    // Psi stores amplitudes in Pitzer (symmetry-blocked) order. Construct such spaces for active occupied and virtual spaces
    RefSCMatrix psi_coefs = reference_->coefs();
    RefDiagSCMatrix psi_evals = reference_->evals();
    Ref<MOIndexSpace> occ_act_sb_psi = occ_act_sb(Alpha);
    Ref<MOIndexSpace> vir_act_sb_psi = vir_act_sb(Alpha);
    
    // map MPQC (energy-ordered) orbitals to Psi (symmetry-ordered) orbitals
    SparseMOIndexMap MPQC2PSI_occ_act;
    SparseMOIndexMap MPQC2PSI_vir_act;
    bool can_map_occ_act = true;
    bool can_map_vir_act = true;
    try {
	SparseMOIndexMap tmpmap_o(sparsemap(*occ_act_sb_psi,*occ_act,1e-9));
	std::swap(MPQC2PSI_occ_act,tmpmap_o);
    }
    catch (CannotConstructMap& e) {
	can_map_occ_act = false;
    }
    try {
	SparseMOIndexMap tmpmap_v(sparsemap(*vir_act_sb_psi,*vir_act,1e-9));
	std::swap(MPQC2PSI_vir_act,tmpmap_v);
    }
    catch (CannotConstructMap& e) {
	can_map_vir_act = false;
    }

    if (use_sparsemap_only_) {
	if (!can_map_occ_act)
	    throw ProgrammingError("PsiCCSD_PT2R12::compute() -- cannot map MPQC occ act. orbitals to Psi occ act. orbitals",__FILE__,__LINE__);
	if (!can_map_vir_act)
	    throw ProgrammingError("PsiCCSD_PT2R12::compute() -- cannot map MPQC vir act. orbitals to Psi vir act. orbitals",__FILE__,__LINE__);
    }

    const Ref<MOIndexSpace>& occ1_act = r12eval->occ_act(Alpha);
    const Ref<MOIndexSpace>& occ2_act = r12eval->occ_act(Beta);
    const Ref<MOIndexSpace>& vir1_act = r12eval->vir_act(Alpha);
    const Ref<MOIndexSpace>& vir2_act = r12eval->vir_act(Beta);
    const bool use_sparsemap = can_map_occ_act && can_map_vir_act;
    if (use_sparsemap) {
	T1[AlphaBeta] = transform_T1(MPQC2PSI_occ_act,MPQC2PSI_vir_act,T1_psi,localkit);
	T2[AlphaBeta] = transform_T2(MPQC2PSI_occ_act,MPQC2PSI_occ_act,MPQC2PSI_vir_act,MPQC2PSI_vir_act,T2_psi,localkit);
	Tau2[AlphaBeta] = transform_T2(MPQC2PSI_occ_act,MPQC2PSI_occ_act,MPQC2PSI_vir_act,MPQC2PSI_vir_act,Tau2_psi,localkit);
	if (test_t2_phases_) {
	    compare_T2(T2[AlphaBeta],T2_MP1[AlphaBeta],
		       occ1_act->rank(),occ2_act->rank(),
		       vir1_act->rank(),vir2_act->rank());
	}
	T1[AlphaBeta].print("CCSD T1 amplitudes from Psi3 (in MPQC orbitals):");
	T2[AlphaBeta].print("CCSD T2 amplitudes from Psi3 (in MPQC orbitals):");
	Tau2[AlphaBeta].print("CCSD Tau2 amplitudes from Psi3 (in MPQC orbitals):");
    }
    else {
	RefSCMatrix MPQC2PSI_tform_oa = transform(*occ_act_sb_psi,*occ_act,localkit);
	RefSCMatrix MPQC2PSI_tform_va = transform(*vir_act_sb_psi,*vir_act,localkit);
	T1[AlphaBeta] = transform_T1(MPQC2PSI_tform_oa,MPQC2PSI_tform_va,T1_psi,localkit);
	T2[AlphaBeta] = transform_T2(MPQC2PSI_tform_oa,MPQC2PSI_tform_oa,MPQC2PSI_tform_va,MPQC2PSI_tform_va,T2_psi,localkit);
	Tau2[AlphaBeta] = transform_T2(MPQC2PSI_tform_oa,MPQC2PSI_tform_oa,MPQC2PSI_tform_va,MPQC2PSI_tform_va,Tau2_psi,localkit);
	if (test_t2_phases_) {
	    compare_T2(T2[AlphaBeta],T2_MP1[AlphaBeta],
		       occ1_act->rank(),occ2_act->rank(),
		       vir1_act->rank(),vir2_act->rank());
	}
	T1[AlphaBeta].print("CCSD T1 amplitudes from Psi3 (in MPQC orbitals, obtained by transform):");
	T2[AlphaBeta].print("CCSD T2 amplitudes from Psi3 (in MPQC orbitals, obtained by transform):");
	Tau2[AlphaBeta].print("CCSD Tau2 amplitudes from Psi3 (in MPQC orbitals, obtained by transform):");
    }

    // Compute Hamiltonian matrix elements
    RefSCMatrix H1_0R[NSpinCases2];
    RefSCMatrix H1_R0[NSpinCases2];
    RefSymmSCMatrix H0_RR[NSpinCases2];
    for(int s=0; s<nspincases2; s++) {
	const SpinCase2 spincase2 = static_cast<SpinCase2>(s);
	const SpinCase1 spin1 = case1(spincase2);
	const SpinCase1 spin2 = case2(spincase2);

	const Ref<MOIndexSpace>& occ1_act = r12eval->occ_act(spin1);
	const Ref<MOIndexSpace>& occ2_act = r12eval->occ_act(spin2);
	const Ref<MOIndexSpace>& vir1_act = r12eval->vir_act(spin1);
	const Ref<MOIndexSpace>& vir2_act = r12eval->vir_act(spin2);

	// H1_0R is just Vij
	H1_0R[s] = Vij[s].clone();
	H1_0R[s].assign(Vij[s]);
	H1_0R[s].print(prepend_spincase(spincase2,"<0|Hb|R>").c_str());
	
	// Vij is also the leading term in H1_R0
	H1_R0[s] = Vij[s].clone();
	H1_R0[s].assign(Vij[s]);

	// if not testing MP1, compute other terms in <R|Hb|0>
	if (!mp2_only_) {

	    // the leading term in <R|(HT)|0> is Tau2.Vab
	    const bool use_tau2 = true;
	    RefSCMatrix HT = Vab[s] * (use_tau2 ? Tau2[s].t() : T2[s].t());
	    H1_R0[s].accumulate(HT);

	    // the next term is T1.Vai
	    {
		// store V_xy_ia as V_xyi_a
		double* tmp = new double[Via[s].rowdim().n() * Via[s].coldim().n()];
		Via[s].convert(tmp);
		RefSCMatrix V_xyi_a = Via[s].kit()->matrix( new SCDimension(Via[s].rowdim().n() * T1[s].rowdim().n()), T1[s].coldim());
		V_xyi_a.assign(tmp);
		delete[] tmp;

		RefSCMatrix V_xyi_J = V_xyi_a * T1[s].t();
		// store V_xyi_J as V_xy_iJ
		tmp = new double[HT.rowdim().n() * HT.coldim().n()];
		V_xyi_J.convert(tmp);
		HT.assign(tmp);
		symmetrize<false>(HT,HT,r12eval->xspace(spin1),occ1_act);
		HT.print("Via.T1");
		delete[] tmp;
		// NOTE: V_{xy}^{ia} T_a^j + V_{xy}^{aj} T_a^i = 2 * symm(V_{xy}^{ia} T_a^j
		HT.scale(2.0);
		H1_R0[s].accumulate(HT);

#if TEST_ViaT1
		// test against V evaluator
		RefSCMatrix vir2_act_coefs = vir2_act->coefs();
		RefSCMatrix to_t1;
		{
		    RefSCDimension rowdim = vir2_act->dim();
		    RefSCDimension coldim = occ2_act->dim();
		    to_t1 = vir2_act_coefs.kit()->matrix(rowdim,coldim);
		    double* tmp = new double[rowdim.n() * coldim.n()];
		    T1[s].t().convert(tmp);
		    to_t1.assign(tmp);
		    delete[] tmp;
		}
		RefSCMatrix Tocc2_act_coefs = vir2_act_coefs * to_t1;
		Ref<MOIndexSpace> Tocc2_act = new MOIndexSpace("i(t1)", "active occupied orbitals (weighted by T1)",
							       occ2_act, Tocc2_act_coefs, occ2_act->basis());
		RefSCMatrix ViJ = r12eval->V(spincase2,occ1_act,Tocc2_act);
		symmetrize<false>(ViJ,ViJ,r12eval->xspace(spin1),occ1_act);
		ViJ.print("Via.T1 computed with R12IntEvalInfo::V()");
		(HT - ViJ).print("Via.T1 - Via.T1 (test): should be zero");
#endif
	    }

	    // ebc term is Tau2.A
	    if (!ebc)
		HT.accumulate_product(A[s],Tau2[s].t());

	    HT.print(prepend_spincase(spincase2,"<R|(H*T)|0>").c_str());
	    
	}

	H1_R0[s].print(prepend_spincase(spincase2,"<R|Hb|0>").c_str());

	// H0_RR is just B
	H0_RR[s] = B[s].clone();
	H0_RR[s].assign(B[s]);
    }
    // Make H1_R0[AlphaAlpha]
    {
	H1_R0[AlphaAlpha] = r12eval->V(AlphaAlpha).clone();
	antisymmetrize(H1_R0[AlphaAlpha],H1_R0[AlphaBeta],r12eval->xspace(Alpha),r12eval->occ_act(Alpha));
    }

    // compute the second-order correction: E2 = - H1_0R . H0_RR^{-1} . H1_R0 = C_MP1 . H1_R0
    const int debug = 0;
    Ref<MP2R12Energy> r12energy = new MP2R12Energy(r12eval,r12eval->r12info()->r12tech()->stdapprox(),debug);
    double E2[NSpinCases2];
    // WARNING only RHF and UHF are considered
    const int num_unique_spincases2 = (reference_->spin_polarized() ? 3 : 2);
    for(int s=0; s<num_unique_spincases2; s++) {
	const SpinCase2 spincase2 = static_cast<SpinCase2>(s);

	RefSCMatrix C_MP1 = r12energy->C(spincase2);
	RefSCMatrix E2_mat = C_MP1.t() * H1_R0[s];
	E2[s] = E2_mat.trace();
    }
    ExEnv::out0() << "E2(AB) = " << E2[AlphaBeta] << endl;
    ExEnv::out0() << "E2(AA) = " << 2.0*E2[AlphaAlpha] << endl;
    ExEnv::out0() << "E2(s) = " << E2[AlphaBeta] - E2[AlphaAlpha] << endl;
    ExEnv::out0() << "E2(t) = " << 3.0*E2[AlphaAlpha] << endl;
}

//////////////////////////////////////////////////////////////////////////

}

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:

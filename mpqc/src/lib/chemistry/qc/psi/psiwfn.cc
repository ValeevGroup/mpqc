
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
#include <chemistry/qc/psi/psiwfn.h>
#include <chemistry/qc/psi/psiinput.timpl.h>

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

  int bytes = keyval->intvalue("memory");
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

static ClassDesc PsiCC_cd(
  typeid(PsiCC),"PsiCC",1,"public PsiCorrWavefunction",
  0, create<PsiCC>, create<PsiCC>);

PsiCC::PsiCC(const Ref<KeyVal>&keyval):
  PsiCorrWavefunction(keyval)
{
}

PsiCC::~PsiCC()
{
}

PsiCC::PsiCC(StateIn&s):
  PsiCorrWavefunction(s)
{
}

void
PsiCC::save_data_state(StateOut&s)
{
  PsiCorrWavefunction::save_data_state(s);
}

const RefSCMatrix&
PsiCC::T(unsigned int rank)
{
    if (rank < 1)
	throw ProgrammingError("PsiCC::T(rank) -- rank must be >= 1", __FILE__, __LINE__);
    if (rank > 2)
	throw FeatureNotImplemented("PsiCC::T(rank) -- rank must be > 2 is not supported by Psi yet", __FILE__, __LINE__);
    if (T_.empty() || T_.size() < rank)
	T_.resize(rank);
    if (T_[rank-1].nonnull())
	return T_[rank-1];

    // For now only handle closed-shell case
    PsiSCF::RefType reftype = reference_->reftype();
    if (reftype != PsiSCF::rhf)
	throw FeatureNotImplemented("PsiCC::T() -- only closed-shell case is implemented",__FILE__,__LINE__);
    // For now only handle C1 case
    //if (nirrep_ != 1)
    //  throw FeatureNotImplemented("PsiCC::T() -- only C1 symmetry can be handled",__FILE__,__LINE__);

    psi::PSIO& psio = exenv()->psio();
    // grab orbital info
    int* doccpi = new int[nirrep_];
    int* mopi = new int[nirrep_];
    psio.open(PSIF_CHKPT,PSIO_OPEN_OLD);
    psio.read_entry(PSIF_CHKPT,"::Closed shells per irrep",reinterpret_cast<char*>(doccpi),nirrep_*sizeof(int));
    psio.read_entry(PSIF_CHKPT,"::MO's per irrep",reinterpret_cast<char*>(mopi),nirrep_*sizeof(int));
    psio.close(PSIF_CHKPT,1);
    std::vector<int> actdoccpi(nirrep_);
    std::vector<int> actuoccpi(nirrep_);
    std::vector<int> actdoccioff(nirrep_);
    std::vector<int> actuoccioff(nirrep_);
    unsigned int ndocc_act = 0;
    unsigned int nuocc_act = 0;
    for(unsigned int irrep=0; irrep<nirrep_; ++irrep) {
	actdoccpi[irrep] = doccpi[irrep] - frozen_docc_[irrep];
	actuoccpi[irrep] = mopi[irrep] - doccpi[irrep] - frozen_uocc_[irrep];
	ndocc_act += actdoccpi[irrep];
	nuocc_act += actuoccpi[irrep];
    }
    actdoccioff[0] = 0;
    actuoccioff[0] = 0;
    for(unsigned int irrep=1; irrep<nirrep_; ++irrep) {
	actdoccioff[irrep] = actdoccioff[irrep-1] + actdoccpi[irrep-1];
	actuoccioff[irrep] = actuoccioff[irrep-1] + actuoccpi[irrep-1];
    }

    // Grab T matrices
    if (rank == 1) {
	// read in the i by a matrix in DPD format
	unsigned int nia_dpd = 0;
	for(unsigned int h=0; h<nirrep_; ++h)  nia_dpd += actdoccpi[h] * actuoccpi[h];
	double* T1 = new double[nia_dpd];
	psio.open(CC_OEI,PSIO_OPEN_OLD);
	psio.read_entry(CC_OEI,"tIA",reinterpret_cast<char*>(T1),nia_dpd*sizeof(double));
	psio.close(CC_OEI,1);

	// form the full matrix
	T_[rank-1] = matrixkit()->matrix(new SCDimension(ndocc_act),new SCDimension(nuocc_act));
	RefSCMatrix& T = T_[rank-1];
	T.assign(0.0);
	unsigned int ia = 0;
	for(unsigned int h=0; h<nirrep_; ++h) {
	    const unsigned int i_offset = actdoccioff[h];
	    const unsigned int a_offset = actuoccioff[h];
	    for(int i=0; i<actdoccpi[h]; ++i)
		for(int a=0; a<actuoccpi[h]; ++a, ++ia)
		    T.set_element(i+i_offset,a+a_offset,T1[ia]);
	}
	delete[] T1;
	T_[rank-1].print("T1 amplitudes");
    }
    else if (rank == 2) {
	// DPD of orbital product spaces
	std::vector<int> ijpi(nirrep_);
	std::vector<int> abpi(nirrep_);
	unsigned int nijab_dpd = 0;
	for(unsigned int h=0; h<nirrep_; ++h) {
	    unsigned int nij = 0;
	    unsigned int nab = 0;
	    for(unsigned int g=0; g<nirrep_; ++g) {
		nij += actdoccpi[g] * actdoccpi[h^g];
		nab += actuoccpi[g] * actuoccpi[h^g];
	    }
	    ijpi[h] = nij;
	    abpi[h] = nab;
	    nijab_dpd += nij*nab;
	}

	// read in T2 in DPD form
	double* T2 = new double[nijab_dpd];
	psio.open(CC_TAMPS,PSIO_OPEN_OLD);
	psio.read_entry(CC_TAMPS,"tIjAb",reinterpret_cast<char*>(T2),nijab_dpd*sizeof(double));
	psio.close(CC_TAMPS,1);

	// convert to the full form
	const unsigned int nij = ndocc_act*ndocc_act;
	const unsigned int nab = nuocc_act*nuocc_act;
	T_[rank-1] = matrixkit()->matrix(new SCDimension(nij),new SCDimension(nab));
	RefSCMatrix& T = T_[rank-1];
	T.assign(0.0);
	unsigned int ijab = 0;
	unsigned int ij_offset = 0;
	unsigned int ab_offset = 0;
	for(unsigned int h=0; h<nirrep_;
		ij_offset+=ijpi[h],
		ab_offset+=abpi[h],
		++h
	    ) {
	    for(unsigned int g=0; g<nirrep_; ++g) {
		unsigned int gh = g^h;
		for(int i=0; i<actdoccpi[g]; ++i) {
		    const unsigned int ii = i + actdoccioff[g];

		    for(int j=0; j<actdoccpi[gh]; ++j) {
			const unsigned int jj = j + actdoccioff[gh];

			const unsigned int ij = ii * ndocc_act + jj;
		
			for(unsigned int f=0; f<nirrep_; ++f) {
			    unsigned int fh = f^h;
			    for(int a=0; a<actuoccpi[f]; ++a) {
				const unsigned int aa = a + actuoccioff[f];

				for(int b=0; b<actuoccpi[fh]; ++b, ++ijab) {
				    const unsigned int bb = b + actuoccioff[fh];
				    const unsigned int ab = aa*nuocc_act + bb;
				
				    T.set_element(ij,ab,T2[ijab]);
				}
			    }
			}
		    }
		}
	    }
	}

	delete[] T2;
	T_[rank-1].print("T2 amplitudes");
    }

    return T_[rank-1];
}

const RefSCMatrix&
PsiCC::Lambda(unsigned int rank)
{
    if (rank < 1)
	throw ProgrammingError("PsiCC::Lambda(rank) -- rank must be >= 1", __FILE__, __LINE__);
    if (rank > 2)
	throw FeatureNotImplemented("PsiCC::Lambda(rank) -- rank must be > 2 is not supported by Psi yet", __FILE__, __LINE__);
    if (Lambda_.empty() || Lambda_.size() < rank)
	Lambda_.resize(rank);
    if (Lambda_[rank-1].nonnull())
	return Lambda_[rank-1];

    throw FeatureNotImplemented("PsiCC::Lambda() -- cannot read Lambda amplitudes yet",__FILE__,__LINE__);
    return Lambda_[rank-1];
}

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
    mbptr12_ =  require_dynamic_cast<MBPT2_R12*>(
	keyval->describedclassvalue("mbpt2r12").pointer(),
	"PsiCCSD_PT2R12::PsiCCSD_PT2R12\n"
	);
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
  input->write_keyword("psi:wfn","ccsd");
  input->close();
}

namespace {
    bool gtzero(double a) { return a > 0.0; }
}

void
PsiCCSD_PT2R12::compute()
{
    const bool test_t2_phases = true;

    // compute CCSD wave function
    PsiWavefunction::compute();

    // grab amplitudes
    RefSCMatrix T1_psi = T(1);
    RefSCMatrix T2_psi = T(2);

    // Compute intermediates
    const double mp2r12_energy = mbptr12_->value();
    Ref<R12IntEval> r12eval = mbptr12_->r12eval();
    RefSCMatrix Vpq[NSpinCases2];
    RefSCMatrix Vab[NSpinCases2];
    RefSCMatrix Vij[NSpinCases2];
    RefSymmSCMatrix B[NSpinCases2];
    RefSymmSCMatrix X[NSpinCases2];
    RefSCMatrix T2_MP1[NSpinCases2];
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

	Vpq[s] = r12eval->V(spincase2,p1,p2);
	Vij[s] = r12eval->V(spincase2);
	X[s] = r12eval->X(spincase2);
	B[s] = r12eval->B(spincase2);
	T2_MP1[s] = r12eval->T2(spincase2);

	// extract Vab from Vpq
	Vab[s] = Vpq[s].kit()->matrix(Vpq[s].rowdim(),T2_MP1[s].coldim());
	const Ref<MOIndexSpace>& v1 = r12eval->vir_act(spin1);
	const Ref<MOIndexSpace>& v2 = r12eval->vir_act(spin2);
	const unsigned int nv1 = v1->rank();
	const unsigned int nv2 = v2->rank();
	const unsigned int nxy = Vpq[s].rowdim().n();

	typedef MOIndexMap OrbMap;
	OrbMap v1_to_p1 (*p1<<*v1);
	OrbMap v2_to_p2 (*p2<<*v2);

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

	Vpq[s].print("Vpq matrix");
	Vab[s].print("Vab matrix");
	Vij[s].print("Vij matrix");
	T2_MP1[s].print("MP1 T2 amplitudes");
    }

    // print out MPQC orbitals to compare to Psi orbitals below;
    Ref<MOIndexSpace> orbs_sb_mpqc = r12eval->r12info()->refinfo()->orbs_sb(Alpha);
    orbs_sb_mpqc->coefs().print("MPQC eigenvector");
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
    try {
	SparseMOIndexMap tmpmap_o(sparsemap(*occ_act_sb_psi,*occ_act,1e-9));
	std::swap(MPQC2PSI_occ_act,tmpmap_o);
	SparseMOIndexMap tmpmap_v(sparsemap(*vir_act_sb_psi,*vir_act,1e-9));
	std::swap(MPQC2PSI_vir_act,tmpmap_v);
    }
    catch (CannotConstructMap& e) {
	throw ProgrammingError("PsiCCSD_PT2R12::compute() -- cannot map MPQC orbitals to Psi orbitals",__FILE__,__LINE__);
    }

    // map Psi T amplitudes to the MPQC orbitals
    RefSCMatrix T1[NSpinCases2];
    RefSCMatrix T2[NSpinCases2];
    T1[AlphaBeta] = Vab[AlphaBeta].kit()->matrix(T1_psi.rowdim(),T1_psi.coldim());  T1[AlphaBeta].assign(0.0);
    T2[AlphaBeta] = Vab[AlphaBeta].kit()->matrix(T2_psi.rowdim(),T2_psi.coldim());  T2[AlphaBeta].assign(0.0);
    for(int s=0; s<nspincases2; s++) {
	const SpinCase2 spincase2 = static_cast<SpinCase2>(s);
	const SpinCase1 spin1 = case1(spincase2);
	const SpinCase1 spin2 = case2(spincase2);

	Ref<R12IntEvalInfo> r12info = r12eval->r12info();
	const Ref<MOIndexSpace>& p1 = r12info->refinfo()->orbs(spin1);
	const Ref<MOIndexSpace>& p2 = r12info->refinfo()->orbs(spin2);
	const unsigned int np1 = p1->rank();
	const unsigned int np2 = p2->rank();

	const Ref<MOIndexSpace>& v1 = r12eval->vir_act(spin1);
	const Ref<MOIndexSpace>& v2 = r12eval->vir_act(spin2);
	const unsigned int nv1 = v1->rank();
	const unsigned int nv2 = v2->rank();

	const Ref<MOIndexSpace>& o1 = r12eval->occ_act(spin1);
	const Ref<MOIndexSpace>& o2 = r12eval->occ_act(spin2);
	const unsigned int no1 = o1->rank();
	const unsigned int no2 = o2->rank();

	typedef MOIndexMap OrbMap;
	OrbMap v1_to_p1 (*p1<<*v1);
	OrbMap v2_to_p2 (*p2<<*v2);
	OrbMap o1_to_p1 (*p1<<*o1);
	OrbMap o2_to_p2 (*p2<<*o2);

	// convert T1 to MPQC orbitals
	for(unsigned int i=0; i<no1; ++i) {
	    const unsigned int ii = MPQC2PSI_occ_act[i].first;
	    const double ii_coef = MPQC2PSI_occ_act[i].second;

	    for(unsigned int a=0; a<nv1; ++a) {
		const unsigned int aa = MPQC2PSI_vir_act[a].first;
		const double aa_coef = MPQC2PSI_vir_act[a].second;

		const double elem = T1_psi.get_element(ii,aa);
		T1[s].set_element(i,a,elem*ii_coef*aa_coef);
	    }
	}

	// convert T2 to MPQC orbitals
	for(unsigned int i=0; i<no1; ++i) {
	    const unsigned int ii = MPQC2PSI_occ_act[i].first;
	    const double ii_coef = MPQC2PSI_occ_act[i].second;

	    for(unsigned int j=0; j<no2; ++j) {
		const unsigned int jj = MPQC2PSI_occ_act[j].first;
		const double jj_coef = MPQC2PSI_occ_act[j].second;

		const unsigned int ij = i*no2+j;
		const unsigned int iijj = ii*no2+jj;

		const double ij_coef = ii_coef * jj_coef;

		for(unsigned int a=0; a<nv1; ++a) {
		    const unsigned int aa = MPQC2PSI_vir_act[a].first;
		    const double aa_coef = MPQC2PSI_vir_act[a].second;

		    const double ija_coef = ij_coef * aa_coef;

		    for(unsigned int b=0; b<nv2; ++b) {
			const unsigned int bb = MPQC2PSI_vir_act[b].first;
			const double bb_coef = MPQC2PSI_vir_act[b].second;

			const unsigned int ab = a*nv2+b;
			const unsigned int aabb = aa*nv2+bb;

			const double t2 = T2_psi.get_element(iijj,aabb);
			const double t2_mpqc = t2 * ija_coef * bb_coef;
			if (test_t2_phases) {
			    const double zero = 1e-8;
			    const double t2_mp1_mpqc = T2_MP1[s].get_element(ij,ab);

			    bool mismapped = false;
			    if (fabs(t2_mp1_mpqc) > zero &&
				fabs(t2_mpqc) < zero)
				mismapped = true;
			    if (fabs(t2_mp1_mpqc) < zero &&
				fabs(t2_mpqc) > zero)
				mismapped = true;
			    if (fabs(t2_mp1_mpqc) > zero &&
				fabs(t2_mpqc) > zero &&
				gtzero(t2_mp1_mpqc) != gtzero(t2_mpqc))
				mismapped = true;
			    if (mismapped)
				throw ProgrammingError("CCSD T2 from Psi3 do not appear to have same structure as MP1 T2 from MPQC",__FILE__,__LINE__);
			}
			T2[s].set_element(ij,ab,t2_mpqc);
		    }
		}
	    }
	}

    }
    T1[AlphaBeta].print("CCSD T1 amplitudes from Psi3 (in MPQC orbitals):");
    T2[AlphaBeta].print("CCSD T2 amplitudes from Psi3 (in MPQC orbitals):");

    // Compute Hamiltonian matrix elements
    RefSCMatrix H1_0R[NSpinCases2];
    RefSCMatrix H1_R0[NSpinCases2];
    RefSCMatrix H0_RR[NSpinCases2];
    for(int s=0; s<nspincases2; s++) {
	const SpinCase2 spincase2 = static_cast<SpinCase2>(s);

	// H1_0R is just Vij
	H1_0R[s] = Vij[s].clone();
	H1_0R[s].assign(Vij[s]);
	
	// Vij is also the leading term in H1_R0
	H1_R0[s] = Vij[s].clone();
	H1_R0[s].assign(Vij[s]);

	// the leading term in <R|(HT)|0> is 1/2 T2.Vab
	RefSCMatrix Vpq_T2 = Vab[s] * (T2[s].t());  Vpq_T2.scale(0.5);
	H1_R0[s].accumulate(Vpq_T2);
	H1_0R[s].print(prepend_spincase(spincase2,"<0|Hb|R>").c_str());
	Vpq_T2[s].print(prepend_spincase(spincase2,"<R|(H*T)|0>").c_str());
	H1_R0[s].print(prepend_spincase(spincase2,"<R|Hb|0>").c_str());

	// H0_RR is just B
	H0_RR[s] = B[s].clone();
	H0_RR[s].assign(B[s]);
    }

    // compute the second-order correction
    
}

//////////////////////////////////////////////////////////////////////////

}

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:

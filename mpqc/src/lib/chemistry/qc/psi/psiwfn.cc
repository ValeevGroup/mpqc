//
// psiwfn.cc
//
// Copyright (C) 2002 Edward Valeev
//
// Author: Edward Valeev <evaleev@vt.edu>
// Maintainer: EV
//
// This file is part of the SC Toolkit.
//
// The SC Toolkit is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The SC Toolkit is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the SC Toolkit; see the file COPYING.LIB.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//
// The U.S. Government is granted a limited license as per AL 91-7.
//

#ifdef __GNUC__
#pragma implementation
#endif

#include <stdexcept>
#include <cmath>

#include <psifiles.h>

#include <util/keyval/keyval.h>
#include <util/misc/formio.h>
#include <util/class/scexception.h>
#include <util/state/stateio.h>
#include <math/scmat/matrix.h>
#include <math/symmetry/pointgrp.h>
#include <chemistry/molecule/molecule.h>
#include <chemistry/qc/mbptr12/mbptr12.h>
#include <chemistry/qc/mbptr12/spin.h>
#include <chemistry/qc/mbptr12/print.h>
#include <chemistry/qc/psi/psiwfn.h>
#include <chemistry/qc/psi/psiinput.timpl.h>

#define TEST_V 0

using namespace std;

namespace sc {
  
  //////////////////////////////////////////////////////////////////////////

  static ClassDesc PsiWavefunction_cd(typeid(PsiWavefunction),
                                      "PsiWavefunction", 2,
                                      "public Wavefunction", 0, 0, 0);
  
  PsiWavefunction::PsiWavefunction(const Ref<KeyVal>&keyval) :
    Wavefunction(keyval) {
    exenv_ << keyval->describedclassvalue("psienv");
    if (exenv_.null()) {
      ExEnv::err0()
          << "PsiWavefunction::PsiWavefunction: no Psi execution environment object (psienv)"
          << endl;
      abort();
    }
    
    nirrep_ = molecule()->point_group()->char_table().order();
    docc_ = read_occ(keyval, "docc", nirrep_);
    socc_ = read_occ(keyval, "socc", nirrep_);
    
    size_t bytes = keyval->sizevalue("memory");
    if (bytes <= 2000000)
      bytes = 2000000;
    int bytes_str_len = (int)ceil(log10((long double)bytes));
    memory_ = new char[bytes_str_len+5];
    sprintf(memory_, "(%ld B)", bytes);
  }
  
  PsiWavefunction::~PsiWavefunction() {
    exenv_->run_psi_module("psiclean");
  }
  
  PsiWavefunction::PsiWavefunction(StateIn&s) :
    SavableState(s), Wavefunction(s) {
    throw std::runtime_error("PsiWavefunction::PsiWavefunction(StateIn&) -- cannot restore state of Psi wave functions");
  }
  
  void PsiWavefunction::save_data_state(StateOut&s) {
    throw std::runtime_error("PsiWavefunction::save_data_state -- cannot save state of Psi wave functions, set savestate = no in your input file");
  }
  
  void PsiWavefunction::print(ostream&o) const {
    Wavefunction::print(o);
    exenv_->print(o);
  }
  
  int PsiWavefunction::debug() const {
    return debug_;
  }
  
  void PsiWavefunction::compute() {
    if (gradient_needed() && !gradient_implemented()) {
      ExEnv::out0()
          << scprintf("Gradient is not implemented for this Psi wavefunction")
          << endl;
      abort();
    }
    double energy_acc = desired_value_accuracy();
    double grad_acc = desired_gradient_accuracy();
    if (energy_acc > 1.0e-6)
      energy_acc = 1.0e-6;
    if (grad_acc > 1.0e-7)
      grad_acc = 1.0e-7;
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
        ExEnv::out0()
            << scprintf("Number of atoms in MPQC and Psi3 do not match")
            << endl;
        abort();
      }
      RefSCVector gradientvec = basis()->matrixkit()->vector(moldim());
      for (int atom=0; atom<natom; atom++) {
        gradientvec[3*atom] = file11->get_grad(0, atom, 0);
        gradientvec[3*atom+1] = file11->get_grad(0, atom, 1);
        gradientvec[3*atom+2] = file11->get_grad(0, atom, 2);
      }
      set_gradient(gradientvec);
      file11->close();
      file11->remove();
    } else {
      double energy = exenv_->chkpt().rd_etot();
      set_energy(energy);
      set_actual_value_accuracy(energy_acc);
    }
  }
  
  RefSymmSCMatrix PsiWavefunction::density() {
    abort();
    return 0;
  }
  
  int PsiWavefunction::nelectron() {
    abort();
    return 0;
  }
  
  void PsiWavefunction::write_basic_input(int conv) {
    const char *dertype = gradient_needed() ? "first" : "none";
    
    Ref<PsiInput> psiinput = get_psi_input();
    psiinput->write_defaults(exenv_, dertype);
    psiinput->write_keyword("psi:memory", memory_);
    psiinput->begin_section("input");
    psiinput->write_keyword("no_reorient", "true");
    psiinput->write_keyword("keep_ref_frame", "true");
    psiinput->write_basis(basis());
    if (basis()->has_pure())
      psiinput->write_keyword("puream", "true");
    psiinput->write_geom(molecule());
    psiinput->end_section();
    psiinput->write_basis_sets(basis());
  }
  
  // Shamelessly borrowed from class SCF
  std::vector<int> PsiWavefunction::read_occ(const Ref<KeyVal> &keyval,
                                             const char *name, int nirrep) {
    std::vector<int> occ;
    if (keyval->exists(name)) {
      if (keyval->count(name) != nirrep) {
        ExEnv::err0() << indent<< "ERROR: PsiWavefunction: have "<< nirrep
            << " irreps but "<< name << " vector is length "
            << keyval->count(name)<< endl;
        abort();
      }
      occ.resize(nirrep);
      for (int i=0; i<nirrep; i++) {
        occ[i] = keyval->intvalue(name, i);
      }
    }
    return occ;
  }
  
  //////////////////////////////////////////////////////////////////////////

  static ClassDesc PsiSCF_cd(typeid(PsiSCF), "PsiSCF", 1,
                             "public PsiWavefunction", 0, 0, 0);
  
  PsiSCF::PsiSCF(const Ref<KeyVal>&keyval) :
    PsiWavefunction(keyval) {
    if (keyval->exists("total_charge"))
      charge_ = keyval->intvalue("total_charge");
    if (keyval->exists("multiplicity")) {
      multp_ = keyval->intvalue("multiplicity");
      if (multp_ < 1) {
        ExEnv::err0() << indent
                      << "ERROR: PsiSCF: valid multiplicity has to be >= 1"<< endl;
        abort();
      }
    }
    else {
      ExEnv::err0() << indent
          << "ERROR: PsiSCF: multiplicity and total_charge need "
          << "to be specified when docc (socc) are missing"<< endl;
      abort();
    }
  }
  
  PsiSCF::~PsiSCF() {
  }
  
  PsiSCF::PsiSCF(StateIn&s) :
    PsiWavefunction(s) {
  }
  
  void PsiSCF::save_data_state(StateOut&s) {
    PsiWavefunction::save_data_state(s);
  }
  
  unsigned int PsiSCF::nmo() {
    int num_mo = exenv()->chkpt().rd_nmo();
    return num_mo;
  }
  
  unsigned int PsiSCF::nocc(SpinCase1 spin) {
    int* doccpi = exenv()->chkpt().rd_clsdpi();
    unsigned int nocc = 0;
    for (unsigned int h=0; h<nirrep_; ++h)
      nocc += doccpi[h];
    psi::Chkpt::free(doccpi);

    if (reftype() != rhf && spin == Alpha) {
      int* soccpi = exenv()->chkpt().rd_openpi();
      for (unsigned int h=0; h<nirrep_; ++h)
        nocc += soccpi[h];
      psi::Chkpt::free(soccpi);
    }

    return nocc;
  }
  
  const RefDiagSCMatrix&PsiSCF::evals(SpinCase1 spin) {
    if (evals_[spin].nonnull())
      return evals_[spin];
    
    PsiSCF::RefType ref = reftype();
    if (ref == rhf && spin == Beta)
      return evals(Alpha);
    
    psi::PSIO& psio = exenv()->psio();
    // grab orbital info
    int num_mo = exenv()->chkpt().rd_nmo();
    int* mopi = exenv()->chkpt().rd_orbspi();
    // get the eigenvalues
    double* E;
    switch (ref) {
      case rhf:
        E = exenv()->chkpt().rd_evals(); break;
      default:
        E = (spin == Alpha) ? exenv()->chkpt().rd_alpha_evals() : exenv()->chkpt().rd_beta_evals();
    }
    
    // convert raw matrices to SCMatrices
    RefSCDimension modim = new SCDimension(num_mo,nirrep_,mopi);
    for (unsigned int h=0; h<nirrep_; ++h)
      modim->blocks()->set_subdim(h, new SCDimension(mopi[h]));
    evals_[spin] = basis_matrixkit()->diagmatrix(modim);
    evals_[spin].assign(E);
    if (debug() >= DefaultPrintThresholds::mostN)
      evals_[spin].print(prepend_spincase(spin,"Psi3 SCF eigenvalues").c_str());
    
    psi::Chkpt::free(E);
    psi::Chkpt::free(mopi);

    return evals_[spin];
  }
  
  const RefSCMatrix&PsiSCF::coefs(SpinCase1 spin) {
    if (coefs_[spin].nonnull())
      return coefs_[spin];
    
    PsiSCF::RefType ref = reftype();
    if (ref == rhf && spin == Beta)
      return coefs(Alpha);
    
    psi::PSIO& psio = exenv()->psio();
    // grab orbital info
    int num_so = exenv()->chkpt().rd_nso();
    int num_mo = exenv()->chkpt().rd_nmo();
    int* mopi = exenv()->chkpt().rd_orbspi();
    int* sopi = exenv()->chkpt().rd_sopi();
    // get MO coefficients in SO basis
    double** C;
    switch (ref) {
      case rhf:
        C = exenv()->chkpt().rd_scf(); break;
      default:
        C = (spin == Alpha) ? exenv()->chkpt().rd_alpha_scf() : exenv()->chkpt().rd_beta_scf();
    }
    // get AO->SO matrix (MPQC AO equiv PSI3 BF)
    double** ao2so = exenv()->chkpt().rd_usotbf();
    
    // convert raw matrices to SCMatrices
    RefSCDimension sodim_nb = new SCDimension(num_so,1);
    sodim_nb->blocks()->set_subdim(0, new SCDimension(num_so));
    RefSCDimension sodim = new SCDimension(num_so,nirrep_,sopi);
    for (unsigned int h=0; h<nirrep_; ++h)
      sodim->blocks()->set_subdim(h, new SCDimension(sopi[h]));
    RefSCDimension modim = new SCDimension(num_mo,nirrep_,mopi);
    for (unsigned int h=0; h<nirrep_; ++h)
      modim->blocks()->set_subdim(h, new SCDimension(mopi[h]));
    RefSCMatrix C_so = basis_matrixkit()->matrix(sodim, modim);
    C_so.assign(C[0]);
    if (debug() >= DefaultPrintThresholds::allN2)
      C_so.print(prepend_spincase(spin,"Psi3 eigenvector in SO basis").c_str());
    RefSCMatrix aotoso = basis_matrixkit()->matrix(sodim, sodim_nb);
    aotoso.assign(ao2so[0]);
    if (debug() >= DefaultPrintThresholds::allN2)
      aotoso.print("Psi3 SO->AO matrix");
    coefs_[spin] = aotoso.t() * C_so;
    if (debug() >= DefaultPrintThresholds::allN2)
      coefs_[spin].print(prepend_spincase(spin,"Psi3 eigenvector in AO basis").c_str());
    
    using psi::Chkpt;
    Chkpt::free(mopi);
    Chkpt::free(sopi);
    Chkpt::free(C);
    Chkpt::free(ao2so);
    
    return coefs_[spin];
  }
  
  //////////////////////////////////////////////////////////////////////////

  static ClassDesc PsiCLHF_cd(typeid(PsiCLHF), "PsiCLHF", 1, "public PsiSCF",
                              0, create<PsiCLHF>, create<PsiCLHF>);
  
  PsiCLHF::PsiCLHF(const Ref<KeyVal>&keyval) :
    PsiSCF(keyval) {
    if (docc_.empty() && multp_ != 1) {
      ExEnv::err0() << indent
          << "ERROR: PsiCLHF: multiplicity should be 1 for CLHF wave function"
          << endl;
      abort();
    }
  }
  
  PsiCLHF::~PsiCLHF() {
  }
  
  PsiCLHF::PsiCLHF(StateIn&s) :
    PsiSCF(s) {
  }
  
  void PsiCLHF::write_basic_input(int convergence) {
    Ref<PsiInput> input = get_psi_input();
    input->write_keyword("psi:reference", "rhf");
    if (!docc_.empty())
      input->write_keyword_array("psi:docc", nirrep_, docc_);
    else {
      input->write_keyword("psi:multp", multp_);
      input->write_keyword("psi:charge", charge_);
      input->write_keyword("psi:reset_occupations", true);
    }
  }
  
  void PsiCLHF::write_input(int convergence) {
    Ref<PsiInput> input = get_psi_input();
    input->open();
    PsiWavefunction::write_basic_input(convergence);
    write_basic_input(convergence);
    input->write_keyword("psi:wfn", "scf");
    input->close();
  }
  
  //////////////////////////////////////////////////////////////////////////

  static ClassDesc PsiHSOSHF_cd(typeid(PsiHSOSHF), "PsiHSOSHF", 1,
                                "public PsiSCF", 0, create<PsiHSOSHF>, create<
                                    PsiHSOSHF>);
  
  PsiHSOSHF::PsiHSOSHF(const Ref<KeyVal>&keyval) :
    PsiSCF(keyval) {
    if ((docc_.empty() || socc_.empty()) && multp_ == 1) {
      ExEnv::err0() << indent
          << "ERROR: PsiHSOSHF: multiplicity should be > 1 for HSOSHF wave function"
          << endl;
      abort();
    }
  }
  
  PsiHSOSHF::~PsiHSOSHF() {
  }
  
  PsiHSOSHF::PsiHSOSHF(StateIn&s) :
    PsiSCF(s) {
  }
  
  void PsiHSOSHF::write_basic_input(int convergence) {
    Ref<PsiInput> input = get_psi_input();
    input->write_keyword("psi:reference", "rohf");
    if (!docc_.empty())
      input->write_keyword_array("psi:docc", nirrep_, docc_);
    if (!socc_.empty())
      input->write_keyword_array("psi:socc", nirrep_, socc_);
    if (docc_.empty() && socc_.empty()) {
      input->write_keyword("psi:multp", multp_);
      input->write_keyword("psi:charge", charge_);
      input->write_keyword("psi:reset_occupations", true);
    }
  }
  
  void PsiHSOSHF::write_input(int convergence) {
    Ref<PsiInput> input = get_psi_input();
    input->open();
    PsiWavefunction::write_basic_input(convergence);
    write_basic_input(convergence);
    input->write_keyword("psi:wfn", "scf");
    input->close();
  }
  
  //////////////////////////////////////////////////////////////////////////

  static ClassDesc PsiUHF_cd(typeid(PsiUHF), "PsiUHF", 1, "public PsiSCF", 0,
                             create<PsiUHF>, create<PsiUHF>);
  
  PsiUHF::PsiUHF(const Ref<KeyVal>&keyval) :
    PsiSCF(keyval) {
  }
  
  PsiUHF::~PsiUHF() {
  }
  
  PsiUHF::PsiUHF(StateIn&s) :
    PsiSCF(s) {
  }
  
  void PsiUHF::write_basic_input(int convergence) {
    Ref<PsiInput> input = get_psi_input();
    input->write_keyword("psi:reference", "uhf");
    if (!docc_.empty())
      input->write_keyword_array("psi:docc", nirrep_, docc_);
    if (!socc_.empty())
      input->write_keyword_array("psi:socc", nirrep_, socc_);
    input->write_keyword("psi:multp", multp_);
    if (docc_.empty() && socc_.empty()) {
      input->write_keyword("psi:charge", charge_);
      input->write_keyword("psi:reset_occupations", true);
    }
  }
  
  void PsiUHF::write_input(int convergence) {
    Ref<PsiInput> input = get_psi_input();
    input->open();
    PsiWavefunction::write_basic_input(convergence);
    write_basic_input(convergence);
    input->write_keyword("psi:wfn", "scf");
    input->close();
  }
  
  //////////////////////////////////////////////////////////////////////////

  static ClassDesc PsiCorrWavefunction_cd(typeid(PsiCorrWavefunction),
                                          "PsiCorrWavefunction", 1,
                                          "public PsiWavefunction", 0, create<
                                              PsiCorrWavefunction>, create<
                                              PsiCorrWavefunction>);
  
  PsiCorrWavefunction::PsiCorrWavefunction(const Ref<KeyVal>&keyval) :
    PsiWavefunction(keyval) {
    frozen_docc_ = read_occ(keyval, "frozen_docc", nirrep_);
    frozen_uocc_ = read_occ(keyval, "frozen_uocc", nirrep_);
    if (frozen_docc_.empty()) {
      frozen_docc_.resize(nirrep_);
      for (unsigned int h=0; h<nirrep_; ++h)
        frozen_docc_[h] = 0;
    }
    if (frozen_uocc_.empty()) {
      frozen_uocc_.resize(nirrep_);
      for (unsigned int h=0; h<nirrep_; ++h)
        frozen_uocc_[h] = 0;
    }
    
    reference_ << keyval->describedclassvalue("reference");
    if (reference_.null()) {
      ExEnv::err0()
          << "PsiCorrWavefunction::PsiCorrWavefunction: no reference wavefunction"
          << endl;
      abort();
    }
  }
  
  PsiCorrWavefunction::~PsiCorrWavefunction() {
  }
  
  PsiCorrWavefunction::PsiCorrWavefunction(StateIn&s) :
    PsiWavefunction(s) {
    reference_ << SavableState::restore_state(s);
    s.get(frozen_docc_);
    s.get(frozen_uocc_);
  }
  
  void PsiCorrWavefunction::save_data_state(StateOut&s) {
    PsiWavefunction::save_data_state(s);
    SavableState::save_state(reference_.pointer(), s);
    s.put(frozen_docc_);
    s.put(frozen_uocc_);
  }
  
  void PsiCorrWavefunction::write_input(int convergence) {
    if (gradient_needed())
      reference_->do_gradient(1);
    else
      reference_->do_gradient(0);
    
    Ref<PsiInput> input = get_psi_input();
    PsiWavefunction::write_basic_input(convergence);
    reference_->write_basic_input(convergence);
    if (!frozen_docc_.empty())
      input->write_keyword_array("psi:frozen_docc", nirrep_, frozen_docc_);
    if (!frozen_uocc_.empty())
      input->write_keyword_array("psi:frozen_uocc", nirrep_, frozen_uocc_);
  }
  
  const Ref<MOIndexSpace>&PsiCorrWavefunction::occ_act_sb(SpinCase1 spin) {
    if (occ_act_sb_[spin].nonnull())
      return occ_act_sb_[spin];
    if (reference_->reftype() == PsiSCF::rhf && spin==Beta)
      return occ_act_sb(Alpha);
    
    const int nmo = reference_->nmo();
    const int nocc = reference_->nocc(spin);
    int nfzc = 0;
    for (unsigned int h=0; h<nirrep_; ++h)
      nfzc += frozen_docc_[h];
    
    const std::string id(spin==Alpha ? "I" : "i");
    occ_act_sb_[spin] = new MOIndexSpace(id,prepend_spincase(spin,"active occupied MOs (Psi3)"),
        reference_->coefs(spin),basis(),integral(),
        reference_->evals(spin),nfzc,nmo-nocc,MOIndexSpace::symmetry);
    
    return occ_act_sb_[spin];
  }
  
  const Ref<MOIndexSpace>&
  PsiCorrWavefunction::vir_act_sb(SpinCase1 spin) {
    if (vir_act_sb_[spin].nonnull())
      return vir_act_sb_[spin];
    if (reference_->reftype() == PsiSCF::rhf && spin==Beta)
      return vir_act_sb(Alpha);

    const int nmo = reference_->nmo();
    const int nocc = reference_->nocc(spin);
    int nfzv = 0;
    for (unsigned int h=0; h<nirrep_; ++h)
      nfzv += frozen_uocc_[h];

    const std::string id(spin==Alpha ? "A" : "a");    
    vir_act_sb_[spin] = new MOIndexSpace(id,prepend_spincase(spin,"active virtual MOs (Psi3)"),
        reference_->coefs(spin),basis(),integral(),
        reference_->evals(spin),nocc,nfzv,MOIndexSpace::symmetry);
    
    return vir_act_sb_[spin];
  }
  
  double
  PsiCorrWavefunction::reference_energy()
  {
    return exenv()->chkpt().rd_eref();
  }

//////////////////////////////////////////////////////////////////////////

}

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:

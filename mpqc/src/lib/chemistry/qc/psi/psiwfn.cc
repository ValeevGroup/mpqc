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
#include <sstream>

#include <psifiles.h>

#include <util/keyval/keyval.h>
#include <util/misc/formio.h>
#include <util/class/scexception.h>
#include <util/state/stateio.h>
#include <math/scmat/matrix.h>
#include <math/scmat/blocked.h>
#include <math/symmetry/pointgrp.h>
#include <chemistry/molecule/molecule.h>
#include <chemistry/qc/basis/petite.h>
#include <chemistry/qc/basis/shellrot.h>
#include <chemistry/qc/mbptr12/mbptr12.h>
#include <chemistry/qc/mbptr12/spin.h>
#include <chemistry/qc/mbptr12/print.h>
#include <chemistry/qc/psi/psiwfn.h>
#include <chemistry/qc/psi/psiinput.timpl.h>

#define TEST_V 0

using namespace std;

namespace {
  template <typename T> T sum(const std::vector<T>& vec) {
    T result(0);
    size_t size = vec.size();
    for(size_t i=0; i<size; ++i)
      result += vec[i];
    return result;
  }
}

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
    
    size_t bytes = keyval->sizevalue("memory");
    if (bytes <= 2000000)
      bytes = 2000000;
    {
      std::ostringstream oss;
      oss << "(" << bytes << " B)";
      const int nchar = oss.str().size();
      memory_ = new char[nchar+1];
      strcpy(memory_,oss.str().c_str());
    }
  }
  
  PsiWavefunction::~PsiWavefunction() {
    exenv_->run_psi_module("psiclean");
  }
  
  PsiWavefunction::PsiWavefunction(StateIn&s) :
    SavableState(s), Wavefunction(s) {
    throw std::runtime_error("PsiWavefunction::PsiWavefunction(StateIn&) -- cannot restore state of Psi wave functions");
  }
  
  Integral::CartesianOrdering
  PsiWavefunction::cartesian_ordering() {
    return Integral::CCACartesianOrdering;
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
    if (debug_ > 1)
      exenv()->get_psi_input()->print();
    exenv()->run_psi();
    
    // read output
    if (gradient_needed()) {
      Ref<PsiFile11> file11 = exenv()->get_psi_file11();
      file11->open();
      
      set_energy(file11->get_energy(0));
      set_actual_value_accuracy(energy_acc);
      
      int natom_mpqc = molecule()->natom();
      int natom = file11->get_natom(0);
      if (natom != natom_mpqc) {
        throw ProgrammingError("Number of atoms in MPQC and Psi3 do not match",__FILE__,__LINE__);
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
      const double energy = exenv()->chkpt().rd_etot();
      set_energy(energy);
      set_actual_value_accuracy(energy_acc);
    }
  }
  
  RefSymmSCMatrix PsiWavefunction::density() {
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
  
  std::vector< std::pair<unsigned int, unsigned int> > PsiWavefunction::shell_map() {
    typedef std::vector< std::pair<unsigned int, unsigned int> > ShellMap;
    const Ref<GaussianBasisSet>& bs = basis();
    // # of MPQC contractions = # of Psi3 shells
    const unsigned int ncontr = exenv()->chkpt().rd_nshell();
    // # of MPQC shells
    const unsigned int nshells = basis()->nshell();
    ShellMap map(ncontr);
    int* snuc = exenv()->chkpt().rd_snuc();
    // ordering of shells on an atom is the same in Psi3 and MPQC
    // but shells si and sj from different atoms (i<j) may be ordered differently (si > sj)
    int first_shell_on_curr_atom = 0;
    int atom_curr = snuc[0] - 1;
    int contr = 0;
    for(unsigned int s=0; s<nshells; ++s) {
      int atom = snuc[contr] - 1;
      if (atom != atom_curr) {
        atom_curr = atom;
        first_shell_on_curr_atom = s;
      }
      const unsigned int shell_mpqc = bs->shell_on_center(atom,s-first_shell_on_curr_atom);
      const unsigned int ncontr_in_shell = bs->shell(shell_mpqc).ncontraction();
      for(unsigned int c=0; c<ncontr_in_shell; ++c, ++contr) {
        map[contr] = make_pair(shell_mpqc,c);
      }
    }
    psi::Chkpt::free(snuc);
    return map;
  }

  std::vector<unsigned int> PsiWavefunction::ao_map() {
    typedef std::pair<unsigned int, unsigned int> ShellContrPair;
    typedef std::vector<ShellContrPair> ShellMap;
    const Ref<GaussianBasisSet>& bs = basis();
    const unsigned int nao = bs->nbasis();
    psi::Chkpt& chkpt = exenv()->chkpt();
    ShellMap smap = shell_map();

    std::vector<unsigned int> map(nao);
    const unsigned int nshells_psi = chkpt.rd_nshell();
    int* first_ao_from_shell_psi = chkpt.rd_puream() ? chkpt.rd_sloc_new() : chkpt.rd_sloc();
    for(unsigned int spsi=0; spsi<nshells_psi; ++spsi) {
      const ShellContrPair& shellcontr = smap[spsi];
      unsigned int smpqc = shellcontr.first;
      unsigned int cmpqc = shellcontr.second;
      const GaussianShell& shell = bs->shell(smpqc);
      unsigned int nbf = shell.nfunction(cmpqc);
      int first_bf_psi = first_ao_from_shell_psi[spsi] - 1;
      int first_bf_mpqc = bs->shell_to_function(smpqc) + shell.contraction_to_function(cmpqc);
      for(unsigned int bf=0; bf<nbf; ++bf) {
        map[first_bf_psi + bf] = first_bf_mpqc + bf;
      }
    }
    psi::Chkpt::free(first_ao_from_shell_psi);
    return map;
  }

  //////////////////////////////////////////////////////////////////////////

  static ClassDesc PsiSCF_cd(typeid(PsiSCF), "PsiSCF", 1,
                             "public PsiWavefunction", 0, 0, 0);
  
  PsiSCF::PsiSCF(const Ref<KeyVal>&keyval) :
    PsiWavefunction(keyval) {
    docc_ = read_occ(keyval, "docc", nirrep_);
    socc_ = read_occ(keyval, "socc", nirrep_);
  }
  
  PsiSCF::~PsiSCF() {
  }
  
  PsiSCF::PsiSCF(StateIn&s) :
    PsiWavefunction(s) {
  }
  
  void PsiSCF::save_data_state(StateOut&s) {
    PsiWavefunction::save_data_state(s);
  }
  
  int PsiSCF::nelectron() {
    return nocc(Alpha) + nocc(Beta);
  }
  
  unsigned int PsiSCF::nmo() {
    int num_mo = exenv()->chkpt().rd_nmo();
    return num_mo;
  }
  
  unsigned int PsiSCF::nocc(SpinCase1 spin) {
    const std::vector<unsigned int>& occpi = this->occpi(spin);
    unsigned int nocc = 0;
    for (unsigned int h=0; h<nirrep_; ++h)
      nocc += occpi[h];
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
    if (debug() >= DefaultPrintThresholds::allN2) {
      aotoso.print("Psi3 SO->AO matrix");
      integral()->petite_list()->sotoao().print("MPQC SO->AO matrix");
      integral()->petite_list()->aotoso().print("MPQC AO->SO matrix");
    }
    coefs_[spin] = aotoso.t() * C_so;
    if (debug() >= DefaultPrintThresholds::allN2)
      coefs_[spin].print(prepend_spincase(spin,"Psi3 eigenvector in AO basis (Psi-ordered)").c_str());

    // shells in Psi3 do not have to follow same order as atoms, as they do in MPQC
    // resort AOs from Psi to MPQC order using Psi3->MPQC shell map
    {
      std::vector<unsigned int> aomap = ao_map();
      const int nao = aomap.size();
      const char* name = "PsiSCF::evecs";
      RefSCMatrix coefs_mpqc = coefs_[spin].clone();
      BlockedSCMatrix* coefs_psi_blkd = require_dynamic_cast<BlockedSCMatrix*>(coefs_[spin].pointer(),name);
      BlockedSCMatrix* coefs_mpqc_blkd = require_dynamic_cast<BlockedSCMatrix*>(coefs_mpqc.pointer(),name);
      for (unsigned int h=0; h<nirrep_; ++h) {
        RefSCMatrix coefs_mpqc_blk = coefs_mpqc_blkd->block(h);
        if (coefs_mpqc_blk.null()) continue;
        RefSCMatrix coefs_psi_blk = coefs_psi_blkd->block(h);
        
        for (unsigned int aopsi=0; aopsi<nao; ++aopsi) {
          RefSCVector row = coefs_psi_blk.get_row(aopsi);
          const unsigned int aompqc = aomap[aopsi];
          coefs_mpqc_blk.assign_row(row, aompqc);
        }
      }
      coefs_[spin] = coefs_mpqc;
    }
    if (debug() >= DefaultPrintThresholds::allN2)
      coefs_[spin].print(prepend_spincase(spin,"Psi3 eigenvector in AO basis (MPQC-ordered)").c_str());

    // Psi3 also uses the symmetry frame for the molecule, whereas MPQC may use a different frame
    // rotate the Psi3 eigenvector from the symmetry frame to the MPQC frame
    {
      SymmetryOperation rr = molecule()->point_group()->symm_frame();  rr.transpose();
      const int nshell = basis()->nshell();
      const char* name = "PsiSCF::evecs";
      BlockedSCMatrix* coefs_blkd = require_dynamic_cast<BlockedSCMatrix*>(coefs_[spin].pointer(),name);
      double* tmpvec_orig = new double[basis()->nbasis()];
      double* tmpvec_tformed = new double[basis()->nbasis()];
      for(int s=0; s<nshell; ++s) {
        const GaussianShell& shell = basis()->shell(s);
        const int ncontr = shell.ncontraction();
        // transform each contraction separately
        for(int c=0; c<ncontr; ++c) {
          const int am = shell.am(c);
          // am=0 functions are invariant to rotations
          if (am == 0) continue;
          ShellRotation sr = integral()->shell_rotation(am,rr,shell.is_pure(c));
          const int nf = sr.dim();
          const int foff = basis()->shell_to_function(s) + shell.contraction_to_function(c);
          // in each block
          for (unsigned int h=0; h<nirrep_; ++h) {
            RefSCMatrix coefs_blk = coefs_blkd->block(h);
            if (coefs_blk.null()) continue;
            const int ncol = coefs_blk.coldim().n();
            // transform each vector
            for(int col=0; col<ncol; ++col) {
              // initialize original vector
              for(int f=0; f<nf; ++f)
                tmpvec_orig[f] = coefs_blk.get_element(f+foff,col);
              // transform
              for(int f=0; f<nf; ++f) {
                double tmp = 0.0;
                for(int g=0; g<nf; ++g) {
                  tmp += sr(f,g) * tmpvec_orig[g];
                }
                tmpvec_tformed[f] = tmp;
              }
              // copy to the original location
              for(int f=0; f<nf; ++f)
                coefs_blk.set_element(f+foff,col,tmpvec_tformed[f]);
            }
          }
        }
      }
      delete[] tmpvec_orig; delete[] tmpvec_tformed;
    }
    if (debug() >= DefaultPrintThresholds::allN2)
      coefs_[spin].print(prepend_spincase(spin,"Psi3 eigenvector in AO basis (MPQC-ordered, in MPQC frame)").c_str());
    
    using psi::Chkpt;
    Chkpt::free(mopi);
    Chkpt::free(sopi);
    Chkpt::free(C);
    Chkpt::free(ao2so);
    
    return coefs_[spin];
  }

  const std::vector<unsigned int>& PsiSCF::occpi(SpinCase1 spin) {
    using psi::Chkpt;
    if (!occpi_[spin].empty())
      return occpi_[spin];
    if (spin == Beta && reftype() == rhf)
      return occpi(Alpha);
    
    occpi_[spin].resize(nirrep_);
    {
      int* doccpi = exenv()->chkpt().rd_clsdpi();
      for (unsigned int h=0; h<nirrep_; ++h)
        occpi_[spin][h] = doccpi[h];
      Chkpt::free(doccpi);
    }
    if (reftype() != rhf && spin == Alpha) {
      int* soccpi = exenv()->chkpt().rd_openpi();
      for (unsigned int h=0; h<nirrep_; ++h)
        occpi_[spin][h] += soccpi[h];
      Chkpt::free(soccpi);
    }
    
    return occpi_[spin];
  }

  const std::vector<unsigned int>& PsiSCF::mopi() {
    using psi::Chkpt;
    if (!mopi_.empty())
      return mopi_;
    
    mopi_.resize(nirrep_);
    {
      int* mopi = exenv()->chkpt().rd_orbspi();
      for (unsigned int h=0; h<nirrep_; ++h)
        mopi_[h] = mopi[h];
      Chkpt::free(mopi);
    }
    
    return mopi_;
  }

  const std::vector<unsigned int>& PsiSCF::uoccpi(SpinCase1 spin) {
    using psi::Chkpt;
    if (!uoccpi_[spin].empty())
      return uoccpi_[spin];
    if (spin == Beta && reftype() == rhf)
      return uoccpi(Alpha);
    
    const std::vector<unsigned int>& mopi = this->mopi();
    const std::vector<unsigned int>& occpi = this->occpi(spin);
    uoccpi_[spin].resize(nirrep_);
    for (unsigned int h=0; h<nirrep_; ++h)
      uoccpi_[spin][h] = mopi[h] - occpi[h];
    
    return uoccpi_[spin];
  }

  void PsiSCF::import_occupations(const Ref<OneBodyWavefunction>& obwfn) {

    RefSCDimension osodim = oso_dimension();
    const bool spin_unrestricted = obwfn->spin_unrestricted();

    // make sure that the number of electrons is the same
    const int nuclear_charge = static_cast<int>(molecule()->nuclear_charge());
    if (charge_ != (nuclear_charge - obwfn->nelectron()) )
      throw InputError("PsiSCF::import_occupations(obwfn) -- number of electrons in obwfn does not match this");
    
    // extract occupations
    std::vector<int> docc_obwfn(nirrep_);
    std::vector<int> socc_obwfn(nirrep_);
    for (int h=0; h<nirrep_; ++h) {
      const int nmo = osodim->blocks()->size(h);
      for (int mo=0; mo<nmo; ++mo) {
        int occ;
        if (spin_unrestricted) {
          const int aocc = (obwfn->alpha_occupation(h, mo) == 1.0) ? 1 : 0;
          const int bocc = (obwfn->beta_occupation(h, mo) == 1.0) ? 1 : 0;
          occ = aocc + bocc;
        } else {
          occ = static_cast<int>(obwfn->occupation(h, mo));
        }
        switch (occ) {
          case 2:
            docc_obwfn[h] += 1;
            break;
          case 1:
            socc_obwfn[h] += 1;
            break;
        }
      }
    }
    
    // if necessary, compare to the existing occupations
    if ( !(docc_.empty() && socc_.empty()) ) {
      if (!docc_.empty() && docc_ != docc_obwfn)
        throw InputError("PsiSCF::import_occupations(obwfn) -- provided docc does not match that from obwfn");
      if (!socc_.empty() && socc_ != socc_obwfn)
        throw InputError("PsiSCF::import_occupations(obwfn) -- provided socc does not match that from obwfn");
    }
    // or just copy
    else {
      docc_ = docc_obwfn;
      socc_ = socc_obwfn;
    }
    
    // lastly, update multp_
    multp_ = 1 + sum(socc_);
    
  }
  
  //////////////////////////////////////////////////////////////////////////

  static ClassDesc PsiCLHF_cd(typeid(PsiCLHF), "PsiCLHF", 1, "public PsiSCF",
                              0, create<PsiCLHF>, create<PsiCLHF>);
  
  PsiCLHF::PsiCLHF(const Ref<KeyVal>&keyval) :
    PsiSCF(keyval) {
    
    multp_ = 1;
    const int nuclear_charge = static_cast<int>(molecule()->nuclear_charge());
    if (docc_.empty()) {
      charge_ = keyval->intvalue("total_charge",KeyValValueint(0));
      if ( (nuclear_charge + charge_) % 2 != 0) {
        throw InputError("PsiCLHF::PsiCLHF -- odd number of electrons, charge keyword is not given or incorrect");
      }

      // try guess wavefunction
      guess_wfn_ << keyval->describedclassvalue("guess_wavefunction");
      if (guess_wfn_.nonnull()) {
        import_occupations(guess_wfn_);
      }
    }
    else {
      const int nelectron = sum(docc_) * 2;
      charge_ = nuclear_charge - nelectron;
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
    input->write_keyword("psi:multp", multp_);
    input->write_keyword("psi:charge", charge_);
    if (!docc_.empty())
      input->write_keyword_array("psi:docc", nirrep_, docc_);
    else {
      input->write_keyword("psi:reset_occupations", true);
      input->write_keyword("psi:hcore_guess", "new");
    }
    input->write_keyword("scf:maxiter", maxiter);
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
    
    if ( (docc_.empty() && !socc_.empty()) ||
         (!docc_.empty() && socc_.empty()) ) {
      throw InputError("PsiHSOSHF::PsiHSOSHF -- must give both docc and socc keywords, or neither");
    }
    const int nuclear_charge = static_cast<int>(molecule()->nuclear_charge());
    if (docc_.empty() || socc_.empty()) {
      charge_ = keyval->intvalue("total_charge",KeyValValueint(0));
      // lowest possible multp is 1 (when # of electrons is even) or 2 (odd)
      const int lowest_possible_multp = 1 + (nuclear_charge - charge_)%2;
      multp_ = keyval->intvalue("multiplicity",KeyValValueint(lowest_possible_multp));
      if ((nuclear_charge + charge_)%2 != (multp_ - 1)%2) {
        throw InputError("PsiHSOSHF::PsiHSOSHF -- inconsistent total_charge and multiplicty");
      }
      // try guess wavefunction
      guess_wfn_ << keyval->describedclassvalue("guess_wavefunction");
      if (guess_wfn_.nonnull()) {
        import_occupations(guess_wfn_);
      }
    }
    else {
      const int nelectron = sum(docc_) * 2 + sum(socc_);
      charge_ = nuclear_charge - nelectron;
      multp_ = sum(socc_) + 1;
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
    input->write_keyword("psi:multp", multp_);
    input->write_keyword("psi:charge", charge_);
    if (docc_.empty() && socc_.empty()) {
      input->write_keyword("psi:reset_occupations", true);
      input->write_keyword("psi:hcore_guess", "new");
    }
    input->write_keyword("scf:maxiter", maxiter);
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

    if ( (docc_.empty() && !socc_.empty()) ||
         (!docc_.empty() && socc_.empty()) ) {
      throw InputError("PsiUHF::PsiUHF -- must give both docc and socc keywords, or neither");
    }
    const int nuclear_charge = static_cast<int>(molecule()->nuclear_charge());
    if (docc_.empty() || socc_.empty()) {
      charge_ = keyval->intvalue("total_charge",KeyValValueint(0));
      // lowest possible multp is 1 (when # of electrons is even) or 2 (odd)
      const int lowest_possible_multp = 1 + (nuclear_charge - charge_)%2;
      multp_ = keyval->intvalue("multiplicity",KeyValValueint(lowest_possible_multp));
      if ((nuclear_charge + charge_)%2 != (multp_ - 1)%2) {
        throw InputError("PsiUHF::PsiUHF -- inconsistent total_charge and multiplicity");
      }
      // try guess wavefunction
      guess_wfn_ << keyval->describedclassvalue("guess_wavefunction");
      if (guess_wfn_.nonnull()) {
        import_occupations(guess_wfn_);
      }
    }
    else {
      const int nelectron = sum(docc_) * 2 + sum(socc_);
      charge_ = nuclear_charge - nelectron;
      multp_ = sum(socc_) + 1;
    }
  
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
    input->write_keyword("psi:charge", charge_);
    if (docc_.empty() && socc_.empty()) {
      input->write_keyword("psi:reset_occupations", true);
      input->write_keyword("psi:hcore_guess", "new");
    }
    input->write_keyword("scf:maxiter", maxiter);
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
                                          "public PsiWavefunction", 0, create<PsiCorrWavefunction>,
                                          create<PsiCorrWavefunction>);
  
  PsiCorrWavefunction::PsiCorrWavefunction(const Ref<KeyVal>&keyval) :
    PsiWavefunction(keyval) {
    
    std::string nfzc_str = keyval->stringvalue("nfzc",KeyValValuestring("0"));
    if (nfzc_str == "auto")
      nfzc_ = molecule()->n_core_electrons()/2;
    else if (nfzc_str == "no" || nfzc_str == "false")
      nfzc_ = 0;
    else
      nfzc_ = atoi(nfzc_str.c_str());
    std::string nfzv_str = keyval->stringvalue("nfzv",KeyValValuestring("0"));
    if (nfzv_str == "no" || nfzv_str == "false")
      nfzv_ = 0;
    else
      nfzv_ = atoi(nfzv_str.c_str());
    
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
    int nfzc; s.get(nfzc); nfzc_ = static_cast<unsigned int>(nfzc);
    int nfzv; s.get(nfzv); nfzv_ = static_cast<unsigned int>(nfzv);
  }
  
  void PsiCorrWavefunction::save_data_state(StateOut&s) {
    PsiWavefunction::save_data_state(s);
    SavableState::save_state(reference_.pointer(), s);
    s.put(static_cast<int>(nfzc_));
    s.put(static_cast<int>(nfzv_));
  }

  void PsiCorrWavefunction::set_desired_value_accuracy(double acc) {
    Function::set_desired_value_accuracy(acc);
    reference()->set_desired_value_accuracy( acc /
                                             valacc_to_refacc() );
  }
  
  void PsiCorrWavefunction::write_input(int convergence) {
    if (gradient_needed())
      reference_->do_gradient(1);
    else
      reference_->do_gradient(0);
    
    Ref<PsiInput> input = get_psi_input();
    PsiWavefunction::write_basic_input(convergence);
    reference_->write_basic_input(convergence);
    input->write_keyword("psi:freeze_core", static_cast<int>(nfzc_));
    input->write_keyword("psi:freeze_virt", static_cast<int>(nfzv_));
  }
  
  int PsiCorrWavefunction::nelectron() {
    return reference()->nelectron();
  }
  
  const Ref<MOIndexSpace>&PsiCorrWavefunction::occ_act_sb(SpinCase1 spin) {
    if (occ_act_sb_[spin].nonnull())
      return occ_act_sb_[spin];
    if (reference_->reftype() == PsiSCF::rhf && spin==Beta)
      return occ_act_sb(Alpha);
    
    const int nmo = reference_->nmo();
    const int nocc = reference_->nocc(spin);
    
    const std::string id(spin==Alpha ? "I" : "i");
    occ_act_sb_[spin] = new MOIndexSpace(id,prepend_spincase(spin,"active occupied MOs (Psi3)"),
        reference_->coefs(spin),basis(),integral(),
        reference_->evals(spin),nfzc_,nmo-nocc,MOIndexSpace::symmetry);
    
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

    const std::string id(spin==Alpha ? "A" : "a");    
    vir_act_sb_[spin] = new MOIndexSpace(id,prepend_spincase(spin,"active virtual MOs (Psi3)"),
        reference_->coefs(spin),basis(),integral(),
        reference_->evals(spin),nocc,nfzv_,MOIndexSpace::symmetry);
    
    return vir_act_sb_[spin];
  }

  const std::vector<unsigned int>&
  PsiCorrWavefunction::frozen_docc() {
    if (frozen_docc_.empty()) {
      frozen_docc_.resize(nirrep_);
      int* frzcpi = exenv()->chkpt().rd_frzcpi();
      std::copy(frzcpi,frzcpi+nirrep_,frozen_docc_.begin());
      psi::Chkpt::free(frzcpi);
    }
    return frozen_docc_;
  }
  
  const std::vector<unsigned int>&
  PsiCorrWavefunction::frozen_uocc() {
    if (frozen_uocc_.empty()) {
      frozen_uocc_.resize(nirrep_);
      int* frzvpi = exenv()->chkpt().rd_frzvpi();
      std::copy(frzvpi,frzvpi+nirrep_,frozen_uocc_.begin());
      psi::Chkpt::free(frzvpi);
    }
    return frozen_uocc_;
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

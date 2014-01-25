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

#include <stdexcept>
#include <cmath>
#include <sstream>
#include <numeric>
#include <functional>
#include <algorithm>
#include <cassert>

#include <psifiles.h>

#include <util/keyval/keyval.h>
#include <util/misc/formio.h>
#include <util/misc/scexception.h>
#include <util/state/stateio.h>
#include <math/scmat/local.h>
#include <math/scmat/util.h>
#include <math/scmat/matrix.h>
#include <math/scmat/blocked.h>
#include <math/symmetry/pointgrp.h>
#include <chemistry/molecule/molecule.h>
#include <chemistry/qc/basis/petite.h>
#include <chemistry/qc/basis/shellrot.h>
#include <chemistry/qc/mbptr12/mbptr12.h>
#include <chemistry/qc/lcao/wfnworld.h>
#include <chemistry/qc/wfn/spin.h>
#include <util/misc/print.h>
#include <chemistry/qc/mbptr12/r12int_eval.h>
#include <chemistry/qc/mbptr12/fixedcoefficient.h>
#include <chemistry/qc/wfn/orbitalspace_utils.h>
#include <chemistry/qc/lcao/fockbuilder.h>

#include <chemistry/qc/mbptr12/creator.h>
#include <chemistry/qc/mbptr12/container.h>
#include <chemistry/qc/mbptr12/compute_tbint_tensor.h>

#include <chemistry/qc/psi/psiwfn.h>
#include <chemistry/qc/psi/psiqtorder.h>
#include <math/mmisc/pairiter.h>
#include <chemistry/qc/scf/clscf.h>

#include <libiwl/iwl.h>
#include <libciomr/libciomr.h>

#define TEST_V 0

/// for indexing of density matrices.
//#define triang_half_INDEX(i,j) ((i>j) ? (ioff[(i)]+(j)) : (ioff[(j)]+(i)))
//#define ordinary_INDEX(i,j,coldim) ((i)*(coldim)+(j))

using namespace std;

namespace sc {

  namespace {

    template <typename T>
    T triang_half_INDEX_ordered(T i, T j) {
      return i*(i+1)/2+j;
    }

    template <typename T>
    T triang_half_INDEX(T i, T j) {
      return (i>j) ? triang_half_INDEX_ordered(i,j) : triang_half_INDEX_ordered(j,i);
    }

    template <typename T>
    T ordinary_INDEX(T i, T j, T coldim) {
      return i*coldim+j;
    }

    /// tpdm_index and init_ioff: for indexing of density matrices.
    int tpdm_index(int i, int j, int k, int l,int coldim){
      //int ind_half1=triang_half_INDEX(i,j);
      //int ind_half2=triang_half_INDEX(k,l);
      int ind_half1=ordinary_INDEX(i,j,coldim);
      int ind_half2=ordinary_INDEX(k,l,coldim);
      return(triang_half_INDEX(ind_half1,ind_half2));
    }

    void vector_to_symmmatrix(RefSymmSCMatrix &matrix, const RefSCVector &vector) {
      int dim = matrix.dim().n();
      for(int i=0; i<dim; i++){
        for(int j=0; j<=i; j++) {
          matrix.set_element(i,j,vector.get_element(triang_half_INDEX(i,j)));
        }
      }
    }

    void symmmatrix_to_vector(RefSCVector &vector, const RefSymmSCMatrix &matrix) {
      int dim = matrix.dim().n();
      for(int i=0; i<dim; i++){
        for(int j=0; j<=i; j++) {
          vector.set_element(triang_half_INDEX(i,j),matrix.get_element(i,j));
        }
      }
    }

    void vector_to_matrix(RefSCMatrix &matrix, const RefSCVector &vector) {
      int dim1 = matrix.rowdim().n();
      int dim2 = matrix.coldim().n();
      for(int i=0; i<dim1; i++) {
        for(int j=0; j<dim2; j++) {
          matrix.set_element(i,j,vector.get_element(ordinary_INDEX(i,j,dim2)));
        }
      }
    }

    template <typename T>
    T lowerupper_index(T p, T q) {
      MPQC_ASSERT(p!=q);
      T result = (p>q) ? p*(p-1)/2+q : q*(q-1)/2+p;
      return result;
    }

    void vector_to_matrix(RefSCMatrix &matrix,const RefSCVector &vector,const SpinCase2 &pairspin) {
      int dim1 = matrix.rowdim().n();
      int dim2 = matrix.coldim().n();
      if(pairspin==AlphaBeta) {
        for(int i=0; i<dim1; i++) {
          for(int j=0; j<dim2; j++) {
            matrix.set_element(i,j,vector.get_element(ordinary_INDEX(i,j,dim2)));
          }
        }
      }
      else {  // pairspin==AlphaAlpha || pairspin==BetaBeta
        matrix->assign(0.0);
        for(int i=0; i<dim1; i++) {
          for(int j=0; j<i; j++) {
            const double value = vector.get_element(lowerupper_index(i,j));
            matrix.set_element(i,j,value);
            matrix.set_element(j,i,-value);
          }
        }
      }
    }

    void matrix_to_vector(RefSCVector &vector, const RefSCMatrix &matrix) {
      int dim1 = matrix.rowdim().n();
      int dim2 = matrix.coldim().n();
      for(int i=0; i<dim1; i++) {
        for(int j=0; j<dim2; j++) {
          vector.set_element(ordinary_INDEX(i,j,dim2),matrix.get_element(i,j));
        }
      }
    }

    void matrix_to_vector(RefSCVector &vector, const RefSCMatrix &matrix,const SpinCase2 &pairspin) {
      int dim1 = matrix.rowdim().n();
      int dim2 = matrix.coldim().n();
      if(pairspin==AlphaBeta) {
        for(int i=0; i<dim1; i++) {
          for(int j=0; j<dim2; j++) {
            vector.set_element(ordinary_INDEX(i,j,dim2),matrix.get_element(i,j));
          }
        }
      }
      else {  // pairspin==AlphaAlpha || pairspin==BetaBeta
        for(int i=0; i<dim1; i++) {
          for(int j=0; j<i; j++) {
            vector.set_element(lowerupper_index(i,j),matrix.get_element(i,j));
          }
        }
      }
    }

    double indexsizeorder_sign(int p,int q) {
      if(p>q) {
        return(1.0);
      }
      else if(q>p) {
        return(-1.0);
      }
      else {
        return(0.0);
      }
    }

  }

  //////////////////////////////////////////////////////////////////////////

  static ClassDesc PsiWavefunction_cd(typeid(PsiWavefunction),
                                      "PsiWavefunction", 2,
                                      "public Wavefunction", 0, 0, 0);

  PsiWavefunction::PsiWavefunction(const Ref<KeyVal>&keyval) :
    Wavefunction(keyval) {

    prerequisite_ << keyval->describedclassvalue("prerequisite");

    exenv_ << keyval->describedclassvalue("psienv");
    if (exenv_ == 0) {
      exenv_ = PsiExEnv::get_default_instance();
    }

    nirrep_ = molecule()->point_group()->char_table().order();

    size_t bytes = keyval->sizevalue("memory");
    if (bytes <= 32000000)
      bytes = 32000000;
    {
      std::ostringstream oss;
      oss << "(" << bytes << " B)";
      const int nchar = oss.str().size();
      memory_str_ = new char[nchar+1];
      strcpy(memory_str_,oss.str().c_str());
    }
    memory_ = bytes;

    compute_1rdm_ = keyval->booleanvalue("compute_1rdm", KeyValValueboolean(false));
    dertype_ = keyval->pcharvalue("dertype", KeyValValuepchar("none"));
  }

  PsiWavefunction::~PsiWavefunction() {
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

    this->print_desired_accuracy();

    const bool have_prerequisite = prerequisite_;
    if (have_prerequisite) {
      prerequisite_->compute();
      this->exenv()->run_psiclean(false);  // do partial cleanup
      exenv()->keep_output();
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
    const bool keep_initial_psi_state = have_prerequisite;
    const bool skip_input = keep_initial_psi_state;
    exenv()->run_psi(skip_input);

    // operation of PsiWavefunction assumes that its number of orbitals
    // exactly matches that of Wavefunction
    {
      const int num_mo_psi = exenv()->chkpt().rd_nmo();
      const int num_mo_mpqc = this->oso_dimension().n();
      if (num_mo_psi != num_mo_mpqc) {
        std::ostringstream oss;
        oss << "number of orbitals in Psi and MPQC do not match" << std::endl;
        oss
            << "this is likely due to the different linear dependency thresholds/algorithms"
            << std::endl;
        oss << "nmo (Psi)  = " << num_mo_psi << std::endl;
        oss << "nmo (MPQC) = " << num_mo_mpqc << std::endl;
        throw FeatureNotImplemented(oss.str().c_str(), __FILE__, __LINE__);
      }
    }

    // compare nuclear repulsion energies from Psi3 and MPQC -- warn, if do not agree
    {
      const double nucrep_energy_mpqc = this->molecule()->nuclear_repulsion_energy();
      const double nucrep_energy_psi3 = this->nuclear_repulsion_energy();
      const double enuc_diff = std::fabs(nucrep_energy_mpqc - nucrep_energy_psi3);
      const double enuc_tol = this->desired_value_accuracy() * 0.1;
      if (enuc_diff > enuc_tol) {
        ExEnv::out0() << "WARNING: PsiWavefunction::compute -- MPQC and Psi3 nuclear repulsion energies differ by "
                      << enuc_diff << " but expected less than " << enuc_tol << std::endl;
        ExEnv::out0() << "Psi3 E(nucrep) = " << scprintf("%25.15lf", nucrep_energy_psi3)
            << "  MPQC E(nucrep) = " <<  scprintf("%25.15lf", nucrep_energy_mpqc) << std::endl;
      }
    }

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

      // MPQC feeds Psi3 geometry in "native" symmetry frame -- convert Psi3 gradients to the "reference" frame
      RefSCVector gradientvec = basis()->matrixkit()->vector(moldim());
      double grad_sf[3], grad_rf[3];
      const SymmetryOperation rf_to_sf = molecule()->point_group()->symm_frame();
      for (int atom=0; atom<natom; atom++) {
        grad_sf[0] = file11->get_grad(0, atom, 0);
        grad_sf[1] = file11->get_grad(0, atom, 1);
        grad_sf[2] = file11->get_grad(0, atom, 2);
        for(int i=0; i<3; ++i) {
          grad_rf[i] = 0.0;
          for(int j=0; j<3; ++j)
            grad_rf[i] += grad_sf[j] * rf_to_sf[i][j];
        }
        const int offset = 3*atom;
        for(int i=0; i<3; ++i)
          gradientvec[offset+i] = grad_rf[i];
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
    throw FeatureNotImplemented("densities not implemented for Psi-based Wavefunctions",
                                __FILE__, __LINE__, this->class_desc());
    return 0;
  }

  void PsiWavefunction::write_basic_input(int conv) {
    const char *dertype = gradient_needed() ? "first" : dertype_;

    Ref<PsiInput> psiinput = get_psi_input();
    psiinput->write_defaults(exenv_, dertype);
    if (compute_1rdm_) psiinput->write_keyword("psi:jobtype", "oeprop");
    psiinput->write_keyword("psi:memory", memory_str_);
    psiinput->write_keyword("psi:convergence", conv);
    psiinput->write_keyword("psi:lindep_cutoff", this->lindep_tol());
    psiinput->write_keyword("cints:cutoff", conv+5);
    psiinput->begin_section("input");
    psiinput->write_keyword("no_reorient", "true");
    psiinput->write_keyword("keep_ref_frame", "true");
    psiinput->write_basis(basis());
    if (basis()->has_pure())
      psiinput->write_keyword("puream", "true");
    psiinput->write_geom(molecule());
    psiinput->end_section();
    if (electric_field() //&& !compute_1rdm_
        ) {
      double E[3];  for(int xyz=0; xyz<3; ++xyz) { E[xyz] = electric_field().get_element(xyz); }
      psiinput->write_keyword_array("psi:efield", 3, E);
    }
    psiinput->write_basis_sets(basis());
  }

  // Shamelessly borrowed from class SCF
  std::vector<unsigned int> PsiWavefunction::read_occ(const Ref<KeyVal> &keyval,
                                                      const char *name,
                                                      size_t nirrep) {
    std::vector<unsigned int> occ;
    if (keyval->exists(name)) {
      if (keyval->count(name) != nirrep) {
        ExEnv::err0() << indent<< "ERROR: PsiWavefunction: have "<< nirrep
            << " irreps but "<< name << " vector is length "
            << keyval->count(name)<< endl;
        abort();
      }
      occ.resize(nirrep);
      for (int i=0; i<nirrep; i++) {
        const int occ_i = keyval->intvalue(name, i);
        if (occ_i < 0) {
          std::ostringstream oss;
          oss << "negative occupancy found in array " << name;
          throw InputError(oss.str().c_str(),__FILE__,__LINE__);
        }
        else occ[i] = occ_i;
      }
    }
    return occ;
  }

#if 0
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
#endif

  double
  PsiWavefunction::nuclear_repulsion_energy()
  {
    return exenv()->chkpt().rd_enuc();
  }

  void
  PsiWavefunction::obsolete() {
    if (prerequisite_) prerequisite_->obsolete();
    this->exenv()->run_psiclean();  // do full cleanup
    Wavefunction::obsolete();
  }

  void
  PsiWavefunction::symmetry_changed() {
    nirrep_ = molecule()->point_group()->char_table().order();
    Wavefunction::symmetry_changed();
    if (prerequisite_) prerequisite_->symmetry_changed();
  }

  //////////////////////////////////////////////////////////////////////////

  static ClassDesc PsiSCF_cd(typeid(PsiSCF), "PsiSCF", 1,
                             "public PsiWavefunction", 0, 0, 0);

  PsiSCF::PsiSCF(const Ref<KeyVal>&keyval) :
    PsiWavefunction(keyval) {
    docc_ = read_occ(keyval, "docc", nirrep_);
    socc_ = read_occ(keyval, "socc", nirrep_);
    maxiter_ = keyval->intvalue("maxiter",KeyValValueint(default_maxiter));
    diisdamp_ = keyval->doublevalue("diisdamp",KeyValValuefloat(0.00));
    levelshift_ = keyval->doublevalue("levelshift",KeyValValuedouble(1.0));
    diis_ = keyval->booleanvalue("diis",KeyValValueboolean((int)true));

    guess_wfn_ << keyval->describedclassvalue("guess_wavefunction");
    if (guess_wfn_) {
      if (guess_wfn_->desired_value_accuracy_set_to_default())
        guess_wfn_->set_desired_value_accuracy( this->desired_value_accuracy() * this->guess_acc_ratio() );
      // get energy to make sure that it's computed.
      const double energy = guess_wfn_->value();
    }

  }

  PsiSCF::~PsiSCF() {
  }

  PsiSCF::PsiSCF(StateIn&s) :
    PsiWavefunction(s) {
    s.get(maxiter_);
  }

  void PsiSCF::save_data_state(StateOut&s) {
    PsiWavefunction::save_data_state(s);
    s.put(maxiter_);
  }

  double PsiSCF::magnetic_moment() const {
    return (const_cast<PsiSCF*>(this)->nocc(Alpha) -
            const_cast<PsiSCF*>(this)->nocc(Beta));
  }

  int PsiSCF::nelectron() {
    return nocc(Alpha) + nocc(Beta);
  }

  unsigned int PsiSCF::nmo() {
    const int num_mo = exenv()->chkpt().rd_nmo();
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
    PsiSCF::RefType ref = reftype();

    SpinCase1 spin_to_seek = spin;
    if (ref == rhf && (spin == AnySpinCase1 || spin == Beta))  // for closed-shell all orbitals are the same
      spin_to_seek = Alpha;
    if (ref == hsoshf && spin == AnySpinCase1)  // for restricted open-shell and AnySpinCase1 will assume spin-restricted case but write to Alpha
      spin_to_seek = Alpha;
    if (ref == uhf && spin == AnySpinCase1) // doesn't make sense to ask for any-spin for UHF
      ProgrammingError("asked for any spin orbitals but the wavefunction is spin-unrestricted");

    if (evals_[spin_to_seek])
      return evals_[spin_to_seek];

    PsiChkpt chkpt(exenv(), integral(), debug());
    const bool seek_spin_restricted = (spin == AnySpinCase1 || ref == rhf);   // assume spin restricted if RHF or asked for AnySpinCase1
    evals_[spin_to_seek] = chkpt.evals(spin_to_seek, seek_spin_restricted);

#if 0
    // grab orbital info
    int num_mo = exenv()->chkpt().rd_nmo();
    int* mopi = exenv()->chkpt().rd_orbspi();
    // get the eigenvalues
    double* E;
    switch (ref) {
      case rhf:
        E = exenv()->chkpt().rd_evals(); break;
      case hsoshf:
        E = (spin == Alpha) ? exenv()->chkpt().rd_alpha_evals() : exenv()->chkpt().rd_beta_evals();
        if (E == 0)
          E = exenv()->chkpt().rd_evals();
        break;
      case uhf:
        E = (spin == Alpha) ? exenv()->chkpt().rd_alpha_evals() : exenv()->chkpt().rd_beta_evals();
        break;
      default:
        throw ProgrammingError("unknown reference",__FILE__,__LINE__);
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
#endif
    return evals_[spin_to_seek];
  }

  const RefSCMatrix&PsiSCF::coefs(SpinCase1 spin) {
    PsiSCF::RefType ref = reftype();

    SpinCase1 spin_to_seek = spin;
    if (ref == rhf && (spin == AnySpinCase1 || spin == Beta))  // for closed-shell all orbitals are the same
      spin_to_seek = Alpha;
    if (ref == hsoshf && spin == AnySpinCase1)  // for restricted open-shell and AnySpinCase1 will assume spin-restricted case but write to Alpha
      spin_to_seek = Alpha;
    if (ref == uhf && spin == AnySpinCase1) // doesn't make sense to ask for any-spin for UHF
      ProgrammingError("asked for any spin orbitals but the wavefunction is spin-unrestricted");

    if (coefs_[spin_to_seek])
      return coefs_[spin_to_seek];

    PsiChkpt chkpt(exenv(), integral(), debug());
    const bool seek_spin_restricted = (spin == AnySpinCase1 || ref == rhf);
    coefs_[spin_to_seek] = chkpt.coefs(spin_to_seek, seek_spin_restricted);
#if 0
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
      case hsoshf:
        C = (spin == Alpha) ? exenv()->chkpt().rd_alpha_scf() : exenv()->chkpt().rd_beta_scf();
        if (C == 0) // no semicanonical orbitals!
          C = exenv()->chkpt().rd_scf();
        break;
      case uhf:
        C = (spin == Alpha) ? exenv()->chkpt().rd_alpha_scf() : exenv()->chkpt().rd_beta_scf();
        break;
      default:
        throw ProgrammingError("unknown reference",__FILE__,__LINE__);
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
    Ref<PetiteList> plist = integral()->petite_list();
    if (debug() >= DefaultPrintThresholds::allN2) {
      aotoso.print("Psi3 SO->AO matrix");
      plist->sotoao().print("MPQC SO->AO matrix");
      plist->aotoso().print("MPQC AO->SO matrix");
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
        if (coefs_mpqc_blk == 0) continue;
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
            if (coefs_blk == 0) continue;
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

    // lastly, change the dimensions to match those used by SCF classes (AO dimension much have subdimension blocked by shells)
    {
      RefSCMatrix coefs_redim = coefs_[spin].kit()->matrix(plist->AO_basisdim(), coefs_[spin].coldim());
      RefSCMatrix coefs = coefs_[spin];
      const int nrow = coefs.nrow();
      const int ncol = coefs.ncol();
      for(int r=0; r<nrow; ++r) {
        for(int c=0; c<ncol; ++c) {
          coefs_redim.set_element(r, c, coefs.get_element(r,c) );
        }
      }
      coefs_[spin] = coefs_redim;
    }

    if (debug() >= DefaultPrintThresholds::allN2)
      coefs_[spin].print(prepend_spincase(spin,"Psi3 eigenvector in AO basis (MPQC-ordered, in MPQC frame)").c_str());

    using psi::Chkpt;
    Chkpt::free(mopi);
    Chkpt::free(sopi);
    Chkpt::free(C);
    Chkpt::free(ao2so);
#endif

    return coefs_[spin_to_seek];
  }

  const Ref<OrbitalSpace>&  PsiSCF::orbs_sb(SpinCase1 spin) {
    if(orbs_sb_[spin]) {
      return orbs_sb_[spin];
    }
    if (this->reftype() == PsiSCF::rhf && spin==Beta) {
      return this->orbs_sb(Alpha);
    }

    const std::string id = ParsedOrbitalSpaceKey::key(std::string("pp(sym)"),spin);
    orbs_sb_[spin] = new OrbitalSpace(id,prepend_spincase(spin,"MOs (Psi3)"),
                                      this->coefs(spin),basis(),integral(),
                                      this->evals(spin),0,0,OrbitalSpace::symmetry);

    return orbs_sb_[spin];
  }

  RefSymmSCMatrix
  PsiSCF::density() {
    RefSymmSCMatrix result = this->alpha_density();
    if (!this->spin_polarized())
      result.scale(2.0);
    else
      result.accumulate(this->beta_density());
    return result;
  }

  RefSymmSCMatrix
  PsiSCF::alpha_density() {
    Ref<PetiteList> plist = integral()->petite_list();
    RefSCMatrix C_ao = coefs(Alpha);
    RefSCMatrix C_so = plist->evecs_to_SO_basis(C_ao);

    RefSymmSCMatrix P = C_so.kit()->symmmatrix(C_so.rowdim()); P.assign(0.0);
    P.accumulate_transform(C_so, this->mo_density(Alpha));
    return P;
  }

  RefSymmSCMatrix
  PsiSCF::beta_density() {
    Ref<PetiteList> plist = integral()->petite_list();
    RefSCMatrix C_ao = coefs(Beta);
    RefSCMatrix C_so = plist->evecs_to_SO_basis(C_ao);

    RefSymmSCMatrix P = C_so.kit()->symmmatrix(C_so.rowdim()); P.assign(0.0);
    P.accumulate_transform(C_so, this->mo_density(Beta));
    return P;
  }

  RefSymmSCMatrix
  PsiSCF::ao_density(SpinCase1 spin) {
    Ref<PetiteList> plist = integral()->petite_list();
    RefSCMatrix C_ao = coefs(spin);

    RefSymmSCMatrix P = C_ao.kit()->symmmatrix(C_ao.rowdim()); P.assign(0.0);
    P.accumulate_transform(C_ao, this->mo_density(spin));
    return P;
  }

  RefSymmSCMatrix
  PsiSCF::ao_density() {
    RefSymmSCMatrix result = this->ao_density(Alpha).copy();
    if (this->spin_polarized())
      result.accumulate(this->ao_density(Beta));
    else
      result.scale(2.0);
    return result;
  }

  RefSymmSCMatrix
  PsiSCF::mo_density(SpinCase1 spin) {
    if ((spin == Beta || spin == AnySpinCase1) && !this->spin_polarized())
      return mo_density(Alpha);

    if (spin == AnySpinCase1 && this->spin_polarized())
      ProgrammingError("asked for any spin density but the density is spin-polarized",
                       __FILE__, __LINE__);

    if (mo_density_[spin])
      return mo_density_[spin];

    RefSCMatrix C_ao = coefs(spin);
    RefSymmSCMatrix P_mo = C_ao.kit()->symmmatrix(C_ao.coldim());
    P_mo.assign(0.0);
    const int nmo = P_mo.n();
    for(int mo=0; mo<nmo; ++mo)
      P_mo.set_element(mo, mo, (spin == Alpha) ? this->alpha_occupation(mo)
                                               : this->beta_occupation(mo) );
    mo_density_[spin] = P_mo;

    if (debug() >= DefaultPrintThresholds::mostN)
      mo_density_[spin].print(prepend_spincase(spin,"Psi MO density").c_str());

    return mo_density_[spin];
  }

  const std::vector<unsigned int>& PsiSCF::occpi(SpinCase1 spin) {
    const double energy = this->energy();  // make sure it's computed
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
    const double energy = this->energy();  // make sure it's computed
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
    const double energy = this->energy();  // make sure it's computed
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

  void PsiSCF::compute_occupations(SpinCase1 spin) {
    const std::vector<unsigned int>& orbspi = this->mopi();
    const std::vector<unsigned int>& occpi = this->occpi(spin);
    const unsigned int nmo = std::accumulate(orbspi.begin(), orbspi.end(), 0);
    occupation_[spin].resize(nmo);
    std::fill(occupation_[spin].begin(), occupation_[spin].end(), 0.0);

    const unsigned int nirrep = orbspi.size();
    typedef std::vector<double>::iterator iter;
    iter v = occupation_[spin].begin();
    for(int h=0; h<nirrep; ++h) {
      const int nocc = occpi[h];
      for(int o=0; o<nocc; ++o, ++v)
        *v = 1.0;
      v += orbspi[h] - nocc;
    }
  }

  double
  PsiSCF::occupation(int mo) {
    return alpha_occupation(mo) + beta_occupation(mo);
  }

  double
  PsiSCF::alpha_occupation(int mo) {
    if (occupation_[Alpha].empty()) compute_occupations(Alpha);
    return occupation_[Alpha].at(mo);
  }

  double
  PsiSCF::beta_occupation(int mo) {
    if (occupation_[Beta].empty()) compute_occupations(Beta);
    return occupation_[Beta].at(mo);
  }

  void PsiSCF::import_occupations(const Ref<OneBodyWavefunction>& obwfn) {

    RefSCDimension osodim = oso_dimension();
    const bool spin_unrestricted = obwfn->spin_unrestricted();

    // make sure that the number of electrons is the same
    const int nuclear_charge = static_cast<int>(molecule()->total_Z());
    if (charge_ != (nuclear_charge - obwfn->nelectron()) )
      throw InputError("PsiSCF::import_occupations(obwfn) -- number of electrons in obwfn does not match this");

    // extract occupations
    std::vector<unsigned int> docc_obwfn(nirrep_, 0u);
    std::vector<unsigned int> socc_obwfn(nirrep_, 0u);
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
    multp_ = 1 + accumulate(socc_.begin(), socc_.end(), 0);

  }

  void
  PsiSCF::obsolete() {
    orbs_sb_[Alpha] = orbs_sb_[Beta] = 0;
    evals_[Alpha] = evals_[Beta] = 0;
    coefs_[Alpha] = coefs_[Beta] = 0;
    mo_density_[Alpha] = mo_density_[Beta] = 0;
    PsiWavefunction::obsolete();
  }

  void
  PsiSCF::symmetry_changed() {
    PsiWavefunction::symmetry_changed();
    occupation_[Alpha].resize(0); occupation_[Beta].resize(0);
    occpi_[Alpha].resize(0); occpi_[Beta].resize(0);
    uoccpi_[Alpha].resize(0); uoccpi_[Beta].resize(0);
    mopi_.resize(0);
    docc_.resize(0);
    socc_.resize(0);
    if (guess_wfn_) guess_wfn_->symmetry_changed();
  }

  //////////////////////////////////////////////////////////////////////////

  static ClassDesc PsiCLHF_cd(typeid(PsiCLHF), "PsiCLHF", 1, "public PsiSCF",
                              0, create<PsiCLHF>, create<PsiCLHF>);

  PsiCLHF::PsiCLHF(const Ref<KeyVal>&keyval) :
    PsiSCF(keyval) {

    multp_ = 1;
    const int nuclear_charge = molecule()->total_Z();
    if (docc_.empty()) {
      charge_ = keyval->intvalue("total_charge",KeyValValueint(0));
      if ( (nuclear_charge + charge_) % 2 != 0) {
        throw InputError("PsiCLHF::PsiCLHF -- odd number of electrons, total_charge keyword is not given or incorrect");
      }

      // try guess wavefunction
      if (guess_wfn_)
        import_occupations(guess_wfn_);
    }
    else {
      const int nelectron = accumulate(docc_.begin(), docc_.end(), 0) * 2;
      charge_ = nuclear_charge - nelectron;
    }
  }

  PsiCLHF::~PsiCLHF() {
  }

  PsiCLHF::PsiCLHF(StateIn&s) :
    PsiSCF(s) {
  }

  void PsiCLHF::print(std::ostream& os) const {
    os << indent << "PsiCLHF:\n" << incindent
      << indent << "total_charge = " << charge_ << endl;

    std::vector<unsigned int> occ = docc_;
    if (docc_.empty())
      occ = const_cast<PsiCLHF*>(this)->occpi(Alpha);

    os << indent << "docc = [";
    for (int i=0; i < nirrep_; i++)
      os << " " << occ[i];
    os << " ]" << endl;

    PsiWavefunction::print(os);
    os << decindent << endl;
  }

  void PsiCLHF::write_basic_input(int convergence) {
    Ref<PsiInput> input = get_psi_input();
    input->write_keyword("psi:reference", "rhf");
    input->write_keyword("psi:multp", multp_);
    input->write_keyword("psi:charge", charge_);
    if (!docc_.empty())
      input->write_keyword_array("psi:docc", docc_);
    else {
      input->write_keyword("psi:reset_occupations", true);
      input->write_keyword("psi:hcore_guess", "new");
    }
    input->write_keyword("scf:maxiter", maxiter_);
    if(diisdamp_ > 0) input->write_keyword("scf:diisdamp", diisdamp_);
    input->write_keyword("scf:convergence",convergence);
    if(levelshift_ > 1.0) {
      input->write_keyword("scf:levelshift", levelshift_);
      input->write_keyword("scf:levelshift_force", true);
    }
    if(!diis_) input->write_keyword("scf:diis", false);
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
    const int nuclear_charge = molecule()->total_Z();
    if (docc_.empty() || socc_.empty()) {
      charge_ = keyval->intvalue("total_charge",KeyValValueint(0));
      // lowest possible multp is 1 (when # of electrons is even) or 2 (odd)
      const int lowest_possible_multp = 1 + (nuclear_charge - charge_)%2;
      multp_ = keyval->intvalue("multiplicity",KeyValValueint(lowest_possible_multp));
      if ((nuclear_charge + charge_)%2 != (multp_ - 1)%2) {
        throw InputError("PsiHSOSHF::PsiHSOSHF -- inconsistent total_charge and multiplicty");
      }

      // try guess wavefunction
      if (guess_wfn_)
        import_occupations(guess_wfn_);
    }
    else {
      const int nsocc = accumulate(socc_.begin(), socc_.end(), 0);
      const int nelectron = accumulate(docc_.begin(), docc_.end(), 0) * 2 + nsocc;
      charge_ = nuclear_charge - nelectron;
      multp_ = nsocc + 1;
    }

  }

  PsiHSOSHF::~PsiHSOSHF() {
  }

  PsiHSOSHF::PsiHSOSHF(StateIn&s) :
    PsiSCF(s) {
  }

  void PsiHSOSHF::print(std::ostream& os) const {
    os << indent << "PsiHSOSHF:\n" << incindent;
    os << indent << "total_charge = " << charge_ << endl;
    os << indent << "multiplicity = " << multp_ << endl;

    std::vector<unsigned int> docc = docc_;
    std::vector<unsigned int> socc = socc_;
    if (docc_.empty())
      docc = const_cast<PsiHSOSHF*>(this)->occpi(Beta);
    if (socc_.empty()) {
      std::vector<unsigned int> docc_a = const_cast<PsiHSOSHF*>(this)->occpi(Alpha);
      socc.resize(docc_a.size());
      std::transform(docc_a.begin(), docc_a.end(), docc.begin(), socc.begin(),
                     std::minus<unsigned int>());
    }

    os << indent << "docc = [";
    for (int i=0; i < nirrep_; i++)
      os << " " << docc[i];
    os << " ]" << endl;
    os << indent << "socc = [";
    for (int i=0; i < nirrep_; i++)
      os << " " << socc[i];
    os << " ]" << endl;

    PsiWavefunction::print(os);
    os << decindent << endl;
  }

  void PsiHSOSHF::write_basic_input(int convergence) {
    Ref<PsiInput> input = get_psi_input();
    input->write_keyword("psi:reference", "rohf");
    if (!docc_.empty())
      input->write_keyword_array("psi:docc", docc_);
    if (!socc_.empty())
      input->write_keyword_array("psi:socc", socc_);
    input->write_keyword("psi:multp", multp_);
    input->write_keyword("psi:charge", charge_);
    if (docc_.empty() && socc_.empty()) {
      input->write_keyword("psi:reset_occupations", true);
      input->write_keyword("psi:hcore_guess", "new");
    }
    input->write_keyword("scf:maxiter", maxiter_);
    input->write_keyword("scf:convergence",convergence);
    if(levelshift_ > 1.0) input->write_keyword("scf:levelshift", levelshift_);
    if(!diis_) input->write_keyword("scf:diis", false);
  }

  void PsiHSOSHF::write_input(int convergence) {
    Ref<PsiInput> input = get_psi_input();
    input->open();
    PsiWavefunction::write_basic_input(convergence);
    write_basic_input(convergence);
    input->write_keyword("psi:wfn", "scf");
    input->close();
  }

  const RefSCMatrix&
  PsiHSOSHF::coefs_semicanonical(SpinCase1 s) {
    if (coefs_sc_[s])
      return coefs_sc_[s];
    semicanonical();
    return coefs_sc_[s];
  }

  const RefDiagSCMatrix&
  PsiHSOSHF::evals_semicanonical(SpinCase1 s) {
    if (evals_sc_[s])
      return evals_sc_[s];
    semicanonical();
    return evals_sc_[s];
  }

#define DEBUG_SEMICANONICAL 0
  void
  PsiHSOSHF::semicanonical()
  {
    // May need to update this
    if (value_needed())
      compute();

    // if orbitals are available read off disk
    PsiChkpt chkpt(exenv(), integral(), debug());
    if (chkpt.have_spin_unrestricted_mos()) {
      const bool seek_spin_restricted = false;
      evals_sc_[Alpha] = chkpt.evals(Alpha, seek_spin_restricted);
      evals_sc_[Beta] = chkpt.evals(Beta, seek_spin_restricted);
      coefs_sc_[Alpha] = chkpt.coefs(Alpha, seek_spin_restricted);
      coefs_sc_[Beta] = chkpt.coefs(Beta, seek_spin_restricted);
      //coefs_sc_[Alpha].print("Alpha semicanonical MOs (computed in Psi)");
      //coefs_sc_[Beta].print("Beta semicanonical MOs (computed in Psi)");
      return;
    }

    const RefSCDimension modim = evals().dim();
    const RefSCDimension sodim = so_dimension();
    const unsigned int nirrep = molecule()->point_group()->char_table().nirrep();

    Ref<PetiteList> plist = integral()->petite_list();
    RefSCMatrix ao_eigenvector = coefs();
    RefSCMatrix so_eigenvector = plist->evecs_to_SO_basis(ao_eigenvector);
  #if DEBUG_SEMICANONICAL
    evals().print("PsiHSOSSCF orbital energies");
    ao_eigenvector.print("PsiHSOSSCF orbitals (AO basis)");
  #endif

    // initialize dimensions for occupied and virtual blocks
    int* aoccpi = new int[nirrep];  for(unsigned int h=0; h<nirrep; ++h) aoccpi[h] = 0;
    int* boccpi = new int[nirrep];  for(unsigned int h=0; h<nirrep; ++h) boccpi[h] = 0;
    int* avirpi = new int[nirrep];  for(unsigned int h=0; h<nirrep; ++h) avirpi[h] = 0;
    int* bvirpi = new int[nirrep];  for(unsigned int h=0; h<nirrep; ++h) bvirpi[h] = 0;
    for(unsigned int h=0, ii=0; h<nirrep; ++h) {
      const unsigned int n = modim->blocks()->size(h);
      for(unsigned int i=0; i<n; ++i, ++ii) {
        const double o = occupation(ii);
        if (o == 2.0) {
          aoccpi[h] += 1;
          boccpi[h] += 1;
        }
        if (o == 1.0) {
          aoccpi[h] += 1;
          bvirpi[h] += 1;
        }
        if (o == 0.0) {
          avirpi[h] += 1;
          bvirpi[h] += 1;
        }
      }
    }
    int naocc = 0;  for(unsigned int h=0; h<nirrep; ++h) naocc += aoccpi[h];
    int nbocc = 0;  for(unsigned int h=0; h<nirrep; ++h) nbocc += boccpi[h];
    int navir = 0;  for(unsigned int h=0; h<nirrep; ++h) navir += avirpi[h];
    int nbvir = 0;  for(unsigned int h=0; h<nirrep; ++h) nbvir += bvirpi[h];
    const RefSCDimension aoccdim = new SCDimension(naocc,nirrep,aoccpi,"I");
    const RefSCDimension boccdim = new SCDimension(nbocc,nirrep,boccpi,"i");
    const RefSCDimension avirdim = new SCDimension(navir,nirrep,avirpi,"A");
    const RefSCDimension bvirdim = new SCDimension(nbvir,nirrep,bvirpi,"a");
    for(int h=0; h<nirrep; ++h) {
      aoccdim->blocks()->set_subdim(h,new SCDimension(aoccpi[h]));
      boccdim->blocks()->set_subdim(h,new SCDimension(boccpi[h]));
      avirdim->blocks()->set_subdim(h,new SCDimension(avirpi[h]));
      bvirdim->blocks()->set_subdim(h,new SCDimension(bvirpi[h]));
    }

    // get occupied and virtual eigenvectors
    RefSCMatrix aoccvec(so_dimension(), aoccdim, basis_matrixkit()); aoccvec.assign(0.0);
    RefSCMatrix boccvec(so_dimension(), boccdim, basis_matrixkit()); boccvec.assign(0.0);
    RefSCMatrix avirvec(so_dimension(), avirdim, basis_matrixkit()); avirvec.assign(0.0);
    RefSCMatrix bvirvec(so_dimension(), bvirdim, basis_matrixkit()); bvirvec.assign(0.0);
    {
      for(unsigned int h=0; h<nirrep; ++h) {
        const unsigned int nmo = modim->blocks()->size(h);
        const unsigned int nso = sodim->blocks()->size(h);

        if(aoccpi[h]) aoccvec.block(h).assign_subblock(so_eigenvector.block(h), 0, nso-1, 0, aoccpi[h]-1, 0, 0);
        if(boccpi[h]) boccvec.block(h).assign_subblock(so_eigenvector.block(h), 0, nso-1, 0, boccpi[h]-1, 0, 0);
        if(avirpi[h]) avirvec.block(h).assign_subblock(so_eigenvector.block(h), 0, nso-1, 0, avirpi[h]-1, 0, aoccpi[h]);
        if(bvirpi[h]) bvirvec.block(h).assign_subblock(so_eigenvector.block(h), 0, nso-1, 0, bvirpi[h]-1, 0, boccpi[h]);

      }
    }

    // get the Fock matrices in AO basis, then transform to SO basis
    RefSymmSCMatrix Pa = this->ao_density(Alpha);
    RefSymmSCMatrix Pb = this->ao_density(Beta);
    RefSymmSCMatrix P = Pa + Pb;
    RefSymmSCMatrix Po = Pa - Pb;
    Ref<TwoBodyFockMatrixBuilder<true> >fbuild = new TwoBodyFockMatrixBuilder<true>(true, true, true,
                                                                                    basis(), basis(), basis(),
                                                                                    P, Po,
                                                                                    integral(),
                                                                                    MessageGrp::get_default_messagegrp(),
                                                                                    ThreadGrp::get_default_threadgrp(),
                                                                                    this->desired_value_accuracy() / 1000.0);
    RefSymmSCMatrix Fa_ao = fbuild->F(Alpha);
    RefSymmSCMatrix Fb_ao = fbuild->F(Beta);
    fbuild = 0;
    Ref< OneBodyFockMatrixBuilder<true> > hbuild = new OneBodyFockMatrixBuilder<true>(OneBodyFockMatrixBuilder<true>::NonRelativistic,
                                                                                      basis(),
                                                                                      basis(),
                                                                                      0,
                                                                                      integral(),
                                                                                      this->desired_value_accuracy() / 1000.0);
    Fa_ao.accumulate(hbuild->result());
    Fb_ao.accumulate(hbuild->result());

    //
    // alpha case
    //
    // diagonalize occ-occ and virt-virt fock matrices
    RefSymmSCMatrix afock_so = plist->to_SO_basis(Fa_ao);
  #if DEBUG_SEMICANONICAL
    {
      RefSymmSCMatrix Fa(oso_dimension(), basis_matrixkit());
      Fa.assign(0.0);
      Fa.accumulate_transform(so_eigenvector, afock_so,
                              SCMatrix::TransposeTransform);
      Fa.print("Alpha Fock matrix (MO basis)");
    }
  #endif

    RefSymmSCMatrix aoccfock(aoccdim, basis_matrixkit());
    aoccfock.assign(0.0);
    aoccfock.accumulate_transform(aoccvec, afock_so,
                                  SCMatrix::TransposeTransform);
    RefDiagSCMatrix aoccevals(aoccdim, basis_matrixkit());
    RefSCMatrix aoccevecs(aoccdim, aoccdim, basis_matrixkit());
    aoccfock.diagonalize(aoccevals, aoccevecs);

    RefSymmSCMatrix avirfock(avirdim, basis_matrixkit());
    avirfock.assign(0.0);
    avirfock.accumulate_transform(avirvec, afock_so,
                                  SCMatrix::TransposeTransform);
    RefDiagSCMatrix avirevals(avirdim, basis_matrixkit());
    RefSCMatrix avirevecs(avirdim, avirdim, basis_matrixkit());
    avirfock.diagonalize(avirevals, avirevecs);
    // form full eigenvectors and eigenvalues
    if (evals_sc_[Alpha] == 0) {
      evals_sc_[Alpha] = evals(Alpha).clone();  evals_sc_[Alpha].assign(0.0);
    }
    RefSCMatrix aevecs(modim, modim, basis_matrixkit()); aevecs.assign(0.0);
    {
      unsigned int aoccoffset = 0;
      unsigned int aviroffset = 0;
      for(unsigned int h=0; h<nirrep; ++h) {
        const unsigned int mostart = modim->blocks()->start(h);
        const unsigned int moend = modim->blocks()->fence(h);
        const unsigned int nmo = modim->blocks()->size(h);

        if (aoccpi[h]) aevecs.block(h).assign_subblock(aoccevecs.block(h), 0, aoccpi[h]-1, 0, aoccpi[h]-1, 0, 0);
        if (avirpi[h]) aevecs.block(h).assign_subblock(avirevecs.block(h), aoccpi[h], nmo-1, aoccpi[h], nmo-1, 0, 0);

        int i = 0;
        for(unsigned int mo=mostart; mo<mostart+aoccpi[h]; ++mo, ++i) {
          const double e = aoccevals.get_element(aoccoffset + i);
          evals_sc_[Alpha].set_element(mo,e);
        }
        i = 0;
        for(unsigned int mo=mostart+aoccpi[h]; mo<moend; ++mo, ++i) {
          const double e = avirevals.get_element(aviroffset + i);
          evals_sc_[Alpha].set_element(mo,e);
        }

        aoccoffset += aoccpi[h];
        aviroffset += avirpi[h];
      }
    }
    coefs_sc_[Alpha] = ao_eigenvector * aevecs;
    canonicalize_column_phases(coefs_sc_[Alpha]);
  #if DEBUG_SEMICANONICAL
    {
      ao_eigenvector.print("Original eigenvector (AO basis)");
      aevecs.print("Alpha transform matrix");
      {
        RefSymmSCMatrix afock_mo(modim, basis_matrixkit());
        afock_mo.assign(0.0);
        afock_mo.accumulate_transform(coefs_sc_[Alpha], Fa_ao,
                                      SCMatrix::TransposeTransform);
        afock_mo.print("Alpha Fock matrix in the semicanonical orbitals");
      }
      evals_sc_[Alpha].print("Alpha semicanonical orbital energies");
      coefs_sc_[Alpha].print("Alpha semicanonical orbitals (AO basis)");
    }
  #endif
    aevecs = 0;

    //
    // beta case
    //
    // diagonalize occ-occ and virt-virt fock matrices
    RefSymmSCMatrix bfock_so = plist->to_SO_basis(Fb_ao);
  #if DEBUG_SEMICANONICAL
    {
      RefSymmSCMatrix Fb(oso_dimension(), basis_matrixkit());
      Fb.assign(0.0);
      Fb.accumulate_transform(so_eigenvector, bfock_so,
                              SCMatrix::TransposeTransform);
      Fb.print("Beta Fock matrix");
    }
  #endif

    RefSymmSCMatrix boccfock(boccdim, basis_matrixkit());
    boccfock.assign(0.0);
    boccfock.accumulate_transform(boccvec, bfock_so,
                                  SCMatrix::TransposeTransform);
    RefDiagSCMatrix boccevals(boccdim, basis_matrixkit());
    RefSCMatrix boccevecs(boccdim, boccdim, basis_matrixkit());
    boccfock.diagonalize(boccevals, boccevecs);

    RefSymmSCMatrix bvirfock(bvirdim, basis_matrixkit());
    bvirfock.assign(0.0);
    bvirfock.accumulate_transform(bvirvec, bfock_so,
                                  SCMatrix::TransposeTransform);
    RefDiagSCMatrix bvirevals(bvirdim, basis_matrixkit());
    RefSCMatrix bvirevecs(bvirdim, bvirdim, basis_matrixkit());
    bvirfock.diagonalize(bvirevals, bvirevecs);
    // form full eigenvectors and eigenvalues
    if (evals_sc_[Beta] == 0) {
      evals_sc_[Beta] = evals(Beta).clone();  evals_sc_[Beta].assign(0.0);
    }
    RefSCMatrix bevecs(modim, modim, basis_matrixkit()); bevecs.assign(0.0);
    {
      unsigned int boccoffset = 0;
      unsigned int bviroffset = 0;
      for(unsigned int h=0; h<nirrep; ++h) {
        const unsigned int mostart = modim->blocks()->start(h);
        const unsigned int moend = modim->blocks()->fence(h);
        const unsigned int nmo = modim->blocks()->size(h);

        if (boccpi[h]) bevecs.block(h).assign_subblock(boccevecs.block(h), 0, boccpi[h]-1, 0, boccpi[h]-1, 0, 0);
        if (bvirpi[h]) bevecs.block(h).assign_subblock(bvirevecs.block(h), boccpi[h], nmo-1, boccpi[h], nmo-1, 0, 0);

        int i = 0;
        for(unsigned int mo=mostart; mo<mostart+boccpi[h]; ++mo, ++i) {
          const double e = boccevals.get_element(boccoffset + i);
          evals_sc_[Beta].set_element(mo,e);
        }
        i = 0;
        for(unsigned int mo=mostart+boccpi[h]; mo<moend; ++mo, ++i) {
          const double e = bvirevals.get_element(bviroffset + i);
          evals_sc_[Beta].set_element(mo,e);
        }

        boccoffset += boccpi[h];
        bviroffset += bvirpi[h];
      }
    }
    coefs_sc_[Beta] = ao_eigenvector * bevecs;
    canonicalize_column_phases(coefs_sc_[Beta]);
  #if DEBUG_SEMICANONICAL
    {
      bevecs.print("Beta transform matrix");
      {
        RefSymmSCMatrix bfock_mo(modim, basis_matrixkit());
        bfock_mo.assign(0.0);
        bfock_mo.accumulate_transform(coefs_sc_[Beta], Fb_ao,
                                      SCMatrix::TransposeTransform);
        bfock_mo.print("Beta Fock matrix in the semicanonical orbitals");
      }
      evals_sc_[Beta].print("Beta semicanonical orbital energies");
      coefs_sc_[Beta].print("Beta semicanonical orbitals (AO basis)");
    }
  #endif
    bevecs = 0;

    //coefs_sc_[Alpha].print("Alpha semicanonical MOs (computed in MPQC)");
    //coefs_sc_[Beta].print("Beta semicanonical MOs (computed in MPQC)");
  }

  void
  PsiHSOSHF::obsolete() {
    coefs_sc_[Alpha] = coefs_sc_[Beta] = 0;
    evals_sc_[Alpha] = evals_sc_[Beta] = 0;
    PsiSCF::obsolete();
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
    const int nuclear_charge = molecule()->total_Z();
    if (docc_.empty() || socc_.empty()) {
      charge_ = keyval->intvalue("total_charge",KeyValValueint(0));
      // lowest possible multp is 1 (when # of electrons is even) or 2 (odd)
      const int lowest_possible_multp = 1 + (nuclear_charge - charge_)%2;
      multp_ = keyval->intvalue("multiplicity",KeyValValueint(lowest_possible_multp));
      if ((nuclear_charge + charge_)%2 != (multp_ - 1)%2) {
        throw InputError("PsiUHF::PsiUHF -- inconsistent total_charge and multiplicity");
      }

      // try guess wavefunction
      if (guess_wfn_)
        import_occupations(guess_wfn_);
    }
    else {
      const int nsocc = accumulate(socc_.begin(), socc_.end(), 0);
      const int nelectron = accumulate(docc_.begin(), docc_.end(), 0) * 2 + nsocc;
      charge_ = nuclear_charge - nelectron;
      multp_ = nsocc + 1;
    }

  }

  PsiUHF::~PsiUHF() {
  }

  PsiUHF::PsiUHF(StateIn&s) :
    PsiSCF(s) {
  }

  void PsiUHF::print(std::ostream& os) const {
    os << indent << "PsiUHF:\n" << incindent;
    os << indent << "total_charge = " << charge_ << endl;
    os << indent << "multiplicity = " << multp_ << endl;

    std::vector<unsigned int> docc = docc_;
    std::vector<unsigned int> socc = socc_;
    if (docc.empty())
      docc = const_cast<PsiUHF*>(this)->occpi(Beta);
    if (socc.empty()) {
      std::vector<unsigned int> docc_a = const_cast<PsiUHF*>(this)->occpi(Alpha);
      socc.resize(docc_a.size());
      std::transform(docc_a.begin(), docc_a.end(), docc.begin(), socc.begin(),
                     std::minus<unsigned int>());
    }

    os << indent << "docc = [";
    for (int i=0; i < nirrep_; i++)
      os << " " << docc[i];
    os << " ]" << endl;
    os << indent << "socc = [";
    for (int i=0; i < nirrep_; i++)
      os << " " << socc[i];
    os << " ]" << endl;

    PsiWavefunction::print(os);
    os << decindent << endl;
  }

  void PsiUHF::write_basic_input(int convergence) {
    Ref<PsiInput> input = get_psi_input();
    input->write_keyword("psi:reference", "uhf");
    if (!docc_.empty())
      input->write_keyword_array("psi:docc", docc_);
    if (!socc_.empty())
      input->write_keyword_array("psi:socc", socc_);
    input->write_keyword("psi:multp", multp_);
    input->write_keyword("psi:charge", charge_);
    if (docc_.empty() && socc_.empty()) {
      input->write_keyword("psi:reset_occupations", true);
      input->write_keyword("psi:hcore_guess", "new");
    }
    input->write_keyword("scf:maxiter", maxiter_);
    input->write_keyword("scf:convergence",convergence);
    if(levelshift_ > 1.0) input->write_keyword("scf:levelshift", levelshift_);
    if(!diis_) input->write_keyword("scf:diis", false);
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

    if(keyval->exists("frozen_docc")) {
      int len = keyval->count("frozen_docc");
      if(len != nirrep()){
        ExEnv::err0() << "PsiCorrWavefunction::PsiCorrWavefunction: length of frozen_docc does not coincide with number of irreps." << endl;
        abort();
      }
      frozen_docc_ = std::vector<unsigned int>(len);
      for(int i=0; i<len; i++){
        frozen_docc_[i] = keyval->intvalue("frozen_docc",i);
      }
      nfzc_ = std::accumulate(frozen_docc_.begin(), frozen_docc_.end(), 0u);
    }
    else { // frozen_docc is not given -- need nfzc, will compute frozen_docc later
      std::string nfzc_str = keyval->stringvalue("nfzc",KeyValValuestring("0"));
      if (nfzc_str == "auto")
        nfzc_ = molecule()->n_core_electrons()/2;
      else if (nfzc_str == "no" || nfzc_str == "false")
        nfzc_ = 0;
      else
        nfzc_ = atoi(nfzc_str.c_str());
    }

    if(keyval->exists("frozen_uocc")) {
      int len = keyval->count("frozen_uocc");
      if(len != nirrep()){
        ExEnv::err0() << "PsiCorrWavefunction::PsiCorrWavefunction: length of frozen_docc does not coincide with number of irreps." << endl;
        abort();
      }
      frozen_uocc_ = std::vector<unsigned int>(len);
      for(int i=0; i<len; i++){
        frozen_uocc_[i] = keyval->intvalue("frozen_uocc",i);
      }
      nfzv_ = std::accumulate(frozen_uocc_.begin(), frozen_uocc_.end(), 0u);
    }
    else { // frozen_uocc is not given -- need nfzv, will compute frozen_uocc later
      std::string nfzv_str = keyval->stringvalue("nfzv",KeyValValuestring("0"));
      if (nfzv_str == "no" || nfzv_str == "false")
        nfzv_ = 0;
      else
        nfzv_ = atoi(nfzv_str.c_str());
    }
    if (nfzv_ != 0)
      throw FeatureNotImplemented("Implementation of frozen virtuals in Psi is not reliable, set nfzv = 0",
                                  __FILE__,__LINE__);

    reference_ << keyval->describedclassvalue("reference");
    if (reference_ == 0) {
      ExEnv::err0()
          << "PsiCorrWavefunction::PsiCorrWavefunction: no reference wavefunction"
          << endl;
      abort();
    }

    set_desired_value_accuracy(desired_value_accuracy());
  }

  PsiCorrWavefunction::~PsiCorrWavefunction() {
  }

  PsiCorrWavefunction::PsiCorrWavefunction(StateIn&s) :
    PsiWavefunction(s) {

    reference_ << SavableState::restore_state(s);
    s.get(frozen_docc_);
    s.get(frozen_uocc_);
    int nfzc; s.get(nfzc); nfzc_ = static_cast<unsigned int>(nfzc);
    int nfzv; s.get(nfzv); nfzv_ = static_cast<unsigned int>(nfzv);
  }

  void PsiCorrWavefunction::save_data_state(StateOut&s) {
    PsiWavefunction::save_data_state(s);
    SavableState::save_state(reference_.pointer(), s);
    s.put(frozen_docc_);
    s.put(frozen_uocc_);
    s.put(static_cast<int>(nfzc_));
    s.put(static_cast<int>(nfzv_));
  }

  void PsiCorrWavefunction::print(std::ostream& os) const {
    os << indent << "PsiCorrWavefunction:" << endl;
    os << incindent;
    PsiWavefunction::print(os);
    reference_->print(os);

    const std::vector<unsigned int>& fzdocc = this->frozen_docc();
    const std::vector<unsigned int>& fzuocc = this->frozen_uocc();
    const int nirrep = fzdocc.size();
    os << indent << "frozen_docc = [" << fzdocc[0];
    for (int i=1; i < nirrep_; i++)
      os << " " << fzdocc[i];
    os << " ]" << endl;
    os << indent << "frozen_uocc = [" << fzuocc[0];
    for (int i=1; i < nirrep_; i++)
      os << " " << fzuocc[i];
    os << " ]" << endl;

    os << decindent;
  }

  void PsiCorrWavefunction::set_desired_value_accuracy(double acc) {
    Function::set_desired_value_accuracy(acc);
    if (reference()->desired_value_accuracy_set_to_default())
      reference()->set_desired_value_accuracy( acc /
                                               valacc_to_refacc() );
  }

  RefSymmSCMatrix PsiCorrWavefunction::density() {
    // read correlated density from Psi3
    throw FeatureNotImplemented("densities not implemented for Psi-based correlated Wavefunctions",
                                __FILE__, __LINE__, this->class_desc());
    return 0;
  }

  void PsiCorrWavefunction::write_input(int convergence) {
    this->write_input_frozen2restricted(convergence, false);
  }

  namespace {
    int accuracy_to_convergence( double acc ) {
      return -(int)(ceil(log10(acc)));
    }
  }

  void PsiCorrWavefunction::write_input_frozen2restricted(int convergence,
                                                          bool frozen2restricted) {
    if (gradient_needed())
      reference_->do_gradient(1);
    else
      reference_->do_gradient(0);

    Ref<PsiInput> input = get_psi_input();
    PsiWavefunction::write_basic_input(convergence);
    reference_->write_basic_input( accuracy_to_convergence(reference_->desired_value_accuracy()) );
    input->write_keyword("psi:tolerance", convergence + 5);
    // For CCSD density, frozen core can not be used
    if (frozen2restricted) {
        if (compute_1rdm_) {
          // CCSD density can not deal with frozen core
          const std::vector<unsigned int> frozen_docc(nirrep_,0);
          input->write_keyword_array("psi:restricted_docc", frozen_docc);
        } else {
            input->write_keyword_array("psi:restricted_docc", this->frozen_docc());
            //input->write_keyword_array("psi:restricted_uocc", this->frozen_uocc());
        }
    }
    else {
      if (compute_1rdm_) {
        // CCSD density can not deal with frozen core
        const std::vector<unsigned int> frozen_docc(nirrep_,0);
        input->write_keyword_array("psi:frozen_docc", frozen_docc);
        input->write_keyword_array("psi:frozen_uocc", this->frozen_uocc());
      } else {
          input->write_keyword_array("psi:frozen_docc", this->frozen_docc());
          input->write_keyword_array("psi:frozen_uocc", this->frozen_uocc());
      }
    }

  }

  void
  PsiCorrWavefunction::compute() {
    // compute reference NOW so that orbital info is available for making the correlated wfn input
    reference_->compute();
    PsiWavefunction::compute();
    const double escf = exenv()->chkpt().rd_escf();
    reference_->set_energy(escf);
    reference_->set_actual_value_accuracy(
        reference_->desired_value_accuracy()
        );
  }

  void PsiCorrWavefunction::obsolete() {
    reference_->obsolete();
    PsiWavefunction::obsolete();
  }
  void PsiCorrWavefunction::symmetry_changed() {
    PsiWavefunction::symmetry_changed();
    reference_->symmetry_changed();
    frozen_docc_.resize(0);
    frozen_uocc_.resize(0);
  }


  const Ref<PsiSCF>&
  PsiCorrWavefunction::reference() const {
    return reference_;
  }

  double PsiCorrWavefunction::magnetic_moment() const {
    return reference_->magnetic_moment();
  }

  int PsiCorrWavefunction::nelectron() {
    return reference()->nelectron();
  }

  const Ref<OrbitalSpace>&  PsiCorrWavefunction::orbs_sb(SpinCase1 spin) {
    return reference_->orbs_sb(spin);
  }

  unsigned int
  PsiCorrWavefunction::nfzc() const { return nfzc_; }
  unsigned int
  PsiCorrWavefunction::nfzv() const { return nfzv_; }

  const std::vector<unsigned int>&
  PsiCorrWavefunction::frozen_docc() const {
    if (frozen_docc_.empty()) {
      frozen_docc_.resize(nirrep_);

      const std::vector<unsigned int> mopi = reference()->mopi();
      RefDiagSCMatrix evals_a = reference()->evals(Alpha);
      typedef MolecularOrbitalMask<double, RefDiagSCMatrix> FZCMask;
      FZCMask fzcmask(nfzc(), evals_a);
      // count number of frozen_core orbitals in each irrep
      typedef std::vector<bool>::const_iterator iter;
      iter v = fzcmask.mask().begin();
      for(int h=0; h<nirrep_; ++h) {
        const unsigned int nmo = mopi[h];
        const int n = std::count(v, v+nmo, false);
        frozen_docc_[h] = n;
        v += nmo;
      }
    }
    return frozen_docc_;
  }

  const std::vector<unsigned int>&
  PsiCorrWavefunction::frozen_uocc() const {
    if (frozen_uocc_.empty()) {
      frozen_uocc_.resize(nirrep_);

      const std::vector<unsigned int> mopi = reference()->mopi();
      RefDiagSCMatrix evals_a = reference()->evals(Alpha);
      typedef MolecularOrbitalMask<double, RefDiagSCMatrix, std::greater<double> > FZVMask;
      FZVMask fzvmask(nfzv(), evals_a);
      // count number of frozen virtual orbitals in each irrep
      typedef std::vector<bool>::const_iterator iter;
      iter v = fzvmask.mask().begin();
      for(int h=0; h<nirrep_; ++h) {
        const unsigned int nmo = mopi[h];
        const int n = std::count(v, v+nmo, false);
        frozen_uocc_[h] = n;
        v += nmo;
      }
    }
    return frozen_uocc_;
  }

  const std::vector<unsigned int> PsiCorrWavefunction::docc_act() {
    std::vector<unsigned int> occ_act_alpha = reference_->occpi(Alpha);
    unsigned int occ_act_alpha_sum = std::accumulate(occ_act_alpha.begin(),occ_act_alpha.end(),0);
    std::vector<unsigned int> occ_act_beta = reference_->occpi(Beta);
    unsigned int occ_act_beta_sum = std::accumulate(occ_act_beta.begin(),occ_act_beta.end(),0);
    unsigned int nirrep = occ_act_alpha.size();
    std::vector<unsigned int> docc_active(nirrep);
    if(occ_act_alpha_sum > occ_act_beta_sum) {
      // docc_active = occ_act_beta - frozen_docc()
      std::transform(occ_act_beta.begin(), occ_act_beta.end(), frozen_docc().begin(),
                     docc_active.begin(),
                     std::minus<unsigned int>());
    }
    else {
      // docc_active = occ_act_alpha - frozen_docc()
      std::transform(occ_act_alpha.begin(), occ_act_alpha.end(), frozen_docc().begin(),
                     docc_active.begin(),
                     std::minus<unsigned int>());
    }

    return(docc_active);
  }

  const std::vector<unsigned int> PsiCorrWavefunction::socc() {
    std::vector<unsigned int> occ_act_alpha = reference_->occpi(Alpha);
    unsigned int occ_act_alpha_sum = std::accumulate(occ_act_alpha.begin(),occ_act_alpha.end(),0);
    std::vector<unsigned int> occ_act_beta = reference_->occpi(Beta);
    unsigned int occ_act_beta_sum = std::accumulate(occ_act_beta.begin(),occ_act_beta.end(),0);
    unsigned int nirrep = occ_act_alpha.size();
    std::vector<unsigned int> singocc(nirrep);
    if(occ_act_alpha_sum > occ_act_beta_sum) {
      // singocc = occ_act_alpha - occ_act_beta
      std::transform(occ_act_alpha.begin(), occ_act_alpha.end(), occ_act_beta.begin(),
                     singocc.begin(),
                     std::minus<unsigned int>());
    }
    else {
      // singocc = occ_act_beta - occ_act_alpha
      std::transform(occ_act_beta.begin(), occ_act_beta.end(), occ_act_alpha.begin(),
                     singocc.begin(),
                     std::minus<unsigned int>());
    }

    return(singocc);
  }

  const std::vector<unsigned int> PsiCorrWavefunction::uocc_act() {
    std::vector<unsigned int> uocc_act_alpha = reference_->uoccpi(Alpha);
    unsigned int uocc_act_alpha_sum = std::accumulate(uocc_act_alpha.begin(),uocc_act_alpha.end(),0);
    std::vector<unsigned int> uocc_act_beta = reference_->uoccpi(Beta);
    unsigned int uocc_act_beta_sum = std::accumulate(uocc_act_beta.begin(),uocc_act_beta.end(),0);
    int nirrep = uocc_act_alpha.size();
    std::vector<unsigned int> uocc_active(nirrep);
    if(uocc_act_alpha_sum > uocc_act_beta_sum){
      // uocc_active = uocc_act_beta - frozen_uocc()
      std::transform(uocc_act_beta.begin(), uocc_act_beta.end(), frozen_uocc().begin(),
                     uocc_active.begin(),
                     std::minus<unsigned int>());
    }
    else {
      // uocc_active = uocc_act_alpha - frozen_uocc()
      std::transform(uocc_act_alpha.begin(), uocc_act_alpha.end(), frozen_uocc().begin(),
                     uocc_active.begin(),
                     std::minus<unsigned int>());
    }

    return(uocc_active);
  }

  std::vector<unsigned int>
  PsiCorrWavefunction::map_density_to_sb() {
    // maps symm -> QT
    std::vector<unsigned int> fmap = index_map_symmtocorrorder(frozen_docc(),
                                                             docc_act(),
                                                             socc(),
                                                             uocc_act(),
                                                             frozen_uocc());
    return index_map_inverse( fmap );
  }

  double
  PsiCorrWavefunction::reference_energy()
  {
    return exenv()->chkpt().rd_eref();
  }

  namespace {
    /// A modification of Psi3's iwl_rdone. For documentation see that of iwl_rdone.
    int _rdopdm(int itap, char *label, double *ints, int ntri, int erase,
               int printflg, FILE *outfile) {
      int nmo;
      double TOL = 1.0e-6;

      psio_open(itap, PSIO_OPEN_OLD);
      psio_read_entry(itap, label, (char *) ints, ntri * sizeof(double));
      psio_close(itap, !erase);

      if (printflg) {
        nmo = round(sqrt((double) ntri));
        for (int i = 0; i < nmo; i++) {
          for (int j = 0; j < nmo; j++) {
            if (fabs(ints[i * nmo + j]) >= TOL) {
              fprintf(outfile, "%5i%5i%12.7f\n", i + 1, j + 1,
                      ints[i * nmo + j]);
            }
          }
        }
      }
      return (1);
    }
  } // end of anonymous namespace

  RefSymmSCMatrix sc::detail::rdopdm(SpinCase1 spin,
                                     const std::vector<unsigned int>& mopi,
                                     const std::vector<unsigned int>& dmap,
                                     Ref<SCMatrixKit> kit) {
    FILE* outfile;
    string opdm_label_str;
    if(spin==Alpha) {
      opdm_label_str = "MO-basis Alpha OPDM";
      outfile = fopen("psiout_opdm_Alpha","w");
    }
    else { // spin==Beta
      opdm_label_str = "MO-basis Beta OPDM";
      outfile = fopen("psiout_opdm_Beta","w");
    }
    char *opdm_label=(char *)opdm_label_str.c_str();
    const int nmo = dmap.size();
    using psi::Chkpt;
    double** c_opdm = Chkpt::matrix<double>(nmo, nmo);
    _rdopdm(PSIF_MO_OPDM,opdm_label,&(c_opdm[0][0]),nmo*nmo,0,1,outfile);
    fclose(outfile);

    const int nirrep = mopi.size();
    const unsigned int nmo_target = std::accumulate(mopi.begin(), mopi.end(), 0);
    RefSCDimension modim = new SCDimension(nmo_target,
                                           nirrep,
                                           reinterpret_cast<int*>(const_cast<unsigned int*>(&(mopi[0])))
                                          );
    for (unsigned int h=0; h<nirrep; ++h)
      modim->blocks()->set_subdim(h, new SCDimension(mopi[h]));

    if (kit == 0)
      kit = new BlockedSCMatrixKit(SCMatrixKit::default_matrixkit());
    RefSymmSCMatrix opdm = kit->symmmatrix(modim);

    // map density to symmetry-blocked indices
    for(unsigned int r=0; r<nmo; ++r) {
      const unsigned int rr = dmap[r];
      for(unsigned int c=0; c<nmo; ++c) {
        const unsigned int cc = dmap[c];
        opdm.set_element( rr, cc, c_opdm[r][c] );
      }
    }

    Chkpt::free(c_opdm);

    return opdm;
  }

  RefSymmSCMatrix sc::detail::rdtpdm(SpinCase2 pairspin,
                                     const std::vector<unsigned int>& dmap,
                                     bool spinfree,
                                     Ref<SCMatrixKit> kit) {

    if (spinfree)
      MPQC_ASSERT(pairspin == AlphaBeta);

    const string pairspin_str = to_string(pairspin);
    iwlbuf tpdm_buf;
    const int nmo = dmap.size();
    const int npair = nmo*nmo;
    int npair_dirac;
    if(pairspin==AlphaBeta) {
      npair_dirac = nmo*nmo;
    }
    else {  // pairspin==AlphaAlpha || pairspin==BetaBeta
      npair_dirac = nmo*(nmo-1)/2;
    }
    //double *tpdm_arr = (double *)malloc(sizeof(double)*nmo*nmo*nmo*nmo);
    double *tpdm_arr;
    double **tpdm_arr2;
    if(pairspin==AlphaBeta) {
      tpdm_arr2 = new double*[npair];
      for(int i=0; i<npair; i++) {
        tpdm_arr2[i] = new double [npair];
        for(int j=0; j<npair; j++) {
          tpdm_arr2[i][j] = 0.0;
        }
      }
    }
    else {  // pairspin==AlphaAlpha || pairspin==BetaBeta
      tpdm_arr = new double [npair*(npair+1)/2];
      for(int i=0; i<npair*(npair+1)/2; i++) {
        tpdm_arr[i] = 0.0;
      }
    }

    //int *ioff_nonsymm = (int *)malloc(nmo*sizeof(int));
    int *ioff_nonsymm = new int [nmo];
    for(int i=0; i<nmo; i++){
      ioff_nonsymm[i] = i*nmo;
    }
    int *ioff_nonsymm_pair = new int [npair];
    for(int i=0; i<npair; i++) {
      ioff_nonsymm_pair[i] = i*npair;
    }
    const int ioff_dim = nmo*nmo;
    int *ioff = new int[ioff_dim];
    ioff[0] = 0;
    for (int i = 1; i<ioff_dim; i++) {
      ioff[i] = ioff[i-1] + i;
    }

    int TPDM_FILE;
    switch (pairspin) {
      case AlphaBeta:  TPDM_FILE = spinfree ? PSIF_MO_TPDM : PSIF_MO_AB_TPDM; break;
      case AlphaAlpha: TPDM_FILE = PSIF_MO_AA_TPDM; break;
      case BetaBeta:   TPDM_FILE = PSIF_MO_BB_TPDM; break;
      default: MPQC_ASSERT(false);
    }
    iwl_buf_init(&tpdm_buf,TPDM_FILE,0.0,1,1);
    if(pairspin==AlphaBeta) {
      iwl_buf_rd_all2(&tpdm_buf,tpdm_arr2,ioff_nonsymm,ioff_nonsymm,1,ioff_nonsymm_pair,1, stdout);
    }
    else {
      iwl_buf_rd_all( &tpdm_buf,tpdm_arr ,ioff_nonsymm,ioff_nonsymm,1,ioff,             0, NULL);
    }
    iwl_buf_close(&tpdm_buf,1);
    delete [] ioff_nonsymm;
    delete [] ioff_nonsymm_pair;
    delete [] ioff;

    const RefSCDimension tpdm_dim(new SCDimension(npair_dirac));
    if (kit == 0)
      kit = SCMatrixKit::default_matrixkit();
    RefSymmSCMatrix tpdm_mat = kit->symmmatrix(tpdm_dim);
    tpdm_mat.assign(0.0);

    if(pairspin==AlphaBeta) {
#if 1
      for(int i=0; i<nmo; i++) {
        const int ii = dmap[i];
        for(int j=0; j<nmo; j++) {
          const int jj = dmap[j];
          const int ind_ij = ordinary_INDEX(ii,jj,nmo);
          const int ind_ji = spinfree ? ordinary_INDEX(jj,ii,nmo) : 0;
          for(int k=0; k<=i; k++) {
            const int kk = dmap[k];
            int lmax;
            if(k==i) {
              lmax = j+1;
            }
            else {
              lmax = nmo;
            }
            for(int l=0; l<lmax; l++) {
              const int ll = dmap[l];
              //const int index = ordinary_INDEX(ordinary_INDEX(i,k,nmo),ordinary_INDEX(j,l,nmo),npair);
              const int ind_kl = ordinary_INDEX(kk,ll,nmo);
              const int ind_lk = spinfree ? ordinary_INDEX(ll,kk,nmo) : 0;
              const double value = tpdm_arr2[ordinary_INDEX(i,k,nmo)][ordinary_INDEX(j,l,nmo)];
              tpdm_mat->set_element(ind_ij,ind_kl,value);
              if (spinfree)
                tpdm_mat->set_element(ind_ji,ind_lk,value);
            }
          }
        }
      }
#else
      for(int i=0; i<nmo; i++) {
        for(int j=0; j<nmo; j++) {
          const int ind_ij = ordinary_INDEX(i,j,nmo);
          for(int k=0; k<nmo; k++) {
            for(int l=0; l<nmo; l++) {
              const int index = ordinary_INDEX(ordinary_INDEX(i,k,nmo),ordinary_INDEX(j,l,nmo),npair);
              //const int index = triang_half_INDEX(ordinary_INDEX(i,k,nmo),ordinary_INDEX(j,l,nmo));
              const int ind_kl = ordinary_INDEX(k,l,nmo);
              tpdm_mat->set_element(ind_ij,ind_kl,tpdm_arr[index]);
            }
          }
        }
      }
#endif
    }
    else {  // pairspin==AlphaAlpha || pairspin==BetaBeta
      for(int i=0; i<nmo; i++) {
        const int ii = dmap[i];
        for(int j=0; j<i; j++) {
          const int jj = dmap[j];
          const int ind_ij = lowerupper_index(ii,jj);
          const double pfac_ij = (ii > jj) ? 1.0 : -1.0;
          for(int k=0; k<=i; k++) {
            const int kk = dmap[k];
            for(int l=0; l<k; l++) {
              const int ll = dmap[l];
              const int ind_kl = lowerupper_index(kk,ll);
              const double pfac_kl = (kk > ll) ? 1.0 : -1.0;
              const int index = triang_half_INDEX(ordinary_INDEX(i,k,nmo),ordinary_INDEX(j,l,nmo));
              tpdm_mat->set_element(ind_ij,ind_kl,
                                    pfac_ij * pfac_kl * tpdm_arr[index]);
            }
          }
        }
      }
    }

    if(pairspin==AlphaBeta) {
      for(int i=0; i<npair; i++) {
        delete [] tpdm_arr2[i];
      }
      delete [] tpdm_arr2;
    }
    else {
      delete [] tpdm_arr;
    }

    // Psi convention for spin-free 2-RDM is a factor of 2 off
    if (spinfree)
      tpdm_mat.scale(2.0);

    return tpdm_mat;
  }

#if 0
  RefSymmSCMatrix PsiCorrWavefunction::onepdm(const SpinCase1 &spin) {
    string opdm_label_str;
    if(spin==Alpha) {
      opdm_label_str = "MO-basis Alpha OPDM";
      outfile = fopen("psiout_opdm_Alpha","w");
    }
    else { // spin==Beta
      opdm_label_str = "MO-basis Beta OPDM";
      outfile = fopen("psiout_opdm_Beta","w");
    }
    char *opdm_label=(char *)opdm_label_str.c_str();
    const int nmo = reference()->nmo();
    double *opdm_arr = new double [nmo*nmo];
    rdopdm(PSIF_MO_OPDM,opdm_label,opdm_arr,nmo*nmo,0,1,outfile);
    fclose(outfile);

    const Ref<SCMatrixKit> localkit = new LocalSCMatrixKit;
    const RefSCDimension opdm_dim(new SCDimension(nmo));
    RefSymmSCMatrix opdm_mat(localkit->symmmatrix(opdm_dim));
    for(int i=0; i<nmo; i++) {
      for(int j=0; j<=i; j++) {
        opdm_mat.set_element(i,j,opdm_arr[i*nmo+j]);
      }
    }

    delete [] opdm_arr;

    if(debug_>=DefaultPrintThresholds::mostN2) {
      opdm_mat.print(prepend_spincase(spin,"Psi opdm").c_str());
    }

    return(opdm_mat);
  }

  RefSymmSCMatrix PsiCorrWavefunction::onepdm(){
    outfile = fopen("psiout_opdm","w");
    const int nmo = reference()->nmo();
    //double *opdm_arr = (double *)malloc(sizeof(double)*nmo*nmo);
    double *opdm_arr = new double [nmo*nmo];
    char *opdm_label = "MO-basis OPDM";
    rdopdm(PSIF_MO_OPDM,opdm_label,opdm_arr,nmo*nmo,0,1,outfile);
    fclose(outfile);

    const Ref<SCMatrixKit> localkit = new LocalSCMatrixKit;
    const RefSCDimension opdm_dim(new SCDimension(nmo));
    RefSymmSCMatrix opdm_mat(localkit->symmmatrix(opdm_dim));
    for(int i=0; i<nmo; i++) {
      for(int j=0; j<=i; j++) {
        opdm_mat.set_element(i,j,opdm_arr[i*nmo+j]);
      }
    }

    //free(opdm_arr);
    delete [] opdm_arr;

    return(opdm_mat);
  }
#endif

  RefSymmSCMatrix PsiCorrWavefunction::mo_density(SpinCase1 spin){
    if ((spin == Beta || spin == AnySpinCase1) && !this->spin_polarized())
      return mo_density(Alpha);

    if (spin == AnySpinCase1 && this->spin_polarized())
      ProgrammingError("asked for any spin density but the density is spin-polarized");

    if (mo_density_[spin])
      return mo_density_[spin];

    // ensure that this has been computed
    { const double energy = this->value(); }

    const std::vector<unsigned int> dmap = map_density_to_sb();
    mo_density_[spin] = detail::rdopdm(spin,
                                       this->reference()->mopi(),
                                       dmap,
                                       this->basis_matrixkit());//read rdm

    if(debug_>=DefaultPrintThresholds::mostN2) {
      mo_density_[spin].print(prepend_spincase(spin,"Psi opdm").c_str());
    }

    return mo_density_[spin];
  }

#if 0
  RefSymmSCMatrix PsiCorrWavefunction::twopdm(){
    iwlbuf tpdm_buf;
    const int nmo = reference()->nmo();
    double *tpdm_arr = (double *)malloc(sizeof(double)*nmo*nmo*nmo*nmo);
    for(int i=0; i<nmo*nmo*nmo*nmo; i++) {
      tpdm_arr[i] = 0.0;
    }
    int *ioff_nonsymm = (int *)malloc(nmo*sizeof(int));
    for(int i=0; i<nmo; i++){
      ioff_nonsymm[i] = i*nmo;
    }
    const int ioff_dim = nmo*nmo;
    int *ioff = new int[ioff_dim];
    ioff[0] = 0;
    for (int i = 1; i<ioff_dim; i++) {
      ioff[i] = ioff[i-1] + i;
    }
    iwl_buf_init(&tpdm_buf,PSIF_MO_TPDM,0.0,1,1);
    iwl_buf_rd_all(&tpdm_buf,tpdm_arr,ioff_nonsymm,ioff_nonsymm,1,ioff,0,NULL);
    iwl_buf_close(&tpdm_buf,1);
    free(ioff_nonsymm);
    //ExEnv::out0() << "tpdm_arr:" << endl;
    //for(int i=0; i<nmo*nmo*nmo*nmo; i++) {
    //  ExEnv::out0() << setprecision(12) << tpdm_arr[i] << endl;
    //}

    Ref<SCMatrixKit> localkit = new LocalSCMatrixKit;
    RefSCDimension tpdm_dim(new SCDimension(ordinary_INDEX(nmo-1,nmo-1,nmo)+1));
    RefSymmSCMatrix tpdm_mat = localkit->symmmatrix(tpdm_dim);
    //for(int i=0; i<nmo; i++) {
    //  for(int j=0; j<nmo; j++) {
    //    for(int k=0; k<nmo; k++) {
    //      for(int m=0; m<nmo; m++) {
    //        tpdm_mat.set_element(ordinary_INDEX(i,j,nmo), ordinary_INDEX(k,m,nmo),
    //                             tpdm_arr[tpdm_index(i,j,k,m,nmo)]);
    //      }
    //    }
    //  }
    //}
    //int index = 0;
    //for(int halfind1=0; halfind1<tpdm_dim.n(); halfind1++) {
    //  for(int halfind2=0; halfind2<tpdm_dim.n(); halfind2++, index++) {
    //    tpdm_mat.set_element(halfind1,halfind2,tpdm_arr[index]);
    //  }
    //}
    int index = 0;
    for(int i=0; i<nmo; i++) {
      for(int j=0; j<nmo; j++) {
        for(int k=0; k<=i; k++) {
          int lmax;
          if(k==i) {
            lmax = j+1;
          }
          else {
            lmax = nmo;
          }
          for(int l=0; l<lmax; l++, index++) {
            //ExEnv::out0() << "i=" << i << "  j=" << j << "  k=" << k << "  l=" << l << endl;
            tpdm_mat.set_element(ordinary_INDEX(i,j,nmo), ordinary_INDEX(k,l,nmo),
                                 2.0 * tpdm_arr[index]);
          }
        }
      }
    }

    //FILE *testf = fopen("tpdm_test","w");
    //print_twopdm_arr(testf,tpdm_arr,1.0e-6);
    //fclose(testf);

    free(tpdm_arr);

    return(tpdm_mat);
  }

  RefSymmSCMatrix PsiCorrWavefunction::twopdm_dirac() {
    iwlbuf tpdm_buf;
    const int nmo = reference()->nmo();
    const int npair = nmo*nmo;
    //double *tpdm_arr = (double *)malloc(sizeof(double)*nmo*nmo*nmo*nmo);
#if 1
    double *tpdm_arr = new double [npair*(npair+1)/2];
    for(int i=0; i<npair*(npair+1)/2; i++) {
      tpdm_arr[i] = 0.0;
    }
#endif
#if 0
    double *tpdm_arr = new double [nmo*nmo*nmo*nmo];
    for(int i=0; i<nmo*nmo*nmo*nmo; i++) {
      tpdm_arr[i] = 0.0;
    }
#endif

    //int *ioff_nonsymm = (int *)malloc(nmo*sizeof(int));
    int *ioff_nonsymm = new int [nmo];
    for(int i=0; i<nmo; i++){
      ioff_nonsymm[i] = i*nmo;
    }
    //const int ioff_dim = nmo*nmo;
    const int ioff_dim = npair;
    int *ioff = new int[ioff_dim];
    ioff[0] = 0;
    for (int i = 1; i<ioff_dim; i++) {
      ioff[i] = ioff[i-1] + i;
    }
    iwl_buf_init(&tpdm_buf,PSIF_MO_TPDM,0.0,1,1);
    iwl_buf_rd_all(&tpdm_buf,tpdm_arr,ioff_nonsymm,ioff_nonsymm,1,ioff,0,NULL);
    iwl_buf_close(&tpdm_buf,1);
    delete [] ioff_nonsymm;
    delete [] ioff;

    const int npair_dirac = nmo*nmo;
    const Ref<SCMatrixKit> localkit = new LocalSCMatrixKit;
    const RefSCDimension tpdm_dim(new SCDimension(npair_dirac));
    RefSymmSCMatrix tpdm_mat = localkit->symmmatrix(tpdm_dim);
    tpdm_mat.assign(0.0);

#if 0
    int index = 0;
    for(int i=0; i<nmo; i++) {
      for(int j=0; j<nmo; j++) {
        for(int k=0; k<=i; k++) {
          int lmax;
          if(k==i) {
            lmax = j+1;
          }
          else {
            lmax = nmo;
          }
          const int ind_ik = ordinary_INDEX(i,k,nmo);
          const int ind_ki = ordinary_INDEX(k,i,nmo);
          for(int l=0; l<lmax; l++, index++) {
            //ExEnv::out0() << "i=" << i << "  j=" << j << "  k=" << k << "  l=" << l << endl;
            //tpdm_mat.set_element(ordinary_INDEX(i,j,nmo), ordinary_INDEX(k,l,nmo),
            //                     2.0 * tpdm_arr[index]);
            const int ind_jl = ordinary_INDEX(j,l,nmo);
            const int ind_lj = ordinary_INDEX(l,j,nmo);
            tpdm_mat.set_element(ind_ik, ind_jl,
                                 2.0 * tpdm_arr[index]);
            if((k!=i) || (l!=j)) {
              tpdm_mat.set_element(ind_ki, ind_lj,
                                   2.0 * tpdm_arr[index]);
            }
          }
        }
      }
    }
#endif
    for(int i=0; i<nmo; i++) {
      for(int j=0; j<nmo; j++) {
        const int ind_ij = ordinary_INDEX(i,j,nmo);
        for(int k=0; k<=i; k++) {
          int lmax;
          if(k==i) {
            lmax = j+1;
          }
          else {
            lmax = nmo;
          }
          for(int l=0; l<lmax; l++) {
            const int index = triang_half_INDEX(ordinary_INDEX(i,k,nmo),ordinary_INDEX(j,l,nmo));
            const int ind_kl = ordinary_INDEX(k,l,nmo);
            tpdm_mat.set_element(ind_ij,ind_kl,
                                 2.0 * tpdm_arr[index]);
          }
        }
      }
    }

    delete [] tpdm_arr;

    return(tpdm_mat);
  }

  RefSymmSCMatrix PsiCorrWavefunction::twopdm_dirac_from_components() {
    const Ref<SCMatrixKit> localkit = new LocalSCMatrixKit;
    const int nmo = reference()->nmo();
    const int npair = nmo*nmo;
    const RefSCDimension pair_dim = new SCDimension(npair);
    RefSymmSCMatrix tpdm = localkit->symmmatrix(pair_dim);
    tpdm->assign(0.0);
    for(int spincase2=0; spincase2<NSpinCases2; spincase2++) {
      SpinCase2 pairspin = static_cast<SpinCase2>(spincase2);
      //const int npair = (pairspin==AlphaBeta) ? nmo*nmo : nmo*(nmo-1)/2;
      RefSymmSCMatrix tpdm_sc = twopdm_dirac(pairspin);
      if(pairspin==AlphaBeta) {
        for(int p=0; p<nmo; p++) {
          for(int q=0; q<nmo; q++) {
            const int ind_pq = ordinary_INDEX(p,q,nmo);
            const int ind_qp = ordinary_INDEX(q,p,nmo);
            for(int r=0; r<=p; r++) {
              const int smax = (r==p) ? q+1 : nmo;
              for(int s=0; s<smax; s++) {
              const int ind_rs = ordinary_INDEX(r,s,nmo);
                const int ind_sr = ordinary_INDEX(s,r,nmo);
                tpdm->accumulate_element(ind_pq,ind_rs,tpdm_sc->get_element(ind_pq,ind_rs)+tpdm_sc->get_element(ind_qp,ind_sr));
              }
            }
          }
        }
      }
      else {
        for(int p=0; p<nmo; p++) {
          for(int q=0; q<nmo; q++) {
            if(p!=q) {
              const int ind_pq = ordinary_INDEX(p,q,nmo);
              const int ind_pq_lu = lowerupper_index(p,q);
              for(int r=0; r<=p; r++) {
                const int smax = (r==p) ? q+1 : nmo;
                for(int s=0; s<smax; s++) {
                  if(r!=s) {
                    const int ind_rs = ordinary_INDEX(r,s,nmo);
                    const int ind_rs_lu = lowerupper_index(r,s);
                    tpdm->accumulate_element(ind_pq,ind_rs,tpdm_sc->get_element(ind_pq_lu,ind_rs_lu));
                  }
                }
              }
            }
          }
        }
      }
    }
    return(tpdm);
  }
#endif

  RefSymmSCMatrix PsiCorrWavefunction::twopdm_dirac(const SpinCase2 &pairspin) {
    // map orbitals to symmetry-blocked orbitals
    const std::vector<unsigned int> dmap = map_density_to_sb();
    RefSymmSCMatrix result = detail::rdtpdm(pairspin, dmap);
    return result;
  }

  RefSymmSCMatrix PsiCorrWavefunction::twopdm_dirac() {
    // map orbitals to symmetry-blocked orbitals
    const std::vector<unsigned int> dmap = map_density_to_sb();
    RefSymmSCMatrix result = detail::rdtpdm(AlphaBeta, dmap, true);
    return result;
  }

  void PsiCorrWavefunction::print_onepdm_vec(FILE *output,const RefSCVector &opdm,double TOL) {
    int nmo = reference()->nmo();
    for(int i=0; i<nmo; i++){
      for(int j=0; j<nmo; j++){
        int ind = triang_half_INDEX(i,j);
        if(fabs(opdm.get_element(ind))>=TOL) {
          fprintf(output," %i %i [%i] %.10f\n",i,j,ind,opdm.get_element(ind));
        }
      }
    }
  }

  void PsiCorrWavefunction::print_onepdm_mat(FILE *output,const RefSymmSCMatrix &opdm,double TOL) {
    int nmo = reference()->nmo();
    for(int i=0; i<nmo; i++){
      for(int j=0; j<nmo; j++){
        if(fabs(opdm.get_element(i,j))>=TOL) {
          fprintf(output," %i %i %.10f\n",i,j,opdm.get_element(i,j));
        }
      }
    }
  }

  void PsiCorrWavefunction::print_twopdm_mat(FILE *output,const RefSymmSCMatrix &tpdm, double TOL) {
    int nmo = reference()->nmo();
    for(int i=0; i<nmo; i++){
      for(int j=0; j<nmo; j++){
        for(int k=0; k<nmo; k++){
          for(int l=0; l<nmo; l++){
            int ind_half1 = ordinary_INDEX(i,j,nmo);
            int ind_half2 = ordinary_INDEX(k,l,nmo);
            if(fabs(tpdm.get_element(ind_half1,ind_half2))>=TOL) {
              fprintf(output," %i %i %i %i [%i] [%i] %.10f\n",i,j,k,l,ind_half1,ind_half2,0.5*tpdm.get_element(ind_half1,ind_half2));
            }
          }
        }
      }
    }
  }

  void PsiCorrWavefunction::print_twopdm_arr(FILE *output,double *tpdm,double TOL){
    int nmo = reference()->nmo();
    for(int i=0; i<nmo; i++){
      for(int j=0; j<nmo; j++){
        for(int k=0; k<nmo; k++){
          for(int l=0; l<nmo; l++){
            int ind=tpdm_index(i,j,k,l,nmo);
            if(fabs(tpdm[ind])>=TOL){
              fprintf(output," %i %i %i %i [%i][%i][[%i]] %.10f\n",i,j,k,l,ordinary_INDEX(i,j,nmo),ordinary_INDEX(k,l,nmo),ind,tpdm[ind]);
            }
          }
        }
      }
    }
  }

//////////////////////////////////////////////////////////////////////////

}

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:

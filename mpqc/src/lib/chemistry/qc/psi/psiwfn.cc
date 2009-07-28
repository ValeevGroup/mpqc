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
#include <numeric>

#include <psifiles.h>

#include <util/keyval/keyval.h>
#include <util/misc/formio.h>
#include <util/class/scexception.h>
#include <util/state/stateio.h>
#include <math/scmat/local.h>
#include <math/scmat/matrix.h>
#include <math/scmat/blocked.h>
#include <math/symmetry/pointgrp.h>
#include <chemistry/molecule/molecule.h>
#include <chemistry/qc/basis/petite.h>
#include <chemistry/qc/basis/shellrot.h>
#include <chemistry/qc/mbptr12/mbptr12.h>
#include <chemistry/qc/mbptr12/spin.h>
#include <chemistry/qc/mbptr12/print.h>
#include <chemistry/qc/mbptr12/r12int_eval.h>
#include <chemistry/qc/mbptr12/fixedcoefficient.h>

#include <chemistry/qc/mbptr12/creator.h>
#include <chemistry/qc/mbptr12/container.h>
#include <chemistry/qc/mbptr12/compute_tbint_tensor.h>

#include <chemistry/qc/psi/psiwfn.h>
#include <chemistry/qc/psi/psiinput.timpl.h>
#include <chemistry/qc/psi/psiqtorder.h>
#include <chemistry/qc/mbptr12/pairiter.h>
#include <chemistry/qc/scf/clscf.h>

#include <libiwl/iwl.h>
#include <libciomr/libciomr.h>

#define TEST_V 0

/// for indexing of density matrices.
//#define triang_half_INDEX(i,j) ((i>j) ? (ioff[(i)]+(j)) : (ioff[(j)]+(i)))
//#define ordinary_INDEX(i,j,coldim) ((i)*(coldim)+(j))

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

  namespace {

    int triang_half_INDEX_ordered(int i, int j) {
      return(i*(i+1)/2+j);
    }

    int triang_half_INDEX(int i, int j) {
      return((i>j) ? triang_half_INDEX_ordered(i,j) : triang_half_INDEX_ordered(j,i));
    }

    int ordinary_INDEX(int i, int j, int coldim) {
      return(i*coldim+j);
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

    int lowerupper_index(int p, int q);

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

    int lowertriang_index(int p,int q) {
      if(q>=p){
        throw ProgrammingError("lowertriang_index(p,q) -- q must be smaller than p.",__FILE__,__LINE__);
      }
      int index=p*(p+1)/2+q-p;
      return(index);
    }

    int lowerupper_index(int p, int q) {
      if(p>q) {
        return(lowertriang_index(p,q));
      }
      else if(q>p) {
        return(lowertriang_index(q,p));
      }
      else {
        throw ProgrammingError("lowerupper_index(p,q) -- p and q are not allowed to be equal.",__FILE__,__LINE__);
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

  //void init_ioff(void)
  //{
  //   int i;

  //   /* set offsets for ij-type canonical ordering */
  //   ioff = (int *) malloc (IOFF_MAX * sizeof(int)) ;
  //   ioff[0] = 0;
  //   for (i = 1; i < IOFF_MAX ; i++) {
  //      ioff[i] = ioff[i-1] + i;
  //      }
  //}


  //////////////////////////////////////////////////////////////////////////

  static ClassDesc PsiWavefunction_cd(typeid(PsiWavefunction),
                                      "PsiWavefunction", 2,
                                      "public Wavefunction", 0, 0, 0);

  PsiWavefunction::PsiWavefunction(const Ref<KeyVal>&keyval) :
    Wavefunction(keyval) {
    exenv_ << keyval->describedclassvalue("psienv");
    if (exenv_.null()) {
      exenv_ = new PsiExEnv;
    }

    nirrep_ = molecule()->point_group()->char_table().order();

    size_t bytes = keyval->sizevalue("memory");
    if (bytes <= 2000000)
      bytes = 2000000;
    {
      std::ostringstream oss;
      oss << "(" << bytes << " B)";
      const int nchar = oss.str().size();
      memory_str_ = new char[nchar+1];
      strcpy(memory_str_,oss.str().c_str());
    }
    memory_ = bytes;
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
    psiinput->write_keyword("psi:memory", memory_str_);
    psiinput->write_keyword("psi:convergence", conv);
    psiinput->write_keyword("cints:cutoff", conv+5);
    psiinput->begin_section("input");
    psiinput->write_keyword("no_reorient", "true");
    psiinput->write_keyword("keep_ref_frame", "true");
    psiinput->write_basis(basis());
    if (basis()->has_pure())
      psiinput->write_keyword("puream", "true");
    psiinput->write_geom(molecule());
    psiinput->end_section();
    if (electric_field()) {
      double E[3];  for(int xyz=0; xyz<3; ++xyz) { E[xyz] = electric_field().get_element(xyz); }
      psiinput->write_keyword_array("psi:efield", 3, E);
    }
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

  double
  PsiWavefunction::nuclear_repulsion_energy() const
  {
    return exenv()->chkpt().rd_enuc();
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
        // get energy to make sure that it's computed.
        const double energy = guess_wfn_->value();
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
      input->write_keyword_array("psi:docc", docc_);
    else {
      input->write_keyword("psi:reset_occupations", true);
      input->write_keyword("psi:hcore_guess", "new");
    }
    input->write_keyword("scf:maxiter", maxiter);
    input->write_keyword("scf:convergence",convergence);
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
      input->write_keyword_array("psi:docc", docc_);
    if (!socc_.empty())
      input->write_keyword_array("psi:socc", socc_);
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
      input->write_keyword_array("psi:docc", docc_);
    if (!socc_.empty())
      input->write_keyword_array("psi:socc", socc_);
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
    input->write_keyword("psi:tolerance", convergence + 5);
    input->write_keyword("psi:freeze_core", static_cast<int>(nfzc_));
    input->write_keyword("psi:freeze_virt", static_cast<int>(nfzv_));
  }

  int PsiCorrWavefunction::nelectron() {
    return reference()->nelectron();
  }

  const Ref<OrbitalSpace>&
  PsiCorrWavefunction::occ_act_sb(SpinCase1 spin) {
    if (occ_act_sb_[spin].nonnull())
      return occ_act_sb_[spin];
    if (reference_->reftype() == PsiSCF::rhf && spin==Beta)
      return occ_act_sb(Alpha);

    const unsigned int nmo = reference_->nmo();
    const std::vector<unsigned int>& mopi = reference_->mopi();
    const std::vector<unsigned int>& occpi = reference_->occpi(spin);
    const int nirreps = occpi.size();
    std::vector<bool> occ_mask(nmo, false);
    for(int irrep=0, irrep_offset=0; irrep<nirreps; ++irrep) {
      for(int i=0; i<occpi[irrep]; ++i) {
        occ_mask[i + irrep_offset] = true;
      }
      irrep_offset += mopi[irrep];
    }

    const std::string id(spin==Alpha ? "I" : "i");
    Ref<OrbitalSpace> orbs_sb = new OrbitalSpace("", "", reference_->coefs(spin), basis(), integral(),
                                                 reference_->evals(spin), 0, 0, OrbitalSpace::symmetry);
    Ref<OrbitalSpace> occ_sb = new MaskedOrbitalSpace("", "", orbs_sb, occ_mask);
    occ_act_sb_[spin] = new OrbitalSpace(id, prepend_spincase(spin,"active occupied MOs (Psi3)"),
                                         occ_sb->coefs(), occ_sb->basis(), occ_sb->integral(),
                                         occ_sb->evals(), nfzc_, 0, OrbitalSpace::symmetry);

    return occ_act_sb_[spin];
  }

  const Ref<OrbitalSpace>&
  PsiCorrWavefunction::vir_act_sb(SpinCase1 spin) {
    if (vir_act_sb_[spin].nonnull())
      return vir_act_sb_[spin];
    if (reference_->reftype() == PsiSCF::rhf && spin==Beta)
      return vir_act_sb(Alpha);

    const unsigned int nmo = reference_->nmo();
    const std::vector<unsigned int>& mopi = reference_->mopi();
    const std::vector<unsigned int>& uoccpi = reference_->uoccpi(spin);
    const int nirreps = uoccpi.size();
    std::vector<bool> uocc_mask(nmo, false);
    for(int irrep=0, irrep_offset=0; irrep<nirreps; ++irrep) {
      const unsigned int nocc = mopi[irrep] - uoccpi[irrep];
      for(int i=0; i<uoccpi[irrep]; ++i) {
        uocc_mask[i + irrep_offset + nocc] = true;
      }
      irrep_offset += mopi[irrep];
    }

    const std::string id(spin==Alpha ? "A" : "a");
    Ref<OrbitalSpace> orbs_sb = new OrbitalSpace("", "", reference_->coefs(spin), basis(), integral(),
                                                 reference_->evals(spin), 0, 0, OrbitalSpace::symmetry);
    Ref<OrbitalSpace> uocc_sb = new MaskedOrbitalSpace("", "", orbs_sb, uocc_mask);
    vir_act_sb_[spin] = new OrbitalSpace(id, prepend_spincase(spin, "active virtual MOs (Psi3)"),
                                         uocc_sb->coefs(), uocc_sb->basis(), uocc_sb->integral(),
                                         uocc_sb->evals(), 0, nfzv_, OrbitalSpace::symmetry);

    return vir_act_sb_[spin];
  }

  const Ref<OrbitalSpace>&  PsiCorrWavefunction::orbs_sb(SpinCase1 spin) {
    if(orbitals_sb_[spin].nonnull()) {
      return(orbitals_sb_[spin]);
    }
    if (reference_->reftype() == PsiSCF::rhf && spin==Beta) {
      return(orbitals_sb_[Alpha]);
    }

    const int nmo = reference_->nmo();
    const int nocc = reference_->nocc(spin);

    const std::string id(spin==Alpha ? "P" : "p");

    orbitals_sb_[spin] = new OrbitalSpace(id,prepend_spincase(spin,"MOs (Psi3)"),
                                                     reference_->coefs(spin),basis(),integral(),
                                                     reference_->evals(spin),0,0,OrbitalSpace::symmetry);

    return(orbitals_sb_[spin]);
  }

  unsigned int
  PsiCorrWavefunction::nfzc() const { return nfzc_; }
  unsigned int
  PsiCorrWavefunction::nfzv() const { return nfzv_; }

  const std::vector<unsigned int>&
  PsiCorrWavefunction::frozen_docc() const {
    if (frozen_docc_.empty()) {
      frozen_docc_.resize(nirrep_);
      int* frzcpi = exenv()->chkpt().rd_frzcpi();
      std::copy(frzcpi,frzcpi+nirrep_,frozen_docc_.begin());
      psi::Chkpt::free(frzcpi);
    }
    return frozen_docc_;
  }

  const std::vector<unsigned int>&
  PsiCorrWavefunction::frozen_uocc() const {
    if (frozen_uocc_.empty()) {
      frozen_uocc_.resize(nirrep_);
      int* frzvpi = exenv()->chkpt().rd_frzvpi();
      std::copy(frzvpi,frzvpi+nirrep_,frozen_uocc_.begin());
      psi::Chkpt::free(frzvpi);
    }
    return frozen_uocc_;
  }

  namespace psiwfnlocal {

    std::vector<unsigned int> vector_binop(const std::vector<unsigned int> &term1,const std::vector<unsigned int> &term2,
                                         const std::string &binop) {
      if(term1.size() != term2.size()){
        throw ProgrammingError("psiwfnlocal::vector_binop -- vectors term1 and term2 don't have the same length.",__FILE__,__LINE__);
      }
      std::vector<unsigned int> result(term1.size());
      std::vector<unsigned int>::iterator result_it = result.begin();
      if(binop==string("plus")) {
        for(std::vector<unsigned int>::const_iterator term1_it = term1.begin(), term2_it = term2.begin();
            term1_it != term1.end() && term2_it != term2.end() && result_it != result.end();
            term1_it++, term2_it++, result_it++) {
          (*result_it) = (*term1_it) + (*term2_it);
        }
      }
      else if(binop=="minus") {
        for(std::vector<unsigned int>::const_iterator term1_it = term1.begin(), term2_it = term2.begin();
            term1_it != term1.end() && term2_it != term2.end() && result_it != result.end();
            term1_it++, term2_it++, result_it++) {
          (*result_it) = (*term1_it) - (*term2_it);
        }
      }
      else {
        throw ProgrammingError("psiwfnlocal::vector_binop -- this binop is not implemented.",__FILE__,__LINE__);
      }
      return(result);
    }

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
      docc_active = psiwfnlocal::vector_binop(occ_act_beta,frozen_docc(),string("minus"));
    }
    else {
      // docc_active = occ_act_alpha - frozen_docc()
      docc_active = psiwfnlocal::vector_binop(occ_act_alpha,frozen_docc(),string("minus"));
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
      singocc = psiwfnlocal::vector_binop(occ_act_alpha,occ_act_beta,string("minus"));
    }
    else {
      // singocc = occ_act_beta - occ_act_alpha
      singocc = psiwfnlocal::vector_binop(occ_act_beta,occ_act_alpha,string("minus"));
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
      uocc_active = psiwfnlocal::vector_binop(uocc_act_beta,frozen_uocc(),string("minus"));
    }
    else {
      // uocc_active = uocc_act_alpha - frozen_uocc()
      uocc_active = psiwfnlocal::vector_binop(uocc_act_alpha,frozen_uocc(),string("minus"));
    }

    return(uocc_active);
  }

  double
  PsiCorrWavefunction::reference_energy()
  {
    return exenv()->chkpt().rd_eref();
  }

  extern "C" {
    FILE *outfile;
  }

  /// A modification of Psi3's iwl_rdone. For documentation see that of iwl_rdone.
  int rdopdm(int itap,char *label, double *ints, int ntri, int erase,
             int printflg, FILE *outfile) {
    int nmo;
    double TOL=1.0e-6;

    psio_open(itap, PSIO_OPEN_OLD);
    psio_read_entry(itap, label, (char *) ints, ntri*sizeof(double));
    psio_close(itap, !erase);

    if (printflg) {
      nmo=round(sqrt((double)ntri));
      for(int i=0; i<nmo; i++) {
        for(int j=0; j<nmo; j++) {
          if(fabs(ints[i*nmo+j])>=TOL){
            fprintf(outfile,"%5i%5i%12.7f\n",i+1,j+1,ints[i*nmo+j]);
          }
        }
      }
    }
    return(1);
  }

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
      opdm_mat.print(prepend_spincase(spin,"opdm").c_str());
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
    outfile = fopen("psiout_tpdm","w");
    iwl_buf_init(&tpdm_buf,PSIF_MO_TPDM,0.0,1,1);
    iwl_buf_rd_all(&tpdm_buf,tpdm_arr,ioff_nonsymm,ioff_nonsymm,1,ioff,1,outfile);
    iwl_buf_close(&tpdm_buf,1);
    fclose(outfile);
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
    outfile = fopen("psiout_tpdm","w");
    iwl_buf_init(&tpdm_buf,PSIF_MO_TPDM,0.0,1,1);
    iwl_buf_rd_all(&tpdm_buf,tpdm_arr,ioff_nonsymm,ioff_nonsymm,1,ioff,1,outfile);
    iwl_buf_close(&tpdm_buf,1);
    fclose(outfile);
    //free(ioff_nonsymm);
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

  RefSymmSCMatrix PsiCorrWavefunction::twopdm_dirac(const SpinCase2 &pairspin) {
    string pairspin_str = (pairspin==AlphaBeta) ? "AlphaBeta" : (pairspin==AlphaAlpha) ? "AlphaAlpha" : "BetaBeta";
    iwlbuf tpdm_buf;
    const int nmo = reference()->nmo();
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
    string outputname = string("psiout_tpdm_")+pairspin_str;
    const int TPDM_FILE = (pairspin==AlphaBeta) ? PSIF_MO_AB_TPDM : (pairspin==AlphaAlpha) ? PSIF_MO_AA_TPDM : PSIF_MO_BB_TPDM;
    //outfile = fopen("psiout_tpdm","w");
    outfile = fopen(outputname.c_str(),"w");
    iwl_buf_init(&tpdm_buf,TPDM_FILE,0.0,1,1);
    if(pairspin==AlphaBeta) {
      iwl_buf_rd_all2(&tpdm_buf,tpdm_arr2,ioff_nonsymm,ioff_nonsymm,1,ioff_nonsymm_pair,1,outfile);
    }
    else {
      iwl_buf_rd_all( &tpdm_buf,tpdm_arr ,ioff_nonsymm,ioff_nonsymm,1,ioff,             1,outfile);
    }
    iwl_buf_close(&tpdm_buf,1);
    fclose(outfile);
    //free(ioff_nonsymm);
    delete [] ioff_nonsymm;
    delete [] ioff_nonsymm_pair;
    delete [] ioff;

    const Ref<SCMatrixKit> localkit = new LocalSCMatrixKit;
    const RefSCDimension tpdm_dim(new SCDimension(npair_dirac));
    RefSymmSCMatrix tpdm_mat = localkit->symmmatrix(tpdm_dim);
    tpdm_mat.assign(0.0);

    if(pairspin==AlphaBeta) {
#if 1
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
              //const int index = ordinary_INDEX(ordinary_INDEX(i,k,nmo),ordinary_INDEX(j,l,nmo),npair);
              const int ind_kl = ordinary_INDEX(k,l,nmo);
              tpdm_mat->set_element(ind_ij,ind_kl,tpdm_arr2[ordinary_INDEX(i,k,nmo)][ordinary_INDEX(j,l,nmo)]);
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
        for(int j=0; j<i; j++) {
          const int ind_ij = lowerupper_index(i,j);
          for(int k=0; k<=i; k++) {
            for(int l=0; l<k; l++) {
              const int index = triang_half_INDEX(ordinary_INDEX(i,k,nmo),ordinary_INDEX(j,l,nmo));
              const int ind_kl = lowerupper_index(k,l);
              tpdm_mat->set_element(ind_ij,ind_kl,tpdm_arr[index]);
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

    return(tpdm_mat);
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

  //static ClassDesc PsiCorrWavefunction_PT2R12_cd(typeid(PsiCorrWavefunction_PT2R12),"PsiCorrWavefunction_PT2R12",
  //                                               1,"public PsiCorrWavefunction",
  //                                               0,create<PsiCorrWavefunction_PT2R12>,create<PsiCorrWavefunction_PT2R12>);
  static ClassDesc PsiCorrWavefunction_PT2R12_cd(typeid(PsiCorrWavefunction_PT2R12),"PsiCorrWavefunction_PT2R12",
                                                 1,"public PsiCorrWavefunction",0,0,0);

  PsiCorrWavefunction_PT2R12::PsiCorrWavefunction_PT2R12(const Ref<KeyVal> &keyval)
    : PsiCorrWavefunction(keyval) {
    reference_mpqc_ << keyval->describedclassvalue("reference_mpqc");
    if (reference_mpqc_.null()) {
        ExEnv::err0() << "PsiCorrWavefunction_PT2R12::PsiCorrWavefunction_PT2R12: no MPQC reference wavefunction" << endl;
        abort();
      }
    copy_orthog_info(reference_mpqc_);
    bool spinadapted = false;
    if (dynamic_cast<CLSCF*>(reference_mpqc_.pointer()) != 0)
      // Default is to use spin-adapted algorithm
      spinadapted = keyval->booleanvalue("spinadapted",KeyValValueboolean((int)true));

    r12evalinfo_ = new R12IntEvalInfo(keyval,this,reference_mpqc_,nfzc_, nfzv_, spinadapted, true);
    // make sure that r12evalinfo is initialized
    r12evalinfo_->initialize();
    if(r12evalinfo_->r12tech()->corrfactor()->geminaldescriptor()->type()!=string("invalid")) {
      r12eval_ = new R12IntEval(r12evalinfo_);
      r12eval_->set_memory(memory_);
    }
    // test that Psi3 is compatible to the integral factory used by mbptr12_ object
    if (r12evalinfo_->integral()->cartesian_ordering() != PsiWavefunction::cartesian_ordering()) {
      throw InputError("PsiCorrWavefunction_PT2R12 -- ordering of functions in shells differs between the provided mbpt2r12 and Psi3. Use different Integral factory to construct mbpt2r12.",__FILE__,__LINE__);
    }

    std::string opdm_type_str = keyval->stringvalue("opdm_type",KeyValValuestring("correlated"));
    if(opdm_type_str == "correlated") {
      ExEnv::out0() << "PsiCorrWavefunction_PT2R12 -- Using correlated onepdm." << endl;
      opdm_type_ = correlated;
    }
    else if(opdm_type_str == "HF") {
      ExEnv::out0() << "PsiCorrWavefunction_PT2R12 -- Using HF onepdm." << endl;
      opdm_type_ = HF;
    }
    else {
      throw InputError("PsiCorrWavefunction_PT2R12::PsiCorrWavefunction_PT2R12 -- opdm_type not known.",__FILE__,__LINE__);
    }

    tpdm_from_opdms_ = keyval->booleanvalue("tpdm_from_opdms",KeyValValueboolean((int)false));
    if(tpdm_from_opdms_==false && opdm_type_==HF) {
      throw InputError("PsiCorrWavefunction_PT2R12::PsiCorrWavefunction_PT2R12(const Ref<KeyVal> &) -- HF opdm usage only makes sense if tpdms are constructed from opdms.",__FILE__,__LINE__);
    }

    exactgamma_phi_twoelec_ = keyval->booleanvalue("exactgamma_phi_twoelec",KeyValValueboolean((int)true));

    if(keyval->exists("debug")) {
      debug_ = keyval->intvalue("debug");
    }
    else {
      debug_ = 0;
    }
    if(r12evalinfo_->r12tech()->corrfactor()->geminaldescriptor()->type()!=string("invalid")) {
      r12eval_->set_debug(debug_);
    }
  }

  PsiCorrWavefunction_PT2R12::PsiCorrWavefunction_PT2R12(StateIn &s)
    : PsiCorrWavefunction(s) {
    reference_mpqc_ << SavableState::restore_state(s);
    r12evalinfo_ << SavableState::restore_state(s);
    // make sure that r12evalinfo is initialized
    r12evalinfo_->initialize();
    if(r12evalinfo_->r12tech()->corrfactor()->geminaldescriptor()->type()!=string("invalid")) {
      r12eval_ << SavableState::restore_state(s);
    }
    int o; s.get(o); opdm_type_ = (onepdm_type)o;
    int t; s.get(t); tpdm_from_opdms_ = (bool)t;
    int e; s.get(e); exactgamma_phi_twoelec_ = (bool)e;
    s.get(debug_);
    if(r12evalinfo_->r12tech()->corrfactor()->geminaldescriptor()->type()!=string("invalid")) {
      r12eval_->set_debug(debug_);
    }
  }

  PsiCorrWavefunction_PT2R12::~PsiCorrWavefunction_PT2R12() {}

  void PsiCorrWavefunction_PT2R12::save_data_state(StateOut &s) {
    PsiCorrWavefunction::save_data_state(s);
    SavableState::save_state(reference_mpqc_, s);
    SavableState::save_state(r12evalinfo_, s);
    if(r12evalinfo_->r12tech()->corrfactor()->geminaldescriptor()->type()!=string("invalid")) {
      SavableState::save_state(r12eval_, s);
    }
    s.put((int)opdm_type_);
    s.put((int)tpdm_from_opdms_);
    s.put((int)exactgamma_phi_twoelec_);
    s.put(debug_);
  }

  void PsiCorrWavefunction_PT2R12::write_input(int convergence) {
    PsiCorrWavefunction::write_input(convergence);
  }

  RefSymmSCMatrix PsiCorrWavefunction_PT2R12::hcore_mo() {
    Ref<SCMatrixKit> localkit = new LocalSCMatrixKit;
    Ref<OrbitalSpace> space = r12evalinfo_->refinfo()->orbs(Alpha);
    RefSCMatrix coeffs = space->coefs();
    RefSCDimension nmodim = coeffs.rowdim();
    RefSCDimension naodim = coeffs.coldim();
    int nmo = nmodim.n();
    int nao = naodim.n();
    //ExEnv::out0() << "nao = " << nao << "    nmo = " << nmo << endl;
    RefSCMatrix coeffs_nb = localkit->matrix(naodim,nmodim);
    //RefSCMatrix coeffs_nb = localkit->matrix(nmodim,naodim);
    for(int i=0; i<nao; i++) {
      for(int j=0; j<nmo; j++) {
        coeffs_nb.set_element(i,j,coeffs.get_element(i,j));
      }
    }

    Ref<OneBodyInt> Hcore=integral()->hcore();
    const double *buffer_Hcore = Hcore->buffer();
    Ref<GaussianBasisSet> basis = this->basis();
    int nshell = basis->nshell();

    RefSymmSCMatrix Hcore_ao = localkit->symmmatrix(naodim);

    for(int P=0; P<nshell; P++) {
      int nump = basis->shell(P).nfunction();
      for(int Q=0; Q<nshell; Q++) {
        int numq = basis->shell(Q).nfunction();

        Hcore->compute_shell(P,Q);
        int index = 0;
        for(int p=0; p<nump; p++) {
          int op = basis->shell_to_function(P) + p;
          for(int q=0; q<numq; q++) {
            int oq = basis->shell_to_function(Q) + q;
            index = p * numq + q;

            Hcore_ao.set_element(op,oq,buffer_Hcore[index]);
          }
        }
      }
    }

    RefSymmSCMatrix Hcore_mo = localkit->symmmatrix(nmodim);
    Hcore_mo.assign(0.0);
    Hcore_mo.accumulate_transform(coeffs_nb,Hcore_ao,SCMatrix::TransposeTransform);
    //Hcore_mo.accumulate_transform(coeffs_nb,Hcore_ao);

    if(debug_>=DefaultPrintThresholds::mostN2) {
      Hcore_mo->print("total hcore_mo");
    }

    return(Hcore_mo);
  }

  RefSymmSCMatrix PsiCorrWavefunction_PT2R12::hcore_mo(SpinCase1 spin) {
    const Ref<SCMatrixKit> localkit = new LocalSCMatrixKit;
    const Ref<OrbitalSpace> space = r12evalinfo_->refinfo()->orbs(spin);
    const int nmo = space->rank();
    const RefSCMatrix coeffs = space->coefs();
    const RefSCDimension nmodim = new SCDimension(nmo);
    const RefSCDimension naodim = coeffs.coldim();
    const int nao = naodim.n();
    RefSCMatrix coeffs_nb = localkit->matrix(naodim,nmodim);
    for(int i=0; i<nao; i++) {
      for(int j=0; j<nmo; j++) {
        coeffs_nb.set_element(i,j,coeffs.get_element(i,j));
      }
    }

    const Ref<OneBodyInt> Hcore=integral()->hcore();
    const double *buffer_Hcore = Hcore->buffer();
    const Ref<GaussianBasisSet> basis = this->basis();
    const int nshell = basis->nshell();

    RefSymmSCMatrix Hcore_ao = localkit->symmmatrix(naodim);

    for(int P=0; P<nshell; P++) {
      int nump = basis->shell(P).nfunction();
      for(int Q=0; Q<nshell; Q++) {
        int numq = basis->shell(Q).nfunction();

        Hcore->compute_shell(P,Q);
        int index = 0;
        for(int p=0; p<nump; p++) {
          int op = basis->shell_to_function(P) + p;
          for(int q=0; q<numq; q++) {
            int oq = basis->shell_to_function(Q) + q;
            index = p * numq + q;

            Hcore_ao.set_element(op,oq,buffer_Hcore[index]);
          }
        }
      }
    }

    RefSymmSCMatrix Hcore_mo = localkit->symmmatrix(nmodim);
    Hcore_mo.assign(0.0);
    Hcore_mo.accumulate_transform(coeffs_nb,Hcore_ao,SCMatrix::TransposeTransform);

    if(debug_>=DefaultPrintThresholds::mostN2) {
      Hcore_mo.print(prepend_spincase(spin,"hcore_mo").c_str());
    }

    return(Hcore_mo);
  }

  RefSymmSCMatrix PsiCorrWavefunction_PT2R12::moints() {
    Ref<SCMatrixKit> localkit = new LocalSCMatrixKit;
    Ref<OrbitalSpace> space = r12evalinfo_->refinfo()->orbs(Alpha);

    RefSCMatrix coeffs = space->coefs();

    RefSCDimension nmodim = coeffs.rowdim();
    RefSCDimension naodim = coeffs.coldim();
    int nmo = nmodim.n();
    int nao = naodim.n();
    RefSCDimension nmodim_nb = new SCDimension(nmo);
    RefSCDimension naodim_nb = new SCDimension(nao);

    RefSCMatrix coeffs_nb = localkit->matrix(naodim,nmodim);
    for(int i=0; i<nao; i++) {
      for(int j=0; j<nmo; j++) {
        coeffs_nb.set_element(i,j,coeffs.get_element(i,j));
      }
    }

    RefSCDimension triangdim = new SCDimension(triang_half_INDEX(nmo-1,nmo-1)+1);

    /// getting AO integrals
    RefSymmSCMatrix aoints = localkit->symmmatrix(triangdim);

    Ref<TwoBodyInt> twoint = integral()->electron_repulsion();
    const double *buffer = twoint->buffer();

    Ref<GaussianBasisSet> basis = this->basis();

    int nshell = basis->nshell();
    for(int P = 0; P < nshell; P++) {
      int nump = basis->shell(P).nfunction();

      for(int Q = 0; Q < nshell; Q++) {
        int numq = basis->shell(Q).nfunction();

        for(int R = 0; R < nshell; R++) {
          int numr = basis->shell(R).nfunction();

          for(int S = 0; S < nshell; S++) {
            int nums = basis->shell(S).nfunction();

            twoint->compute_shell(P,Q,R,S);

            int index = 0;
            for(int p=0; p < nump; p++) {
              int op = basis->shell_to_function(P)+p;

              for(int q = 0; q < numq; q++) {
                int oq = basis->shell_to_function(Q)+q;

                for(int r = 0; r < numr; r++) {
                  int oor = basis->shell_to_function(R)+r;

                  for(int s = 0; s < nums; s++,index++) {
                    int os = basis->shell_to_function(S)+s;

                    aoints.set_element(triang_half_INDEX(op,oq),triang_half_INDEX(oor,os),buffer[index]);

                  }
                }
              }
            }

          }
        }
      }
    }

    twoint = 0;

    /// tranformation of rs indices to MOs.
    RefSCMatrix moints_pqRS = localkit->matrix(triangdim,triangdim); // moints in chemist's notation order
    RefSymmSCMatrix mat_rs = localkit->symmmatrix(naodim);
    RefSymmSCMatrix mat_RS = localkit->symmmatrix(nmodim);
    for(int p=0; p<nmo; p++) {
      for(int q=0; q<=p; q++) {
        int ind_pq = triang_half_INDEX(p,q);
        RefSCVector vec_rs = aoints.get_row(ind_pq);
        vector_to_symmmatrix(mat_rs,vec_rs);
        mat_RS.assign(0.0);
        mat_RS.accumulate_transform(coeffs_nb,mat_rs,SCMatrix::TransposeTransform);
        symmmatrix_to_vector(vec_rs,mat_RS);
        moints_pqRS.assign_row(vec_rs,ind_pq);
      }
    }

    aoints = RefSymmSCMatrix();

    /// transformation of pq indices to MOs
    RefSymmSCMatrix moints_PQRS = localkit->symmmatrix(triangdim);  // moints in chemist's notation order
    RefSymmSCMatrix mat_pq = localkit->symmmatrix(naodim);
    RefSymmSCMatrix mat_PQ = localkit->symmmatrix(nmodim);
    for(int R=0; R<nmo; R++) {
      for(int S=0; S<=R; S++) {
        int ind_RS = triang_half_INDEX(R,S);
        RefSCVector vec_pq = moints_pqRS.get_column(ind_RS);
        vector_to_symmmatrix(mat_pq,vec_pq);
        mat_PQ.assign(0.0);
        mat_PQ.accumulate_transform(coeffs_nb,mat_pq,SCMatrix::TransposeTransform);
        symmmatrix_to_vector(vec_pq,mat_PQ);
        moints_PQRS.assign_row(vec_pq,ind_RS);
      }
    }

    moints_pqRS = RefSCMatrix();

    return(moints_PQRS);
  }

  RefSCMatrix PsiCorrWavefunction_PT2R12::moints(SpinCase2 pairspin) {
    Ref<SCMatrixKit> localkit = new LocalSCMatrixKit;
    Ref<OrbitalSpace> space1;
    Ref<OrbitalSpace> space2;

    if(pairspin == AlphaBeta) {
      space1 = r12evalinfo_->refinfo()->orbs(Alpha);
      space2 = r12evalinfo_->refinfo()->orbs(Beta);
    }
    if(pairspin == AlphaAlpha) {
      space1 = r12evalinfo_->refinfo()->orbs(Alpha);
      space2 = r12evalinfo_->refinfo()->orbs(Alpha);
    }
    if(pairspin == BetaBeta) {
      space1 = r12evalinfo_->refinfo()->orbs(Beta);
      space2 = r12evalinfo_->refinfo()->orbs(Beta);
    }

    RefSCMatrix coeffs1 = space1->coefs();
    RefSCMatrix coeffs2 = space2->coefs();

    RefSCDimension nmodim = coeffs1.rowdim();
    RefSCDimension naodim = coeffs1.coldim();
    int nmo = nmodim.n();
    int nao = naodim.n();
    RefSCDimension nmodim_nb = new SCDimension(nmo);
    RefSCDimension naodim_nb = new SCDimension(nao);

    RefSCMatrix coeffs1_nb = localkit->matrix(naodim,nmodim);
    for(int i=0; i<nao; i++) {
      for(int j=0; j<nmo; j++) {
        coeffs1_nb.set_element(i,j,coeffs1.get_element(i,j));
      }
    }
    RefSCMatrix coeffs2_nb = localkit->matrix(naodim,nmodim);
    for(int i=0; i<nao; i++) {
      for(int j=0; j<nmo; j++) {
        coeffs2_nb.set_element(i,j,coeffs2.get_element(i,j));
      }
    }

    RefSCDimension triangdim = new SCDimension(triang_half_INDEX(nmo-1,nmo-1)+1);

    /// getting AO integrals
    RefSymmSCMatrix aoints = localkit->symmmatrix(triangdim);

    Ref<TwoBodyInt> twoint = integral()->electron_repulsion();
    const double *buffer = twoint->buffer();

    Ref<GaussianBasisSet> basis = this->basis();

    int nshell = basis->nshell();
    for(int P = 0; P < nshell; P++) {
      int nump = basis->shell(P).nfunction();

      for(int Q = 0; Q < nshell; Q++) {
        int numq = basis->shell(Q).nfunction();

        for(int R = 0; R < nshell; R++) {
          int numr = basis->shell(R).nfunction();

          for(int S = 0; S < nshell; S++) {
            int nums = basis->shell(S).nfunction();

            twoint->compute_shell(P,Q,R,S);

            int index = 0;
            for(int p=0; p < nump; p++) {
              int op = basis->shell_to_function(P)+p;

              for(int q = 0; q < numq; q++) {
                int oq = basis->shell_to_function(Q)+q;

                for(int r = 0; r < numr; r++) {
                  int oor = basis->shell_to_function(R)+r;

                  for(int s = 0; s < nums; s++,index++) {
                    int os = basis->shell_to_function(S)+s;


                    aoints.set_element(triang_half_INDEX(op,oq),triang_half_INDEX(oor,os),buffer[index]);

                  }
                }
              }
            }

          }
        }
      }
    }

    twoint = 0;

    /// tranformation of rs indices to MOs.
    RefSCMatrix moints_pqRS = localkit->matrix(triangdim,triangdim); // moints in chemist's notation order
    RefSymmSCMatrix mat_rs = localkit->symmmatrix(naodim);
    RefSymmSCMatrix mat_RS = localkit->symmmatrix(nmodim);
    for(int p=0; p<nmo; p++) {
      for(int q=0; q<=p; q++) {
        int ind_pq = triang_half_INDEX(p,q);
        RefSCVector vec_rs = aoints.get_row(ind_pq);
        vector_to_symmmatrix(mat_rs,vec_rs);
        mat_RS.assign(0.0);
        mat_RS.accumulate_transform(coeffs2_nb,mat_rs,SCMatrix::TransposeTransform);
        symmmatrix_to_vector(vec_rs,mat_RS);
        moints_pqRS.assign_row(vec_rs,ind_pq);
      }
    }

    aoints = RefSymmSCMatrix();

    /// transformation of pq indices to MOs
    RefSCMatrix moints_PQRS = localkit->matrix(triangdim,triangdim);  // moints in chemist's notation order
    RefSymmSCMatrix mat_pq = localkit->symmmatrix(naodim);
    RefSymmSCMatrix mat_PQ = localkit->symmmatrix(nmodim);
    for(int R=0; R<nmo; R++) {
      for(int S=0; S<=R; S++) {
        int ind_RS = triang_half_INDEX(R,S);
        RefSCVector vec_pq = moints_pqRS.get_column(ind_RS);
        vector_to_symmmatrix(mat_pq,vec_pq);
        mat_PQ.assign(0.0);
        mat_PQ.accumulate_transform(coeffs1_nb,mat_pq,SCMatrix::TransposeTransform);
        symmmatrix_to_vector(vec_pq,mat_PQ);
        moints_PQRS.assign_row(vec_pq,ind_RS);
      }
    }

    moints_pqRS = RefSCMatrix();

    return(moints_PQRS);
  }

  RefSCMatrix PsiCorrWavefunction_PT2R12::g(SpinCase2 pairspin) {
    const Ref<SCMatrixKit> localkit = new LocalSCMatrixKit;
    const SpinCase1 spin1 = case1(pairspin);
    const SpinCase1 spin2 = case2(pairspin);
    const Ref<OrbitalSpace> space1 = r12evalinfo_->refinfo()->orbs(spin1);
    const Ref<OrbitalSpace> space2 = r12evalinfo_->refinfo()->orbs(spin2);
    Ref<SpinMOPairIter> PQ_iter = new SpinMOPairIter(space1,space2,pairspin);
    const int nmo1 = space1->rank();
    const int nmo2 = space2->rank();
    int nmo;
    if(nmo1==nmo2) {
      nmo = nmo1;
    }
    else {
      throw ProgrammingError("PsiCorrWavefunction_PT2R12::g(SpinCase2) -- nmo dimension of the spaces does not fit.",__FILE__,__LINE__);
    }
    const int pairrank = PQ_iter->nij();
    const RefSCDimension pairdim = new SCDimension(pairrank);
    RefSCMatrix G = localkit->matrix(pairdim,pairdim);
    G.assign(0.0);

    const bool antisymmetrize = (pairspin==AlphaBeta) ? false : true;
    std::vector<std::string> tforms;
    {
      const std::string tform_key = ParsedTwoBodyFourCenterIntKey::key(space1->id(),space2->id(),
                                                                       space1->id(),space2->id(),
                                                                       std::string("ERI"),
                                                                       std::string(TwoBodyIntLayout::b1b2_k1k2));
      tforms.push_back(tform_key);
    }

    r12eval_->compute_tbint_tensor<ManyBodyTensors::I_to_T,false,false>(G,TwoBodyOper::eri,space1,space1,space2,space2,
        antisymmetrize,tforms);

    return(G);
  }

  RefSCMatrix PsiCorrWavefunction_PT2R12::f(SpinCase1 spin) {
    const Ref<OrbitalSpace> bra_space = r12eval_->ggspace(spin);
    const Ref<OrbitalSpace> ket_space = r12eval_->ggspace(spin);
    RefSCMatrix fmat = r12eval_->fock(bra_space,ket_space,spin);

    return(fmat);
  }

  RefSymmSCMatrix PsiCorrWavefunction_PT2R12::onepdm_transformed(const SpinCase1 &spin) {
    const int nmo = reference()->nmo();
    const RefSymmSCMatrix opdm_pq = onepdm(spin);
    const RefSCMatrix transform = MPQC2PSI_transform_matrix(spin);

    RefSymmSCMatrix opdm_PQ = opdm_pq.clone();
    opdm_PQ.assign(0.0);
    opdm_PQ.accumulate_transform(transform,opdm_pq);

    return(opdm_PQ);
  }

  RefSymmSCMatrix PsiCorrWavefunction_PT2R12::onepdm_transformed() {
    const int nmo = reference()->nmo();
    const RefSymmSCMatrix opdm_pq = onepdm();
    const RefSCMatrix transform = MPQC2PSI_transform_matrix(Alpha);

    RefSymmSCMatrix opdm_PQ = opdm_pq.clone();
    opdm_PQ.assign(0.0);
    opdm_PQ.accumulate_transform(transform,opdm_pq);

    return(opdm_PQ);
  }

  RefSymmSCMatrix PsiCorrWavefunction_PT2R12::twopdm_transformed() {
    const Ref<SCMatrixKit> localkit = new LocalSCMatrixKit;
    const int nmo = reference()->nmo();
    const RefSCDimension nmodim = new SCDimension(nmo);
    RefSymmSCMatrix tpdm_pqrs = twopdm();
    //tpdm_pqrs.print("tpdm_pqrs");
    const RefSCMatrix transform = MPQC2PSI_transform_matrix(Alpha);
    const RefSCDimension halfdim = tpdm_pqrs.dim();

    RefSCMatrix tpdm_pqRS = localkit->matrix(halfdim,halfdim);
    RefSCMatrix mat_rs = localkit->matrix(nmodim,nmodim);
    RefSCMatrix mat_RS = localkit->matrix(nmodim,nmodim);
    for(int p=0; p<nmo; p++) {
      for(int q=0; q<nmo; q++) {
        int ind_pq = ordinary_INDEX(p,q,nmo);
        RefSCVector vec_rs = tpdm_pqrs.get_row(ind_pq);
        vector_to_matrix(mat_rs,vec_rs);
        //mat_RS.assign(0.0);
        //mat_RS.accumulate_transform(transform,mat_rs);
        mat_RS = transform*mat_rs*transform.t();
        matrix_to_vector(vec_rs,mat_RS);
        tpdm_pqRS.assign_row(vec_rs,ind_pq);
      }
    }

    tpdm_pqrs = 0;   // deallocate tpdm_pqrs

    RefSymmSCMatrix tpdm_PQRS = localkit->symmmatrix(halfdim);
    RefSCMatrix mat_pq = localkit->matrix(nmodim,nmodim);
    RefSCMatrix mat_PQ = localkit->matrix(nmodim,nmodim);
    for(int R=0; R<nmo; R++) {
      for(int S=0; S<nmo; S++) {
        int ind_RS = ordinary_INDEX(R,S,nmo);
        RefSCVector vec_pq = tpdm_pqRS.get_column(ind_RS);
        vector_to_matrix(mat_pq,vec_pq);
        //mat_PQ.assign(0.0);
        //mat_PQ.accumulate_transform(transform,mat_pq);
        mat_PQ  = transform*mat_pq*transform.t();
        matrix_to_vector(vec_pq,mat_PQ);
        tpdm_PQRS.assign_row(vec_pq,ind_RS);
      }
    }

    tpdm_pqRS = 0;   // deallocate tpdm_pqRS

    //tpdm_PQRS.print("tpdm_PQRS");

    return(tpdm_PQRS);
  }

  RefSymmSCMatrix PsiCorrWavefunction_PT2R12::twopdm_transformed_dirac() {
    const Ref<SCMatrixKit> localkit = new LocalSCMatrixKit;
    const int nmo = reference()->nmo();
    const RefSCDimension nmodim = new SCDimension(nmo);
    RefSymmSCMatrix tpdm_pqrs = twopdm_dirac();

    const RefSCMatrix transform = MPQC2PSI_transform_matrix(Alpha);
    const RefSCDimension halfdim = tpdm_pqrs.dim();

    RefSCMatrix tpdm_pqRS = localkit->matrix(halfdim,halfdim);
    RefSCMatrix mat_rs = localkit->matrix(nmodim,nmodim);
    RefSCMatrix mat_RS = localkit->matrix(nmodim,nmodim);
    for(int p=0; p<nmo; p++) {
      for(int q=0; q<nmo; q++) {
        int ind_pq = ordinary_INDEX(p,q,nmo);
        RefSCVector vec_rs = tpdm_pqrs.get_row(ind_pq);
        vector_to_matrix(mat_rs,vec_rs);
        mat_RS = transform*mat_rs*transform.t();
        matrix_to_vector(vec_rs,mat_RS);
        tpdm_pqRS.assign_row(vec_rs,ind_pq);
      }
    }

    tpdm_pqrs = 0;   // deallocate tpdm_pqrs

    RefSymmSCMatrix tpdm_PQRS = localkit->symmmatrix(halfdim);
    RefSCMatrix mat_pq = localkit->matrix(nmodim,nmodim);
    RefSCMatrix mat_PQ = localkit->matrix(nmodim,nmodim);
    for(int R=0; R<nmo; R++) {
      for(int S=0; S<nmo; S++) {
        int ind_RS = ordinary_INDEX(R,S,nmo);
        RefSCVector vec_pq = tpdm_pqRS.get_column(ind_RS);
        vector_to_matrix(mat_pq,vec_pq);
        //mat_PQ.assign(0.0);
        //mat_PQ.accumulate_transform(transform,mat_pq);
        mat_PQ  = transform*mat_pq*transform.t();
        matrix_to_vector(vec_pq,mat_PQ);
        tpdm_PQRS.assign_row(vec_pq,ind_RS);
      }
    }

    tpdm_pqRS = 0;   // deallocate tpdm_pqRS

    return(tpdm_PQRS);
  }

  RefSymmSCMatrix PsiCorrWavefunction_PT2R12::twopdm_transformed_dirac(const SpinCase2 &pairspin) {
    const Ref<SCMatrixKit> localkit = new LocalSCMatrixKit;
    const int nmo = reference()->nmo();
    const RefSCDimension nmodim = new SCDimension(nmo);
    RefSymmSCMatrix tpdm_pqrs = twopdm_dirac(pairspin);

#if 0
    if(pairspin==AlphaBeta) {
      double tol=1.0e-6;
      RefSymmSCMatrix tpdm_vgl = twopdm_dirac();
      RefSymmSCMatrix tpdm_vgl1 = twopdm_dirac_from_components();
      //tpdm_pqrs->print("tpdm_pqrs");
      //tpdm_vgl->print("tpdm_vgl");
      ExEnv::out0() << "Begin testing PsiCorrWavefunction_PT2R12::twopdm_transformed_dirac(const SpinCase2 &) ..." << endl;
      for(int p=0; p<nmo; p++) {
        for(int q=0; q<nmo; q++) {
          for(int r=0; r<nmo; r++) {
            for(int s=0; s<nmo; s++) {
#if 1
              if(fabs(tpdm_vgl->get_element(ordinary_INDEX(p,q,nmo),ordinary_INDEX(r,s,nmo))-tpdm_vgl1->get_element(ordinary_INDEX(p,q,nmo),ordinary_INDEX(r,s,nmo)))>tol) {
                ExEnv::out0() << p << "  " << q << "  " << r << "  " << s  << "  "
                              << setprecision(12) << tpdm_vgl->get_element(ordinary_INDEX(p,q,nmo),ordinary_INDEX(r,s,nmo)) << "  "
                              << setprecision(12) << tpdm_vgl1->get_element(ordinary_INDEX(p,q,nmo),ordinary_INDEX(r,s,nmo)) << endl;
              }
#else
              if(fabs(tpdm_vgl->get_element(ordinary_INDEX(p,q,nmo),ordinary_INDEX(r,s,nmo))-2.0*tpdm_pqrs->get_element(ordinary_INDEX(p,q,nmo),ordinary_INDEX(r,s,nmo)))>tol) {
                ExEnv::out0() << p << "  " << q << "  " << r << "  " << s  << "  "
                              << setprecision(12) << tpdm_vgl->get_element(ordinary_INDEX(p,q,nmo),ordinary_INDEX(r,s,nmo)) << "  "
                              << setprecision(12) << 2.0*tpdm_pqrs->get_element(ordinary_INDEX(p,q,nmo),ordinary_INDEX(r,s,nmo)) << endl;
              }
#endif
            }
          }
        }
      }
      ExEnv::out0() << "End testing PsiCorrWavefunction_PT2R12::twopdm_transformed_dirac(const SpinCase2 &) ." << endl;
    }
#endif

    const SpinCase1 spin1 = case1(pairspin);
    const SpinCase1 spin2 = case2(pairspin);
    const RefSCMatrix transform_spin1 = MPQC2PSI_transform_matrix(spin1);
    const RefSCMatrix transform_spin2 = MPQC2PSI_transform_matrix(spin2);
    const RefSCDimension halfdim = tpdm_pqrs.dim();

    RefSCMatrix tpdm_pqRS = localkit->matrix(halfdim,halfdim);
    tpdm_pqRS->assign(0.0);
    RefSCMatrix mat_rs = localkit->matrix(nmodim,nmodim);
    RefSCMatrix mat_RS = localkit->matrix(nmodim,nmodim);
    for(int p=0; p<nmo; p++) {
      const int qmax = (pairspin==AlphaBeta) ? nmo : p;
      for(int q=0; q<qmax; q++) {
        const int ind_pq = (pairspin==AlphaBeta) ? ordinary_INDEX(p,q,nmo) : lowerupper_index(p,q);
        RefSCVector vec_rs = tpdm_pqrs.get_row(ind_pq);
        vector_to_matrix(mat_rs,vec_rs,pairspin);
        mat_RS = transform_spin1*mat_rs*transform_spin2.t();
        matrix_to_vector(vec_rs,mat_RS,pairspin);
        tpdm_pqRS.assign_row(vec_rs,ind_pq);
      }
    }

    tpdm_pqrs = 0;   // deallocate tpdm_pqrs

    RefSymmSCMatrix tpdm_PQRS = localkit->symmmatrix(halfdim);
    tpdm_PQRS->assign(0.0);
    RefSCMatrix mat_pq = localkit->matrix(nmodim,nmodim);
    RefSCMatrix mat_PQ = localkit->matrix(nmodim,nmodim);
    for(int R=0; R<nmo; R++) {
      const int Smax = (pairspin==AlphaBeta) ? nmo : R;
      for(int S=0; S<Smax; S++) {
        const int ind_RS = (pairspin==AlphaBeta) ? ordinary_INDEX(R,S,nmo) : lowerupper_index(R,S);
        RefSCVector vec_pq = tpdm_pqRS.get_column(ind_RS);
        vector_to_matrix(mat_pq,vec_pq,pairspin);
        mat_PQ  = transform_spin1*mat_pq*transform_spin2.t();
        matrix_to_vector(vec_pq,mat_PQ,pairspin);
        tpdm_PQRS.assign_row(vec_pq,ind_RS);  // assign_column does not exist. assign_row can be used because it is a symmetric matrix.
      }
    }

    tpdm_pqRS = 0;   // deallocate tpdm_pqRS

    return(tpdm_PQRS);
  }

  RefSymmSCMatrix PsiCorrWavefunction_PT2R12::onepdm_refmo_HF(const SpinCase1 &spin) {
    int nmo = reference()->nmo();
    int nocc = reference()->nocc(spin);
    Ref<SCMatrixKit> localkit = new LocalSCMatrixKit;
    RefSCDimension opdm_dim(new SCDimension(nmo));
    RefSymmSCMatrix opdm_hf(localkit->symmmatrix(opdm_dim));
    opdm_hf->assign(0.0);

    for(int i=0; i<nocc; i++) {
      opdm_hf->set_element(i,i,1.0);
    }

    return(opdm_hf);
  }

  RefSymmSCMatrix PsiCorrWavefunction_PT2R12::onepdm_refmo(const SpinCase1 &spin) {
    RefSymmSCMatrix opdm_so;
    if(opdm_type_==correlated) {
      //opdm_so = onepdm_refmo_closedshell();
      //ExEnv::out0() << "Calling onepdm_transformed." << endl;
      opdm_so = onepdm_transformed(spin);
    }
    else if(opdm_type_==HF) {
      opdm_so = onepdm_refmo_HF(spin);
    }
    else {
      throw InputError("PsiCorrWavefunction_PT2R12::onepdm_refmo -- opdm_type not known.",__FILE__,__LINE__);
    }

    return(opdm_so);
  }

  RefSymmSCMatrix PsiCorrWavefunction_PT2R12::twopdm_refmo_from_onepdm(const SpinCase2 &pairspin) {
    Ref<SCMatrixKit> localkit = new LocalSCMatrixKit;

    int nall_a = r12evalinfo_->refinfo()->orbs(Alpha)->rank();
    int nall_b = r12evalinfo_->refinfo()->orbs(Beta)->rank();
    RefSCDimension dim_aa;
    if(pairspin==AlphaBeta) {
      dim_aa = new SCDimension(nall_a*nall_b);
    }
    else if(pairspin==AlphaAlpha) {
      dim_aa = new SCDimension((nall_a*(nall_a-1))/2);
    }
    else if(pairspin==BetaBeta) {
      dim_aa = new SCDimension((nall_b*(nall_b-1))/2);
    }

    //RefSymmSCMatrix tpdm_so = localkit->symmmatrix(r12eval_->dim_aa(pairspin));
    RefSymmSCMatrix tpdm_so = localkit->symmmatrix(dim_aa);

    const SpinCase1 spin1 = case1(pairspin);
    const SpinCase1 spin2 = case2(pairspin);
    const RefSymmSCMatrix opdm1 = onepdm_refmo(spin1);
    const RefSymmSCMatrix opdm2  = onepdm_refmo(spin2);
    const Ref<OrbitalSpace> orbs1 = r12evalinfo_->refinfo()->orbs(spin1);
    const Ref<OrbitalSpace> orbs2 = r12evalinfo_->refinfo()->orbs(spin2);
    SpinMOPairIter PQ_iter(orbs1,orbs2,pairspin);
    SpinMOPairIter RS_iter(orbs1,orbs2,pairspin);
    for(PQ_iter.start(); int(PQ_iter); PQ_iter.next()) {
      const int PQ = PQ_iter.ij();
      const int P = PQ_iter.i();
      const int Q = PQ_iter.j();
      for(RS_iter.start(); int(RS_iter); RS_iter.next()) {
        const int RS = RS_iter.ij();
        const int R = RS_iter.i();
        const int S = RS_iter.j();
        if (pairspin == AlphaBeta)
          tpdm_so.set_element(PQ,RS,opdm1.get_element(P,R)*opdm2.get_element(Q,S));
        else
          tpdm_so.set_element(PQ,RS,opdm1.get_element(P,R)*opdm2.get_element(Q,S)
                              -opdm1.get_element(Q,R)*opdm2.get_element(P,S));
      }
    }
    return(tpdm_so);
  }

  RefSymmSCMatrix PsiCorrWavefunction_PT2R12::twopdm_refmo(const SpinCase2 &pairspin) {
    RefSymmSCMatrix tpdm_so;

    if(tpdm_from_opdms_) {
      tpdm_so = twopdm_refmo_from_onepdm(pairspin);
    }
    else {
      tpdm_so = twopdm_transformed_dirac(pairspin);
    }

    return(tpdm_so);
  }

  RefSymmSCMatrix PsiCorrWavefunction_PT2R12::lambda_refmo(const SpinCase2 &pairspin) {
    SpinCase1 spin1 = case1(pairspin);
    SpinCase1 spin2 = case2(pairspin);

    RefSymmSCMatrix opdm1 = onepdm_refmo(spin1);
    RefSymmSCMatrix opdm2 = onepdm_refmo(spin2);
    RefSymmSCMatrix tpdm = twopdm_refmo(pairspin);
    RefSymmSCMatrix lambda = tpdm.clone();
    lambda.assign(0.0);
    const int nmo = opdm1.dim().n();

    SpinMOPairIter pq_iter(r12evalinfo_->refinfo()->orbs(spin1),r12evalinfo_->refinfo()->orbs(spin2),pairspin);
    SpinMOPairIter rs_iter(r12evalinfo_->refinfo()->orbs(spin1),r12evalinfo_->refinfo()->orbs(spin2),pairspin);
    for(pq_iter.start(); int(pq_iter); pq_iter.next()) {
      const int p = pq_iter.i();
      const int q = pq_iter.j();
      const int pq = pq_iter.ij();
      for(rs_iter.start(); int(rs_iter); rs_iter.next()) {
        const int r = rs_iter.i();
        const int s = rs_iter.j();
        const int rs = rs_iter.ij();
        lambda.set_element(pq,rs,tpdm.get_element(pq,rs)-opdm1.get_element(p,r)*opdm2.get_element(q,s));
        if(pairspin!=AlphaBeta) {
          lambda.accumulate_element(pq,rs,opdm1.get_element(q,r)*opdm2.get_element(p,s));
        }
      }
    }

    return(lambda);
  }

  RefSCMatrix PsiCorrWavefunction_PT2R12::C(SpinCase2 S) {
    //return(r12energy_->C(S));
    Ref<LocalSCMatrixKit> local_matrix_kit = new LocalSCMatrixKit();
    RefSCMatrix Cmat = local_matrix_kit->matrix(r12eval_->dim_GG(S),r12eval_->dim_gg(S));
    if(S==AlphaBeta) {
      SpinMOPairIter OW_iter(r12eval_->GGspace(Alpha), r12eval_->GGspace(Beta), S );
      SpinMOPairIter PQ_iter(r12eval_->ggspace(Alpha), r12eval_->GGspace(Beta), S );
      Ref<LinearR12::GeminalDescriptor> geminaldesc = r12evalinfo_->corrfactor()->geminaldescriptor();
      CuspConsistentGeminalCoefficient coeff_gen(S,geminaldesc);
      for(OW_iter.start(); int(OW_iter); OW_iter.next()) {
        for(PQ_iter.start(); int(PQ_iter); PQ_iter.next()) {
          unsigned int O = OW_iter.i();
          unsigned int W = OW_iter.j();
          unsigned int P = PQ_iter.i();
          unsigned int Q = PQ_iter.j();
          int OW = OW_iter.ij();
          int PQ = PQ_iter.ij();
          Cmat.set_element(OW,PQ,coeff_gen.C(O,W,P,Q));
        }
      }
    }
    else {
      SpinCase1 spin = (S==AlphaAlpha) ? Alpha : Beta;
      SpinMOPairIter OW_iter(r12eval_->GGspace(spin), r12eval_->GGspace(spin), S );
      SpinMOPairIter PQ_iter(r12eval_->ggspace(spin), r12eval_->GGspace(spin), S );
      Ref<LinearR12::GeminalDescriptor> geminaldesc = r12evalinfo_->corrfactor()->geminaldescriptor();
      CuspConsistentGeminalCoefficient coeff_gen(S,geminaldesc);
      for(OW_iter.start(); int(OW_iter); OW_iter.next()) {
        for(PQ_iter.start(); int(PQ_iter); PQ_iter.next()) {
          unsigned int O = OW_iter.i();
          unsigned int W = OW_iter.j();
          unsigned int P = PQ_iter.i();
          unsigned int Q = PQ_iter.j();
          int OW = OW_iter.ij();
          int PQ = PQ_iter.ij();
          Cmat.set_element(OW,PQ,coeff_gen.C(O,W,P,Q));
        }
      }
    }
    return(Cmat);
  }

  RefSCMatrix PsiCorrWavefunction_PT2R12::V_genref_projector2(SpinCase2 pairspin) {
    SpinCase1 spin1 = case1(pairspin);
    SpinCase1 spin2 = case2(pairspin);

    const Ref<OrbitalSpace>& GG1_space = r12eval_->GGspace(spin1);
    const Ref<OrbitalSpace>& GG2_space = r12eval_->GGspace(spin2);
    const Ref<OrbitalSpace>& orbs1 = r12evalinfo_->refinfo()->orbs(spin1);
    const Ref<OrbitalSpace>& orbs2 = r12evalinfo_->refinfo()->orbs(spin2);
    const Ref<OrbitalSpace>& cabs1 = r12evalinfo_->ribs_space(spin1);
    const Ref<OrbitalSpace>& cabs2 = r12evalinfo_->ribs_space(spin2);

    const Ref<OrbitalSpace>& f_p_A1 = r12eval_->F_p_A(spin1);
    const Ref<OrbitalSpace>& f_p_A2 = r12eval_->F_p_A(spin2);

    Ref<SCMatrixKit> localkit = new LocalSCMatrixKit;
    RefSCMatrix V_genref = localkit->matrix(r12eval_->dim_GG(pairspin),r12eval_->dim_gg(pairspin));
    V_genref.assign(0.0);

    const bool antisymmetrize = pairspin!=AlphaBeta;
    const bool part1_equiv_part2 = (GG1_space == GG2_space && orbs1 == orbs2);

    /** If COMPUTE_R_F_LAMBDA_TERM is defined, the term
     * "\bar{r}^{r s}_{p_2 \alpha} f^{\alpha}_{p_3} \lambda^{p_2 p_3}_{p q}"
     * is added to "V_{pq}^{rs}". */
//#define COMPUTE_R_F_LAMBDA_TERM
#ifdef COMPUTE_R_F_LAMBDA_TERM
    //
    // ij|A|kl = ij|f12|kl_f, symmetrized if part1_equiv_part2
    //
    RefSCMatrix A = localkit->matrix(r12eval_->dim_GG(pairspin),r12eval_->dim_gg(pairspin));
    A.assign(0.0);

    {
      std::vector<  Ref<TwoBodyMOIntsTransform> > tforms;
      {
        NewTransformCreator tform_creator(r12eval_,GG1_space,orbs1,GG2_space,f_p_A2,true);
        fill_container(tform_creator,tforms);
      }
      std::vector< Ref<TwoBodyIntDescr> > intdescr(1);
      intdescr[0]=tforms[0]->intdescr();

      r12eval_->compute_tbint_tensor<ManyBodyTensors::I_to_T,true,false>(
                         A,r12eval_->corrfactor()->tbint_type_f12(),
                         GG1_space,orbs1,GG2_space,f_p_A2,
                         antisymmetrize,tforms,intdescr);
    }
    if (part1_equiv_part2) {
      // no need to symmetrize if computing antisymmetric matrix -- compute_tbint_tensor takes care of that
      if (!antisymmetrize)
        symmetrize<false>(A,A,GG1_space,orbs1);
      A.scale(2.0);
    }
    else {
      std::vector<  Ref<TwoBodyMOIntsTransform> > tforms;
      {
        NewTransformCreator tform_creator(r12eval_,GG1_space,f_p_A1,GG2_space,orbs2,true);
        fill_container(tform_creator,tforms);
      }
      std::vector< Ref<TwoBodyIntDescr> > intdescr(1);
      intdescr[0]=tforms[0]->intdescr();
      r12eval_->compute_tbint_tensor<ManyBodyTensors::I_to_T,true,false>(
                         A,r12eval_->corrfactor()->tbint_type_f12(),
                         GG1_space,f_p_A1,GG2_space,orbs2,
                         antisymmetrize,tforms,intdescr);
    }

    const RefSymmSCMatrix lambda = lambda_refmo(pairspin);
    V_genref.accumulate(A * lambda);
#endif /* COMPUTE_R_F_LAMBDA_TERM */

    const RefSCMatrix V_intermed = r12eval_->V(pairspin);
    const RefSymmSCMatrix tpdm = twopdm_refmo(pairspin);
    V_genref.accumulate(V_intermed * tpdm);

    return(V_genref);
  }

  RefSCMatrix PsiCorrWavefunction_PT2R12::V_transformed_by_C(SpinCase2 pairspin) {
    RefSCMatrix T = C(pairspin);
    RefSCMatrix V = r12eval_->V(pairspin);
    RefSCMatrix Vtransposed = V.t();
    RefSCMatrix VT = Vtransposed*T;
    return(VT);
  }

  RefSymmSCMatrix PsiCorrWavefunction_PT2R12::X_transformed_by_C(SpinCase2 pairspin) {
    RefSCMatrix T = C(pairspin);
    RefSymmSCMatrix X = r12eval_->X(pairspin);
    RefSCDimension transformed_dim = T.coldim();
    RefSymmSCMatrix XT = T.kit()->symmmatrix(transformed_dim);
    XT.assign(0.0);
    XT.accumulate_transform(T,X,SCMatrix::TransposeTransform);

    return(XT);
  }

  RefSymmSCMatrix PsiCorrWavefunction_PT2R12::B_transformed_by_C(SpinCase2 pairspin) {
    RefSymmSCMatrix B = r12eval_->B(pairspin);
    RefSCMatrix T = C(pairspin);
    RefSCDimension transformed_dim = T.coldim();
    RefSymmSCMatrix BT = T.kit()->symmmatrix(transformed_dim);
    BT.assign(0.0);
    BT.accumulate_transform(T,B,SCMatrix::TransposeTransform);

    return(BT);
  }

  RefSCVector PsiCorrWavefunction_PT2R12::dipolemoments(const Ref<DipoleData> &dipoledata) {
    Ref<SCMatrixKit> localkit = new LocalSCMatrixKit;
    Ref<GaussianBasisSet> basis = this->basis();
    Ref<OrbitalSpace> space = r12evalinfo_->refinfo()->orbs(Alpha);
    const int nmo = reference_->nmo();

    // computing dipole integrals in MO basis
    RefSCMatrix MX, MY, MZ;
    RefSCMatrix MXX, MYY, MZZ;
    r12evalinfo_->compute_multipole_ints(space,space,MX,MY,MZ,MXX,MYY,MZZ);

    // computing oneparticle density matrix in MO basis
    RefSymmSCMatrix opdm = onepdm_transformed();
    //opdm.print("oneparticle density matrix in MO basis");

    RefSCDimension M_dim = opdm.dim();
    RefSCMatrix MX_nb = localkit->matrix(M_dim,M_dim);
    RefSCMatrix MY_nb = localkit->matrix(M_dim,M_dim);
    RefSCMatrix MZ_nb = localkit->matrix(M_dim,M_dim);
    for(int i=0; i<M_dim.n(); i++) {
      for(int j=0; j<M_dim.n(); j++) {
        MX_nb.set_element(i,j,MX.get_element(i,j));
        MY_nb.set_element(i,j,MY.get_element(i,j));
        MZ_nb.set_element(i,j,MZ.get_element(i,j));
      }
    }

    RefSCMatrix opdm_x = MX_nb * opdm;
    RefSCMatrix opdm_y = MY_nb * opdm;
    RefSCMatrix opdm_z = MZ_nb * opdm;

    RefSCVector dipoles = localkit->vector(RefSCDimension(new SCDimension(3)));
    dipoles.assign(0.0);

    double dx=0.0, dy=0.0, dz=0.0;
    for(int p=0; p<nmo; p++) {
      dx += opdm_x.get_element(p,p);
      dy += opdm_y.get_element(p,p);
      dz += opdm_z.get_element(p,p);
    }

    int natom = reference_->molecule()->natom();
    double dx_nuc=0.0, dy_nuc=0.0, dz_nuc=0.0;
    for(int i=0; i<natom; i++) {
      dx_nuc += reference_->molecule()->r(i,0) * reference_->molecule()->Z(i);
      dy_nuc += reference_->molecule()->r(i,1) * reference_->molecule()->Z(i);
      dz_nuc += reference_->molecule()->r(i,2) * reference_->molecule()->Z(i);
    }

    dipoles[0] = dx + dx_nuc;
    dipoles[1] = dy + dy_nuc;
    dipoles[2] = dz + dz_nuc;

    return(dipoles);
  }

  double PsiCorrWavefunction_PT2R12::energy_HF() {
    const int nelec = nelectron();
    const int nocc_trunc = nelec/2;
    const double dnocc = (double)nelec/2.0;
    int nocc;
    if((dnocc-(double)nelec) > 0.4) {
      nocc = nocc_trunc + 1;
    }
    else {
      nocc = nocc_trunc;
    }
    //ExEnv::out0() << "nocc = " << nocc << endl;
    const int nmo = reference_->nmo();

    const RefSymmSCMatrix hcore = hcore_mo();
    const RefSymmSCMatrix intsmo = moints();

    double energy = 0.0;
    double energy_hcore = 0.0;
    for(int i=0; i<nocc; i++) {
      energy_hcore += 2.0 * hcore.get_element(i,i);
    }
    energy += energy_hcore;
    ExEnv::out0() << "hcore contribution to HF energy: " << setprecision(12) << energy_hcore << endl;

    double energy_twoel = 0.0;
    for(int i=0; i<nocc; i++) {
      for(int j=0; j<nocc; j++) {
        energy_twoel += 2.0 * intsmo.get_element(triang_half_INDEX(i,i),triang_half_INDEX(j,j))
                       -intsmo.get_element(triang_half_INDEX(i,j),triang_half_INDEX(j,i));
      }
    }
    energy += energy_twoel;
    ExEnv::out0() << "two-electron contribution to HF energy: " << setprecision(12) << energy_twoel << endl;

    energy += reference_->nuclear_repulsion_energy();

    return(energy);
  }

  double PsiCorrWavefunction_PT2R12::energy_conventional() {
    double twoparticle_energy[NSpinCases2];
    double oneparticle_energy[NSpinCases1];
    const int npurespincases2 = (r12eval_->spin_polarized()) ? 3 : 2;
    const int npurespincases1 = (r12eval_->spin_polarized()) ? 2 : 1;
    const int nmo = r12evalinfo_->refinfo()->orbs()->rank();
    const RefSCDimension nmodim = new SCDimension(nmo);
    const Ref<SCMatrixKit> localkit = new LocalSCMatrixKit;

    for(int spincase1=0; spincase1<npurespincases1; spincase1++) {
      const SpinCase1 spin = static_cast<SpinCase1>(spincase1);
      const RefSymmSCMatrix opdm = onepdm_refmo(spin);
      const RefSymmSCMatrix hcore = hcore_mo(spin);
      oneparticle_energy[spin] = (hcore*opdm)->trace();
    }

    for(int spincase2=0; spincase2<npurespincases2; spincase2++) {
      const SpinCase2 pairspin = static_cast<SpinCase2>(spincase2);
      const SpinCase1 spin1 = case1(pairspin);
      const SpinCase1 spin2 = case2(pairspin);

      const RefSymmSCMatrix tpdm = twopdm_refmo(pairspin);
#if 0
      if(pairspin==AlphaBeta) {
        ExEnv::out0() << "start test in PsiCorrWavefunction_PT2R12::energy_conventional() ..." << endl;
        const double tol = 1.0e-6;
        const RefSymmSCMatrix tpdm_vgl = twopdm_transformed_dirac();
        for(int p=0; p<nmo; p++) {
          for(int q=0; q<nmo; q++) {
            const int ind_pq = ordinary_INDEX(p,q,nmo);
            for(int r=0; r<nmo; r++) {
              for(int s=0; s<nmo; s++) {
                const int ind_rs = ordinary_INDEX(r,s,nmo);
                const double val1 = 2.0 * tpdm->get_element(ind_pq,ind_rs);
                const double val2 = tpdm_vgl->get_element(ind_pq,ind_rs);
                if(fabs(val1-val2)>tol) {
                  ExEnv::out0() << p << "  " << q << "  " << r << "  " << s
                                << "  " << setprecision(12) << val1
                                << "  " << setprecision(12) << val2 << endl;
                }
              }
            }
          }
        }
        ExEnv::out0() << "Test in PsiCorrWavefunction_PT2R12::energy_conventional() finished." << endl;
      }
#endif
      const RefSCMatrix G = g(pairspin);

      twoparticle_energy[pairspin] = (G*tpdm).trace();
    }

    if(!r12eval_->spin_polarized()) {
      twoparticle_energy[BetaBeta] = twoparticle_energy[AlphaAlpha];
      oneparticle_energy[Beta] = oneparticle_energy[Alpha];
    }

//#define ENERGY_CONVENTIONAL_PRINT_CONTRIBUTIONS

    double energy_hcore = 0.0;
    for(int i=0; i<NSpinCases1; i++) {
      const SpinCase1 spin = static_cast<SpinCase1>(i);
#ifdef ENERGY_CONVENTIONAL_PRINT_CONTRIBUTIONS
      ExEnv::out0() << ((spin==Alpha) ? "Alpha" : "Beta") << " hcore energy: "
                    << setprecision(12) << oneparticle_energy[i] << endl;
#endif
      energy_hcore += oneparticle_energy[i];
    }
#ifdef ENERGY_CONVENTIONAL_PRINT_CONTRIBUTIONS
    ExEnv::out0() << "hcore energy: " << setprecision(12) << energy_hcore << endl;
#endif
    double energy_twoelec = 0.0;
    for(int i=0; i<NSpinCases2; i++) {
#ifdef ENERGY_CONVENTIONAL_PRINT_CONTRIBUTIONS
      string pairspin_str = "";
      if(i==0) {
        pairspin_str = "AlphaBeta";
      }
      else if(i==1) {
        pairspin_str = "AlphaAlpha";
      }
      else if(i==2) {
        pairspin_str = "BetaBeta";
      }
      ExEnv::out0() << "pairspin " << pairspin_str << " energy: " << setprecision(12) << twoparticle_energy[i] << endl;
#endif
      energy_twoelec += twoparticle_energy[i];
    }
#ifdef ENERGY_CONVENTIONAL_PRINT_CONTRIBUTIONS
    ExEnv::out0() << "two-electron energy: " << setprecision(12) << energy_twoelec << endl;
#endif
    double energy = energy_hcore + energy_twoelec;
    energy += reference_->nuclear_repulsion_energy();

    return(energy);
  }

  RefSCMatrix PsiCorrWavefunction_PT2R12::phi_HF(SpinCase2 pairspin) {
    RefSCMatrix fmat[NSpinCases1];
    for(int i=0; i<NSpinCases1; i++) {
      SpinCase1 spin = static_cast<SpinCase1>(i);
      fmat[i] = f(spin);

      //fmat[i].print(prepend_spincase(spin,"Fock matrix").c_str());
    }
    int nmo = reference()->nmo();
    RefSCMatrix phi;
    phi = C(pairspin).kit()->matrix(r12eval_->dim_gg(pairspin),r12eval_->dim_gg(pairspin));
    phi.assign(0.0);

    if(pairspin==AlphaBeta) {
      SpinMOPairIter UV_iter(r12evalinfo_->refinfo()->orbs(Alpha),r12evalinfo_->refinfo()->orbs(Beta),pairspin);
      SpinMOPairIter PQ_iter(r12evalinfo_->refinfo()->orbs(Alpha),r12evalinfo_->refinfo()->orbs(Beta),pairspin);
      for(PQ_iter.start(); int(PQ_iter); PQ_iter.next()) {
        const int P = PQ_iter.i();
        const int Q = PQ_iter.j();
        const int PQ = PQ_iter.ij();
        for(UV_iter.start(); int(UV_iter); UV_iter.next()) {
          const int U = UV_iter.i();
          const int V = UV_iter.j();
          const int UV = UV_iter.ij();
          if(U==P) {
            phi->accumulate_element(PQ,UV,-fmat[Beta].get_element(Q,V));
          }
          if(Q==V) {
            phi->accumulate_element(PQ,UV,-fmat[Alpha].get_element(P,U));
          }
        }
      }
    }  // pairspin==AlphaBeta
    else {  // pairspin!=AlphaBeta
      SpinCase1 spin = (pairspin == AlphaAlpha) ? Alpha : Beta;
      SpinCase1 otherspin = (spin==Alpha) ? Beta : Alpha;
      SpinMOPairIter UV_iter(r12evalinfo_->refinfo()->orbs(spin),r12evalinfo_->refinfo()->orbs(spin),pairspin);
      SpinMOPairIter PQ_iter(r12evalinfo_->refinfo()->orbs(spin),r12evalinfo_->refinfo()->orbs(spin),pairspin);

      for(PQ_iter.start(); int(PQ_iter); PQ_iter.next()) {
        const int P = PQ_iter.i();
        const int Q = PQ_iter.j();
        const int PQ = PQ_iter.ij();
        for(UV_iter.start(); int(UV_iter); UV_iter.next()) {
          const int U = UV_iter.i();
          const int V = UV_iter.j();
          const int UV = UV_iter.ij();

          if(V==Q) {
            phi->accumulate_element(PQ,UV,-1.0*fmat[spin].get_element(P,U));
          }
          if(U==Q) {
            phi->accumulate_element(PQ,UV,fmat[spin].get_element(P,V));
          }
          if(V==P) {
            phi->accumulate_element(PQ,UV,fmat[spin].get_element(Q,U));
          }
          if(U==P) {
            phi->accumulate_element(PQ,UV,-1.0*fmat[spin].get_element(Q,V));
          }

        }
      }
    }  // pairspin!=AlphaBeta
    return(phi);
  }

  RefSCMatrix PsiCorrWavefunction_PT2R12::phi_twoelec(SpinCase2 pairspin) {
    double energ_PT2R12[NSpinCases2];
    RefSCMatrix fmat[NSpinCases1];
    RefSymmSCMatrix opdm[NSpinCases1];
    const RefSymmSCMatrix tpdm = twopdm_refmo(pairspin);

    for(int i=0; i<NSpinCases1; i++) {
      SpinCase1 spin = static_cast<SpinCase1>(i);
      fmat[i] = f(spin);
      opdm[i] = onepdm_refmo(spin);
    }
    const int nmo = opdm[0].dim().n();

    double J = 0.0;
    for(int i=0; i<NSpinCases1; i++) {
      for(int z=0; z<nmo; z++) {
        for(int y=0; y<nmo; y++) {
          J += fmat[i].get_element(z,y)*opdm[i].get_element(y,z);
        }
      }
    }

    RefSCMatrix phi;
    phi = C(pairspin).kit()->matrix(r12eval_->dim_gg(pairspin),r12eval_->dim_gg(pairspin));
    phi.assign(0.0);
    phi.accumulate(tpdm*J);

    // redefined the sign for phi to correspond directly to the Fock operator
    phi.scale(-1.0);

    if(debug_>=DefaultPrintThresholds::mostO4) {
      phi->print(prepend_spincase(pairspin,"phi").c_str());
    }

    return(phi);
  }

//#define DEBUG_PHI
  // In this function, the statements commented out with "#if 0" would be more
  // efficient, but are not tested yet.
  RefSCMatrix PsiCorrWavefunction_PT2R12::phi(SpinCase2 pairspin) {
    RefSCMatrix fmat[NSpinCases1];
    RefSymmSCMatrix opdm[NSpinCases1];
    RefSymmSCMatrix tpdm[NSpinCases2];

    for(int i=0; i<NSpinCases1; i++) {
      SpinCase1 spin = static_cast<SpinCase1>(i);
      fmat[i] = f(spin);
      opdm[i] = onepdm_refmo(spin);
    }
    const int nmo = opdm[0].dim().n();

    // R12 intermediates, transformed with geminal amplitude matrix
    for(int i=0; i<NSpinCases2; i++) {
      SpinCase2 S = static_cast<SpinCase2>(i);
      tpdm[i] = twopdm_refmo(S);
    }

    // intermediates
    RefSCMatrix I[NSpinCases1];
    RefSCMatrix K[NSpinCases1];
    double J = 0.0;
    RefSCMatrix L[NSpinCases1];
    RefSCMatrix M[NSpinCases1];
    for(int i=0; i<NSpinCases1; i++) {
#if 1
      I[i] = fmat[i].clone();
      I[i].assign(0.0);
      for(int q=0; q<nmo; q++) {
        for(int v=0; v<nmo; v++) {
          for(int z=0; z<nmo; z++) {
            for(int y=0; y<nmo; y++) {
              I[i].accumulate_element(q,v,opdm[i].get_element(q,z)*fmat[i].get_element(z,y)*opdm[i].get_element(y,v));
            }
          }
        }
      }
#endif
#if 0
      I[i] = opdm[i] * fmat[i] * opdm[i];
#endif

#if 1
      for(int z=0; z<nmo; z++) {
        for(int y=0; y<nmo; y++) {
          J += fmat[i].get_element(z,y)*opdm[i].get_element(y,z);
        }
      }
#endif
#if 0
      J += (fmat[i] * opdm[i]).trace();
#endif

#if 1
      K[i] = fmat[i].clone();
      K[i].assign(0.0);
      for(int z=0; z<nmo; z++) {
        for(int u=0; u<nmo; u++) {
          for(int y=0; y<nmo; y++) {
            K[i].accumulate_element(z,u,opdm[i].get_element(y,u)*fmat[i].get_element(z,y));
          }
        }
      }
#endif
#if 0
      K[i] = opdm[i] * fmat[i].t();
#endif

#if 1
      L[i] = fmat[i].clone();
      L[i].assign(0.0);
      for(int p=0; p<nmo; p++) {
        for(int y=0; y<nmo; y++) {
          for(int z=0; z<nmo; z++) {
            L[i].accumulate_element(p,y,opdm[i].get_element(p,z)*fmat[i].get_element(z,y));
          }
        }
      }
#endif
#if 0
      L[i] = opdm[i] *  fmat[i];
#endif

      M[i] = fmat[i].clone();
      M[i].assign(0.0);
      if(i==Alpha) {
        for(int q=0; q<nmo; q++) {
          for(int y=0; y<nmo; y++) {
            for(int v=0; v<nmo; v++) {
              for(int z=0; z<nmo; z++) {
                int qy_alphabeta = q*nmo + y;
                int vz_alphabeta = v*nmo + z;
                M[i].accumulate_element(q,v,fmat[Beta].get_element(z,y)*tpdm[AlphaBeta].get_element(qy_alphabeta,vz_alphabeta));
                if((q!=y) && (v!=z)) {
                  int qy_alphaalpha = lowerupper_index(q,y);
                  int vz_alphaalpha = lowerupper_index(v,z);
                  int sign = indexsizeorder_sign(q,y)*indexsizeorder_sign(v,z);
                  M[i].accumulate_element(q,v,sign*fmat[Alpha].get_element(z,y)*tpdm[AlphaAlpha].get_element(qy_alphaalpha,vz_alphaalpha));
                }
              }
            }
          }
        }

      }  // i==Alpha
      else {  // i!=Alpha
        for(int q=0; q<nmo; q++) {
          for(int v=0; v<nmo; v++) {
            for(int z=0; z<nmo; z++) {
              for(int y=0; y<nmo; y++) {
                int yq_alphabeta = y*nmo + q;
                int zv_alphabeta = z*nmo + v;
                M[i].accumulate_element(q,v,fmat[Alpha].get_element(z,y)*tpdm[AlphaBeta].get_element(yq_alphabeta,zv_alphabeta));
                if((q!=y) && (v!=z)) {
                  int qy_betabeta = lowerupper_index(q,y);
                  int vz_betabeta = lowerupper_index(v,z);
                  int sign = indexsizeorder_sign(q,y)*indexsizeorder_sign(v,z);
                  M[i].accumulate_element(q,v,sign*fmat[Beta].get_element(z,y)*tpdm[BetaBeta].get_element(qy_betabeta,vz_betabeta));
                }
              }
            }
          }
        }
      }  // i!=Alpha

      const SpinCase1 spin = static_cast<SpinCase1>(i);
      L[i].print(prepend_spincase(spin,"K").c_str());
      I[i].print(prepend_spincase(spin,"I").c_str());
      M[i].print(prepend_spincase(spin,"M").c_str());

    }

    // computing phi
    RefSCMatrix phi;
    phi = C(pairspin).kit()->matrix(r12eval_->dim_gg(pairspin),r12eval_->dim_gg(pairspin));
    phi.assign(0.0);
#ifdef DEBUG_PHI
    RefSCMatrix term1 = phi.clone();
    RefSCMatrix term2 = phi.clone();
    RefSCMatrix term3 = phi.clone();
    RefSCMatrix term4 = phi.clone();
    RefSCMatrix term5 = phi.clone();
    term1->assign(0.0);
    term2->assign(0.0);
    term3->assign(0.0);
    term4->assign(0.0);
    term5->assign(0.0);
#endif
    if(pairspin==AlphaBeta) {
      SpinMOPairIter UV_iter(r12evalinfo_->refinfo()->orbs(Alpha),r12evalinfo_->refinfo()->orbs(Beta),pairspin);
      SpinMOPairIter PQ_iter(r12evalinfo_->refinfo()->orbs(Alpha),r12evalinfo_->refinfo()->orbs(Beta),pairspin);
      for(PQ_iter.start(); int(PQ_iter); PQ_iter.next()) {
        int P = PQ_iter.i();
        int Q = PQ_iter.j();
        int PQ = PQ_iter.ij();
        for(UV_iter.start(); int(UV_iter); UV_iter.next()) {
          int U = UV_iter.i();
          int V = UV_iter.j();
          int UV = UV_iter.ij();

          phi.accumulate_element(PQ,UV,2.0*opdm[Alpha].get_element(P,U)*I[Beta].get_element(Q,V)
                                      +2.0*opdm[Beta].get_element(Q,V)*I[Alpha].get_element(P,U));

          phi.accumulate_element(PQ,UV,-2.0*opdm[Alpha].get_element(P,U)*opdm[Beta].get_element(Q,V)*J);

          for(int Z=0; Z<PQ_iter.ni(); Z++) {
            int ZV = Z*nmo + V;
            int UZ = U*nmo + Z;
            phi.accumulate_element(PQ,UV,-K[Alpha].get_element(Z,U)*tpdm[AlphaBeta].get_element(PQ,ZV)
                                         -K[Beta].get_element(Z,V)*tpdm[AlphaBeta].get_element(PQ,UZ));
          }

          for(int Y=0; Y<PQ_iter.ni(); Y++) {
            int YQ = Y*nmo + Q;
            int PY = P*nmo + Y;
            phi.accumulate_element(PQ,UV,-L[Alpha].get_element(P,Y)*tpdm[AlphaBeta].get_element(YQ,UV)
                                         -L[Beta].get_element(Q,Y)*tpdm[AlphaBeta].get_element(PY,UV));
          }

          phi.accumulate_element(PQ,UV,opdm[Alpha].get_element(P,U)*M[Beta].get_element(Q,V)
                                      +opdm[Beta].get_element(Q,V)*M[Alpha].get_element(P,U));

        }
      }
    }  // pairspin==AlphaBeta
    else {  // pairspin!=AlphaBeta
      SpinCase1 spin;
      SpinCase1 otherspin;
      if(pairspin==AlphaAlpha) {
        spin = Alpha;
        otherspin = Beta;
      }
      else {
        spin = Beta;
        otherspin = Alpha;
      }
      SpinMOPairIter UV_iter(r12evalinfo_->refinfo()->orbs(spin),r12evalinfo_->refinfo()->orbs(spin),pairspin);
      SpinMOPairIter PQ_iter(r12evalinfo_->refinfo()->orbs(spin),r12evalinfo_->refinfo()->orbs(spin),pairspin);
      for(PQ_iter.start(); int(PQ_iter); PQ_iter.next()) {
        int P = PQ_iter.i();
        int Q = PQ_iter.j();
        int PQ = PQ_iter.ij();
        double sign_PQ = indexsizeorder_sign(P,Q);
        int QP = lowerupper_index(Q,P);
        double sign_QP = indexsizeorder_sign(Q,P);

        for(UV_iter.start(); int(UV_iter); UV_iter.next()) {
          int U = UV_iter.i();
          int V = UV_iter.j();
          int UV = UV_iter.ij();

          if((U!=V) && (P!=Q)) {
            phi.accumulate_element(PQ,UV,2.0*opdm[spin].get_element(P,U)*I[spin].get_element(Q,V)
                                        -2.0*opdm[spin].get_element(Q,U)*I[spin].get_element(P,V)
                                        -2.0*opdm[spin].get_element(P,V)*I[spin].get_element(Q,U)
                                        +2.0*opdm[spin].get_element(Q,V)*I[spin].get_element(P,U));
#ifdef DEBUG_PHI
            term1.accumulate_element(PQ,UV,2.0*opdm[spin].get_element(P,U)*I[spin].get_element(Q,V)
                                        -2.0*opdm[spin].get_element(Q,U)*I[spin].get_element(P,V)
                                        -2.0*opdm[spin].get_element(P,V)*I[spin].get_element(Q,U)
                                        +2.0*opdm[spin].get_element(Q,V)*I[spin].get_element(P,U));
#endif

            phi.accumulate_element(PQ,UV,opdm[spin].get_element(P,V)*opdm[spin].get_element(Q,U)*J
                                        -opdm[spin].get_element(Q,V)*opdm[spin].get_element(P,U)*J
                                        -opdm[spin].get_element(P,U)*opdm[spin].get_element(Q,V)*J
                                        +opdm[spin].get_element(Q,U)*opdm[spin].get_element(P,V)*J);
#ifdef DEBUG_PHI
            term2.accumulate_element(PQ,UV,opdm[spin].get_element(P,V)*opdm[spin].get_element(Q,U)*J
                                          -opdm[spin].get_element(Q,V)*opdm[spin].get_element(P,U)*J
                                          -opdm[spin].get_element(P,U)*opdm[spin].get_element(Q,V)*J
                                          +opdm[spin].get_element(Q,U)*opdm[spin].get_element(P,V)*J);
#endif

            for(int Z=0; Z<UV_iter.ni(); Z++) {
              if(V!=Z) {
                int VZ = lowerupper_index(V,Z);
                double sign_VZ = indexsizeorder_sign(V,Z);
                phi.accumulate_element(PQ,UV,sign_VZ*K[spin].get_element(Z,U)*tpdm[pairspin].get_element(PQ,VZ));
#ifdef DEBUG_PHI
                term3.accumulate_element(PQ,UV,sign_VZ*K[spin].get_element(Z,U)*tpdm[pairspin].get_element(PQ,VZ));
#endif
              }
              if(U!=Z) {
                int UZ = lowerupper_index(U,Z);
                double sign_UZ = indexsizeorder_sign(U,Z);
                phi.accumulate_element(PQ,UV,-sign_UZ*K[spin].get_element(Z,V)*tpdm[pairspin].get_element(PQ,UZ));
#ifdef DEBUG_PHI
                term3.accumulate_element(PQ,UV,-sign_UZ*K[spin].get_element(Z,V)*tpdm[pairspin].get_element(PQ,UZ));
#endif
              }
            }

            for(int Y=0; Y<PQ_iter.ni(); Y++) {
              double sign_UV = indexsizeorder_sign(U,V);
              int VU = lowerupper_index(V,U);
              double sign_VU = indexsizeorder_sign(V,U);
              if(Q!=Y) {
                int QY = lowerupper_index(Q,Y);
                double sign_QY = indexsizeorder_sign(Q,Y);
                phi.accumulate_element(PQ,UV,sign_QY*L[spin].get_element(P,Y)*tpdm[pairspin].get_element(QY,UV));
#ifdef DEBUG_PHI
                term4.accumulate_element(PQ,UV,sign_QY*L[spin].get_element(P,Y)*tpdm[pairspin].get_element(QY,UV));
#endif
              }
              if(P!=Y) {
                int PY = lowerupper_index(P,Y);
                double sign_PY = indexsizeorder_sign(P,Y);
                phi.accumulate_element(PQ,UV,-sign_PY*L[spin].get_element(Q,Y)*tpdm[pairspin].get_element(PY,UV));
#ifdef DEBUG_PHI
                term4.accumulate_element(PQ,UV,-sign_PY*L[spin].get_element(Q,Y)*tpdm[pairspin].get_element(PY,UV));
#endif
              }
            }

            phi.accumulate_element(PQ,UV,opdm[spin].get_element(P,U)*M[spin].get_element(Q,V)
                                        -opdm[spin].get_element(Q,U)*M[spin].get_element(P,V)
                                        -opdm[spin].get_element(P,V)*M[spin].get_element(Q,U)
                                        +opdm[spin].get_element(Q,V)*M[spin].get_element(P,U));
#ifdef DEBUG_PHI
            term5.accumulate_element(PQ,UV,opdm[spin].get_element(P,U)*M[spin].get_element(Q,V)
                                        -opdm[spin].get_element(Q,U)*M[spin].get_element(P,V)
                                        -opdm[spin].get_element(P,V)*M[spin].get_element(Q,U)
                                        +opdm[spin].get_element(Q,V)*M[spin].get_element(P,U));
#endif

          }
        }
      }
#ifdef DEBUG_PHI
      term1->print(prepend_spincase(pairspin,"term1 of phi").c_str());
      term2->print(prepend_spincase(pairspin,"term2 of phi").c_str());
      term3->print(prepend_spincase(pairspin,"term3 of phi").c_str());
      term4->print(prepend_spincase(pairspin,"term4 of phi").c_str());
      term5->print(prepend_spincase(pairspin,"term5 of phi").c_str());
      RefSCMatrix term_sum = term1 + term2 + term3 + term4 + term5;
      term_sum->print(prepend_spincase(pairspin,"term_sum of phi").c_str());
#endif
    }  // pairspin!=AlphaBeta

    // redefined sign so that phi corresponds to the Fock matrix
    phi.scale(-1.0);

    if(debug_>=DefaultPrintThresholds::mostO4) {
      phi->print(prepend_spincase(pairspin,"phi").c_str());
    }

    return(phi);
  }

  RefSymmSCMatrix PsiCorrWavefunction_PT2R12::phi_cumulant(SpinCase2 spin12) {
    RefSCMatrix fmat[NSpinCases1];
    RefSymmSCMatrix opdm[NSpinCases1];
    RefSymmSCMatrix tpcm[NSpinCases2];

    for(int s=0; s<NSpinCases1; s++) {
      SpinCase1 spin = static_cast<SpinCase1>(s);
      fmat[s] = f(spin);
      opdm[s] = onepdm_refmo(spin);
    }
    const int nmo = opdm[0].dim().n();

    // R12 intermediates, transformed with geminal amplitude matrix
    for(int i=0; i<NSpinCases2; i++) {
      SpinCase2 S = static_cast<SpinCase2>(i);
      tpcm[i] = lambda_refmo(S);
    }

    // precompute intermediates
    RefSCMatrix K[NSpinCases1];   // \gamma f
    RefSCMatrix I[NSpinCases1];   // \gamma f \gamma
    RefSCMatrix M[NSpinCases1];   // trace(f lambda)
    for(int i=0; i<NSpinCases1; i++) {
      K[i] = fmat[i].clone();
      K[i].assign(0.0);
      for(int u=0; u<nmo; u++) {
        for(int z=0; z<nmo; z++) {
          for(int y=0; y<nmo; y++) {
            K[i].accumulate_element(u,z,opdm[i].get_element(u,y)*fmat[i].get_element(y,z));
          }
        }
      }
      I[i] = fmat[i].clone();
      I[i].assign(0.0);
      for(int q=0; q<nmo; q++) {
        for(int v=0; v<nmo; v++) {
          for(int z=0; z<nmo; z++) {
              I[i].accumulate_element(q,v,K[i].get_element(q,z)*opdm[i].get_element(z,v));
          }
        }
      }
    }
#if 0
    for(int s=0; s<NSpinCases1; ++s) {
      K[s] = opdm[s] * fmat[s];
      I[s] = K[s] * opdm[s];
    }
#endif

    M[Alpha] = fmat[Alpha].clone(); M[Alpha].assign(0.0);
    M[Beta] = fmat[Beta].clone(); M[Beta].assign(0.0);
    // each M has contributions from 2 pairspincases
    // AlphaAlpha contributes to Alpha only
    for(int P=0; P<nmo; ++P) {
      for(int Q=0; Q<nmo; ++Q) {
        int PQ_aa, sign_PQ=+1;
        if (P > Q) {
          PQ_aa = P*(P-1)/2 + Q;
        }
        else {
          PQ_aa = Q*(Q-1)/2 + P;
          sign_PQ = -1;
        }
        const int PQ_ab = P * nmo + Q;
        const int QP_ab = Q * nmo + P;
        for(int U=0; U<nmo; ++U) {
          for(int V=0; V<nmo; ++V) {
            int UV_aa, sign_UV=+1;
            if (U > V) {
              UV_aa = U*(U-1)/2 + V;
            }
            else {
              UV_aa = V*(V-1)/2 + U;
              sign_UV = -1;
            }
            const int UV_ab = U * nmo + V;
            const int VU_ab = V * nmo + U;

            double M_PU_a = 0.0;
            double M_PU_b = 0.0;
            if (P!=Q && U!=V) {
              M_PU_a += sign_PQ * sign_UV * fmat[Alpha].get_element(Q, V) * tpcm[AlphaAlpha].get_element(PQ_aa, UV_aa);
              M_PU_b += sign_PQ * sign_UV * fmat[Beta].get_element(Q, V) * tpcm[BetaBeta].get_element(PQ_aa, UV_aa);
            }
            M_PU_a += fmat[Beta].get_element(Q, V) * tpcm[AlphaBeta].get_element(PQ_ab, UV_ab);
            M_PU_b += fmat[Alpha].get_element(Q, V) * tpcm[AlphaBeta].get_element(QP_ab, VU_ab);
            M[Alpha].accumulate_element(P, U, M_PU_a);
            M[Beta].accumulate_element(P, U, M_PU_b);
          }
        }
      }
    }

    for(int s=0; s<NSpinCases1; ++s) {
      const SpinCase1 spin = static_cast<SpinCase1>(s);
      K[s].print(prepend_spincase(spin,"K new").c_str());
      I[s].print(prepend_spincase(spin,"I new").c_str());
      M[s].print(prepend_spincase(spin,"M new").c_str());
    }

    // compute phi:
    // \phi^{u v}_{p q} = & P(p q) P(u v) (\gamma^u_p \gamma^{q_3}_q f^{q_2}_{q_3} \gamma^v_{q_2}
    //                                     + \frac{1}{2} \gamma^v_{q_2} f^{q_2}_{q_3} \lambda^{u q_3}_{p q}
    //                                     + \frac{1}{2} \gamma^{q_2}_p f^{q_3}_{q_2} \lambda^{u v}_{q_3 q}
    //                                     - \gamma^u_p f^{q_2}_{q_3} \lambda^{q_3 v}_{q_2 q})
    RefSymmSCMatrix phi = tpcm[spin12].clone();
    phi.assign(0.0);
    const SpinCase1 spin1 = case1(spin12);
    const SpinCase1 spin2 = case2(spin12);

    if (spin12 == AlphaBeta) {
      SpinMOPairIter UV_iter(r12evalinfo_->refinfo()->orbs(Alpha),r12evalinfo_->refinfo()->orbs(Beta),spin12);
      SpinMOPairIter PQ_iter(r12evalinfo_->refinfo()->orbs(Alpha),r12evalinfo_->refinfo()->orbs(Beta),spin12);
      for(PQ_iter.start(); int(PQ_iter); PQ_iter.next()) {
        int P = PQ_iter.i();
        int Q = PQ_iter.j();
        int PQ = PQ_iter.ij();
        for(UV_iter.start(); int(UV_iter); UV_iter.next()) {
          int U = UV_iter.i();
          int V = UV_iter.j();
          int UV = UV_iter.ij();

          // first term is easy
          double phi_PQ_UV = (   I[Alpha].get_element(P, U) * opdm[Beta].get_element(Q, V) +
                              opdm[Alpha].get_element(P, U) *    I[Beta].get_element(Q, V)
                             );

          // second and third terms contain contributions from the cumulant of same spin
          for(int Q3=0; Q3<nmo; ++Q3) {
            int UQ3 = U * nmo + Q3;
            int Q3V = Q3 * nmo + V;
            int PQ3 = P * nmo + Q3;
            int Q3Q = Q3 * nmo + Q;
            // +1 instead of +1/2 because there are 2 permutation operators!
            phi_PQ_UV += (tpcm[AlphaBeta].get_element(PQ, UQ3) * K[Beta].get_element(V, Q3) +
                          tpcm[AlphaBeta].get_element(PQ, Q3V) * K[Alpha].get_element(U, Q3) +
                          tpcm[AlphaBeta].get_element(PQ3, UV) * K[Beta].get_element(Q, Q3) +
                          tpcm[AlphaBeta].get_element(Q3Q, UV) * K[Alpha].get_element(P, Q3)
                         );
          }

          // fourth term is also easy
          phi_PQ_UV -= (   M[Alpha].get_element(P, U) * opdm[Beta].get_element(Q, V) +
                        opdm[Alpha].get_element(P, U) *    M[Beta].get_element(Q, V)
                       );

          phi(PQ,UV) = phi_PQ_UV;
        }
      }
    }
    else if (spin12 == AlphaAlpha || spin12 == BetaBeta) {
      const SpinCase1 spin = spin1;
      SpinMOPairIter UV_iter(r12evalinfo_->refinfo()->orbs(spin),r12evalinfo_->refinfo()->orbs(spin),spin12);
      SpinMOPairIter PQ_iter(r12evalinfo_->refinfo()->orbs(spin),r12evalinfo_->refinfo()->orbs(spin),spin12);
      for(PQ_iter.start(); int(PQ_iter); PQ_iter.next()) {
        int P = PQ_iter.i();
        int Q = PQ_iter.j();
        int PQ = PQ_iter.ij();
        for(UV_iter.start(); int(UV_iter); UV_iter.next()) {
          int U = UV_iter.i();
          int V = UV_iter.j();
          int UV = UV_iter.ij();

          // first term is easy
          double phi_PQ_UV = (   I[spin].get_element(P, U) * opdm[spin].get_element(Q, V) +
                              opdm[spin].get_element(P, U) *    I[spin].get_element(Q, V) -
                                 I[spin].get_element(P, V) * opdm[spin].get_element(Q, U) -
                              opdm[spin].get_element(P, V) *    I[spin].get_element(Q, U)
                             );

          // second and third terms contain contributions from the cumulant of same spin
          // +1 instead of +1/2 because there are 2 permutation operators!
          // V
          for(int Q3=0; Q3<U; ++Q3) {
            const int UQ3 = U * (U-1)/2 + Q3;
            phi_PQ_UV += tpcm[spin12].get_element(PQ, UQ3) * K[spin].get_element(V, Q3);
          }
          for(int Q3=U+1; Q3<nmo; ++Q3) {
            const int Q3U = Q3 * (Q3-1)/2 + U;
            // - sign!
            phi_PQ_UV -= tpcm[spin12].get_element(PQ, Q3U) * K[spin].get_element(V, Q3);
          }
          // U
          for(int Q3=0; Q3<V; ++Q3) {
            const int VQ3 = V * (V-1)/2 + Q3;
            // - sign!
            phi_PQ_UV -= tpcm[spin12].get_element(PQ, VQ3) * K[spin].get_element(U, Q3);
          }
          for(int Q3=V+1; Q3<nmo; ++Q3) {
            const int Q3V = Q3 * (Q3-1)/2 + V;
            phi_PQ_UV += tpcm[spin12].get_element(PQ, Q3V) * K[spin].get_element(U, Q3);
          }
          // Q
          for(int Q3=0; Q3<P; ++Q3) {
            const int PQ3 = P * (P-1)/2 + Q3;
            phi_PQ_UV += tpcm[spin12].get_element(PQ3, UV) * K[spin].get_element(Q, Q3);
          }
          for(int Q3=P+1; Q3<nmo; ++Q3) {
            const int Q3P = Q3 * (Q3-1)/2 + P;
            // - sign!
            phi_PQ_UV -= tpcm[spin12].get_element(Q3P, UV) * K[spin].get_element(Q, Q3);
          }
          // P
          for(int Q3=0; Q3<Q; ++Q3) {
            const int QQ3 = Q * (Q-1)/2 + Q3;
            // - sign!
            phi_PQ_UV -= tpcm[spin12].get_element(QQ3, UV) * K[spin].get_element(P, Q3);
          }
          for(int Q3=Q+1; Q3<nmo; ++Q3) {
            const int Q3Q = Q3 * (Q3-1)/2 + Q;
            // -1 instead of -1/2 because there are 2 permutation operators!
            phi_PQ_UV += tpcm[spin12].get_element(Q3Q, UV) * K[spin].get_element(P, Q3);
          }

          // fourth term is also easy
          phi_PQ_UV -= (    M[spin].get_element(P, U) * opdm[spin].get_element(Q, V) +
                         opdm[spin].get_element(P, U) *    M[spin].get_element(Q, V) -
                            M[spin].get_element(P, V) * opdm[spin].get_element(Q, U) -
                         opdm[spin].get_element(P, V) *    M[spin].get_element(Q, U)
                       );

          phi(PQ,UV) = phi_PQ_UV;
        }
      }
    }
    else { // invalid spin12
      assert(false);
    }

    if(debug_>=DefaultPrintThresholds::mostO4) {
      phi->print(prepend_spincase(spin12,"phi (new)").c_str());
    }

    return phi;
  }

  double PsiCorrWavefunction_PT2R12::energy_PT2R12_projector1(SpinCase2 pairspin) {
    const int nelectron = reference_->nelectron();
    SpinCase1 spin1 = case1(pairspin);
    SpinCase1 spin2 = case2(pairspin);
    Ref<OrbitalSpace> gg1space = r12eval_->ggspace(spin1);
    Ref<OrbitalSpace> gg2space = r12eval_->ggspace(spin2);
    SpinMOPairIter gg_iter(gg1space,gg2space,pairspin);

    RefSymmSCMatrix tpdm = twopdm_refmo(pairspin);
    RefSCMatrix Phi;
    if(exactgamma_phi_twoelec_ && (nelectron==2)) {
      Phi = phi_twoelec(pairspin);
    }
    else {
      Phi = phi(pairspin);
    }
    RefSCMatrix VT = V_transformed_by_C(pairspin);
    RefSymmSCMatrix TXT = X_transformed_by_C(pairspin);
    RefSymmSCMatrix TBT = B_transformed_by_C(pairspin);


    ExEnv::out0() << "pairspin "
                  << ((pairspin==AlphaBeta) ? "AlphaBeta" : ((pairspin==AlphaAlpha) ? "AlphaAlpha" : "BetaBeta")) << endl;

    // printing pair energies
    RefSCMatrix VT_t_tpdm = 2.0*VT*tpdm;
    RefSCMatrix TBT_t_tpdm = TBT*tpdm;
    RefSCMatrix TXT_t_Phi = TXT*Phi;
    RefSCMatrix HylleraasMatrix = VT_t_tpdm + TBT_t_tpdm - TXT_t_Phi;

    ExEnv::out0() << "pair energies:" << endl;
    for(gg_iter.start(); int(gg_iter); gg_iter.next()) {
      int i = gg_iter.i();
      int j = gg_iter.j();
      int ij = gg_iter.ij();
      ExEnv::out0() << setw(6) << i << setw(6) << j
                    << setw(20) << setprecision(12) << HylleraasMatrix.get_element(ij,ij) << endl;
    }

    double energy = 0.0;
    if(nfzc_==0) {
      ExEnv::out0() << "Computing energy as trace of HylleraasMatrix." << endl;
      energy = HylleraasMatrix.trace();
    }
    else {
      ExEnv::out0() << "Computing energy as sum of pair energies." << endl;
      for(gg_iter.start(); int(gg_iter); gg_iter.next()) {
        int i = gg_iter.i();
        int j = gg_iter.j();
        int ij = gg_iter.ij();
        if((i>=nfzc_) && (j>=nfzc_)) {
          energy+=HylleraasMatrix.get_element(ij,ij);
        }
      }
    }
    ExEnv::out0() << "energy " << setprecision(12) << energy << endl;

    return(energy);
  }

  double PsiCorrWavefunction_PT2R12::energy_PT2R12_projector2(SpinCase2 pairspin) {
    const int nelectron = reference_->nelectron();
    SpinCase1 spin1 = case1(pairspin);
    SpinCase1 spin2 = case2(pairspin);
    Ref<OrbitalSpace> gg1space = r12eval_->ggspace(spin1);
    Ref<OrbitalSpace> gg2space = r12eval_->ggspace(spin2);
    SpinMOPairIter gg_iter(gg1space,gg2space,pairspin);

    RefSymmSCMatrix TBT = B_transformed_by_C(pairspin);
    RefSymmSCMatrix tpdm = twopdm_refmo(pairspin);
    RefSCMatrix HylleraasMatrix = TBT*tpdm;
    TBT=0;
    tpdm=0;

    RefSymmSCMatrix TXT = X_transformed_by_C(pairspin);
    if(exactgamma_phi_twoelec_ && (nelectron==2)) { // cannot compute exact phi for 2-e systems
      assert(false);
    }
    RefSymmSCMatrix Phi = phi_cumulant(pairspin);
    RefSCMatrix TXT_t_Phi = TXT*Phi; TXT_t_Phi.scale(-1.0);
    HylleraasMatrix.accumulate(TXT_t_Phi);
    TXT=0;
    Phi=0;
    TXT_t_Phi=0;

    RefSCMatrix V_genref = V_genref_projector2(pairspin);
    RefSCMatrix T = C(pairspin);
    HylleraasMatrix.accumulate(2.0*V_genref.t()*T);
    T=0;
    V_genref=0;

    ExEnv::out0() << "pair energies:" << endl;
    for(gg_iter.start(); int(gg_iter); gg_iter.next()) {
      int i = gg_iter.i();
      int j = gg_iter.j();
      int ij = gg_iter.ij();
      ExEnv::out0() << setw(6) << i << setw(6) << j
                    << setw(20) << setprecision(12) << HylleraasMatrix.get_element(ij,ij) << endl;
    }

    double energy = 0.0;
    if(nfzc_==0) {
      ExEnv::out0() << "Computing energy as trace of HylleraasMatrix." << endl;
      energy = HylleraasMatrix.trace();
    }
    else {
      ExEnv::out0() << "Computing energy as sum of pair energies." << endl;
      for(gg_iter.start(); int(gg_iter); gg_iter.next()) {
        int i = gg_iter.i();
        int j = gg_iter.j();
        int ij = gg_iter.ij();
        if((i>=nfzc_) && (j>=nfzc_)) {
          energy+=HylleraasMatrix.get_element(ij,ij);
        }
      }
    }
    ExEnv::out0() << "energy " << setprecision(12) << energy << endl;

    return(energy);
  }

  void PsiCorrWavefunction_PT2R12::print_pair_energies(const RefSCMatrix &Contrib_mat,SpinCase2 pairspin) {
    SpinCase1 spin1 = case1(pairspin);
    SpinCase1 spin2 = case2(pairspin);
    Ref<OrbitalSpace> gg1space = r12eval_->ggspace(spin1);
    Ref<OrbitalSpace> gg2space = r12eval_->ggspace(spin2);
    SpinMOPairIter gg_iter(gg1space,gg2space,pairspin);
    double energy = 0.0;
    ExEnv::out0() << "pair energies:" << endl;
    for(gg_iter.start(); int(gg_iter); gg_iter.next()) {
      int i = gg_iter.i();
      int j = gg_iter.j();
      int ij = gg_iter.ij();
      ExEnv::out0() << setw(6) << i << setw(6) << j
                    << setw(20) << setprecision(12) << Contrib_mat.get_element(ij,ij) << endl;
      if((i>=nfzc_) && (j>=nfzc_)) {
        energy+=Contrib_mat.get_element(ij,ij);
      }
    }
    ExEnv::out0() << "total contribution: " << setprecision(12) << energy << endl;
  }


//////////////////////////////////////////////////////////////////////////

}

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:

//
// mbpt.cc
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Ida Nielsen <ida@kemi.aau.dk>
// Maintainer: LPS
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

#include <util/misc/scexception.h>
#include <util/misc/formio.h>
#include <util/misc/exenv.h>
#include <util/state/stateio.h>
#include <math/scmat/blocked.h>
#include <chemistry/qc/basis/petite.h>
#include <chemistry/qc/mbpt/mbpt.h>

using namespace std;
using namespace sc;

/////////////////////////////////////////////////////////////////
// Function dquicksort performs a quick sort (smaller -> larger)
// of the double data in item by the integer indices in index;
// data in item remain unchanged
/////////////////////////////////////////////////////////////////
static void
dqs(double *item,int *index,int left,int right)
{
  int i,j;
  double x;
  int y;

  i=left; j=right;
  x=item[index[(left+right)/2]];

  do {
    while(item[index[i]]<x && i<right) i++;
    while(x<item[index[j]] && j>left) j--;

    if (i<=j) {
      if (item[index[i]] != item[index[j]]) {
        y=index[i];
        index[i]=index[j];
        index[j]=y;
        }
      i++; j--;
      }
    } while(i<=j);

  if (left<j) dqs(item,index,left,j);
  if (i<right) dqs(item,index,i,right);
}
static void
dquicksort(double *item,int *index,int n)
{
  int i;
  if (n<=0) return;
  for (i=0; i<n; i++) {
    index[i] = i;
    }
  dqs(item,index,0,n-1);
  }

///////////////////////////////////////////////////////////////////////////
// MBPT2

static ClassDesc MBPT2_cd(
  typeid(MBPT2),"MBPT2",9,"public Wavefunction",
  0, create<MBPT2>, create<MBPT2>);

MBPT2::MBPT2(StateIn& s):
  SavableState(s),
  Wavefunction(s)
{
  reference_ << SavableState::restore_state(s);
  s.get(nfzc);
  s.get(nfzv);
  if (s.version(::class_desc<MBPT2>()) >= 8) {
      double dmem_alloc;
      s.get(dmem_alloc);
      mem_alloc = size_t(dmem_alloc);
    }
  else {
      unsigned int imem_alloc;
      s.get(imem_alloc);
      mem_alloc = imem_alloc;
    }
  s.get(method_);
  s.get(algorithm_);

  if (s.version(::class_desc<MBPT2>()) <= 6) {
      int debug_old;
      s.get(debug_old);
    }

  if (s.version(::class_desc<MBPT2>()) >= 2) {
      s.get(do_d1_);
    }
  else {
      do_d1_ = 0;
    }

  if (s.version(::class_desc<MBPT2>()) >= 3) {
      s.get(dynamic_);
    }
  else {
      dynamic_ = 0;
    }

  if (s.version(::class_desc<MBPT2>()) >= 9) {
      s.get(print_percent_);
    }
  else {
      print_percent_ = 10.0;
    }

  if (s.version(::class_desc<MBPT2>()) >= 4) {
      s.get(cphf_epsilon_);
    }
  else {
      cphf_epsilon_ = 1.0e-8;
    }

  if (s.version(::class_desc<MBPT2>()) >= 5) {
      s.get(max_norb_);
    }
  else {
      max_norb_ = 0;
    }

  if (s.version(::class_desc<MBPT2>()) >= 6) {
      s.get(do_d2_);
    }
  else {
      do_d2_ = 1;
    }

  hf_energy_ = 0.0;

  symorb_irrep_ = 0;
  symorb_num_ = 0;

  eliminate_in_gmat_ = 1;

  mem = MemoryGrp::get_default_memorygrp();
  msg_ = MessageGrp::get_default_messagegrp();
  thr_ = ThreadGrp::get_default_threadgrp();

  restart_ecorr_ = 0.0;
  restart_orbital_v1_ = 0;
  restart_orbital_memgrp_ = 0;
}

MBPT2::MBPT2(const Ref<KeyVal>& keyval):
  Wavefunction(keyval)
{
  reference_ << keyval->describedclassvalue("reference");
  if (reference_.null()) {
      ExEnv::err0() << "MBPT2::MBPT2: no reference wavefunction" << endl;
      abort();
    }
  copy_orthog_info(reference_);
  nfzc = keyval->intvalue("nfzc");
  std::string nfzc_charval = keyval->stringvalue("nfzc");
  if (nfzc_charval == "auto") {
      if (molecule()->max_z() > 30) {
          ExEnv::err0()
               << "MBPT2: cannot use \"nfzc = auto\" for Z > 30" << endl;
          abort();
        }
      nfzc = molecule()->n_core_electrons()/2;
      ExEnv::out0() << indent
           << "MBPT2: auto-freezing " << nfzc << " core orbitals" << endl;
    }
  nfzv = keyval->intvalue("nfzv");
  mem_alloc = keyval->sizevalue("memory");
  if (keyval->error() != KeyVal::OK) {
      // by default, take half of the memory
      mem_alloc = ExEnv::memory()/2;
      if (mem_alloc == 0) {
          mem_alloc = DEFAULT_MPQC_MEMORY;
        }
    }
  mem << keyval->describedclassvalue("memorygrp");
  if (mem.null()) {
      mem = MemoryGrp::get_default_memorygrp();
    }
  msg_ = MessageGrp::get_default_messagegrp();
  thr_ = ThreadGrp::get_default_threadgrp();

  method_ = keyval->stringvalue("method", KeyValValuestring("default"));

  algorithm_ = keyval->stringvalue("algorithm");

  do_d1_ = keyval->booleanvalue("compute_d1");
  do_d2_ = keyval->booleanvalue("compute_d2",KeyValValueboolean(1));

  restart_ecorr_ = keyval->doublevalue("restart_ecorr");
  restart_orbital_v1_ = keyval->intvalue("restart_orbital_v1");
  restart_orbital_memgrp_ = keyval->intvalue("restart_orbital_memgrp");

  KeyValValueint default_dynamic(0);
  dynamic_ = keyval->booleanvalue("dynamic", default_dynamic);

  KeyValValuedouble default_print_percent(10.0);
  print_percent_ = keyval->doublevalue("print_percent", default_print_percent);

  cphf_epsilon_ = keyval->doublevalue("cphf_epsilon",KeyValValuedouble(1.e-8));

  max_norb_ = keyval->intvalue("max_norb",KeyValValueint(-1));

  hf_energy_ = 0.0;

  symorb_irrep_ = 0;
  symorb_num_ = 0;

  eliminate_in_gmat_ = 1;
}

MBPT2::~MBPT2()
{
  delete[] symorb_irrep_;
  delete[] symorb_num_;
}

void
MBPT2::save_data_state(StateOut& s)
{
  Wavefunction::save_data_state(s);
  SavableState::save_state(reference_.pointer(),s);
  s.put(nfzc);
  s.put(nfzv);
  double dmem_alloc = mem_alloc;
  s.put(dmem_alloc);
  s.put(method_);
  s.put(algorithm_);
  s.put(do_d1_);
  s.put(dynamic_);
  s.put(print_percent_);
  s.put(cphf_epsilon_);
  s.put(max_norb_);
  s.put(do_d2_);
}

void
MBPT2::print(ostream&o) const
{
  o << indent << "MBPT2:" << endl;
  o << incindent;
  Wavefunction::print(o);
  o << indent << "Reference Wavefunction:" << endl;
  o << incindent; reference_->print(o); o << decindent << endl;
  o << decindent;
}

//////////////////////////////////////////////////////////////////////////////

double
MBPT2::magnetic_moment() const
{
  return reference_->magnetic_moment();
}

RefSymmSCMatrix
MBPT2::density()
{
  return 0;
}

//////////////////////////////////////////////////////////////////////////////

void
MBPT2::compute()
{
  std::string reference_name(reference_->integral()->class_name());
  std::string integral_name(integral()->class_name());
  if ( reference_name != integral_name ) {
      FeatureNotImplemented ex(
          "cannot use a reference with a different Integral specialization",
          __FILE__, __LINE__, class_desc());
      try {
          ex.elaborate()
              << "reference uses " << reference_->integral()->class_name()
              << " but this object uses " << integral()->class_name()
              << std::endl;
        }
      catch (...) {}
      throw ex;
    }

  init_variables();

  reference_->set_desired_value_accuracy(desired_value_accuracy()
                                         / ref_to_mp2_acc());
  if (gradient_needed()) {
      if (nsocc) {
          ExEnv::errn() << "MBPT2: cannot compute open shell gradients" << endl;
          abort();
        }
      compute_cs_grad();
    }
  else {
      if (nsocc && algorithm_ == "memgrp") {
          ExEnv::errn() << "MBPT2: memgrp algorithm cannot compute open shell energy"
               << endl;
          abort();
        }
      if (!nsocc && (algorithm_.empty() || algorithm_ == "memgrp")) {
          compute_cs_grad();
        }
      else if ((algorithm_.empty() && msg_->n() <= 32)
               || (!algorithm_.empty() && algorithm_ == "v1")) {
          compute_hsos_v1();
        }
      else if (algorithm_.empty() || algorithm_ == "v2") {
          compute_hsos_v2();
        }
      else if (algorithm_ == "v2lb") {
          compute_hsos_v2_lb();
        }
      else {
          ExEnv::errn() << "MBPT2: unknown algorithm: " << algorithm_ << endl;
          abort();
        }
    }
}

//////////////////////////////////////////////////////////////////////////////

void
MBPT2::obsolete()
{
  // Solaris 2.7 workshop 5.0 is causing this routine to
  // be incorrectly called in a base class CTOR.  Thus
  // reference_ might be null and it must be tested.
  if (reference_.nonnull()) reference_->obsolete();
  Wavefunction::obsolete();
}

//////////////////////////////////////////////////////////////////////////////

bool
MBPT2::analytic_gradient_implemented() const
{
  int nb = reference_->oso_dimension()->n();
  int n = 0;

  for (int i=0; i<nb; i++) {
    if (reference_->occupation(i) == 1.0) n++;
    }

  if (n) return false;
  return true;
}

//////////////////////////////////////////////////////////////////////////////

int
MBPT2::value_implemented() const
{
  return 1;
}

//////////////////////////////////////////////////////////////////////////////

void
MBPT2::set_desired_value_accuracy(double acc)
{
  Function::set_desired_value_accuracy(acc);
  reference_->set_desired_value_accuracy(acc / ref_to_mp2_acc());
}

//////////////////////////////////////////////////////////////////////////////

void
MBPT2::eigen(RefDiagSCMatrix &vals, RefSCMatrix &vecs, RefDiagSCMatrix &occs)
{
  int i, j;
  if (nsocc) {
      if (reference_->n_fock_matrices() != 2) {
          ExEnv::errn() << "MBPT2: given open reference with"
               << " wrong number of Fock matrices" << endl;
          abort();
        }

      // Notation: oo = orthonormal symmetry orbital basis
      //           ao = atomic orbital basis
      //           so = symmetrized atomic orbital basis
      //           mo1 = SCF molecular orbital basis
      //           mo2 = MBPT molecular orbital basis

      // get the closed shell and open shell AO fock matrices
      RefSymmSCMatrix fock_c_so = reference_->fock(0);
      RefSymmSCMatrix fock_o_so = reference_->fock(1);

      // transform the AO fock matrices to the MO basis
      RefSymmSCMatrix fock_c_mo1
          = basis_matrixkit()->symmmatrix(oso_dimension());
      RefSymmSCMatrix fock_o_mo1
          = basis_matrixkit()->symmmatrix(oso_dimension());
      RefSCMatrix vecs_so_mo1 = reference_->eigenvectors();

      fock_c_mo1.assign(0.0);
      fock_o_mo1.assign(0.0);
      fock_c_mo1.accumulate_transform(vecs_so_mo1.t(), fock_c_so);
      fock_o_mo1.accumulate_transform(vecs_so_mo1.t(), fock_o_so);
      fock_c_so = 0;
      fock_o_so = 0;

     /* Convert to the Guest & Saunders general form.
        This is the form used for an OPT2 calculation.

            C        O         V
        ----------
        |        |
     C  |   fc   |
        |        |
        -------------------
        |        |        |
     O  | 2fc-fo |   fc   |
        |        |        |
        ----------------------------
        |        |        |        |
     V  |   fc   |   fo   |   fc   |
        |        |        |        |
        ----------------------------
      */
      RefSymmSCMatrix fock_eff_mo1 = fock_c_mo1.clone();
      fock_eff_mo1.assign(fock_c_mo1);
      for (i=0; i<oso_dimension()->n(); i++) {
          double occi = reference_->occupation(i);
          for (j=0; j<=i; j++) {
              double occj = reference_->occupation(j);
              if ((occi == 2.0 && occj == 1.0)
                  || (occi == 1.0 && occj == 2.0)) {
                  fock_eff_mo1.accumulate_element(i,j,
                                                  fock_c_mo1(i,j)-fock_o_mo1(i,j));
                }
              else if ((occi == 0.0 && occj == 1.0)
                  || (occi == 1.0 && occj == 0.0)) {
                  fock_eff_mo1.accumulate_element(i,j,
                                                  fock_o_mo1(i,j)-fock_c_mo1(i,j));
                }
            }
        }

      // diagonalize the new fock matrix
      RefDiagSCMatrix vals_mo2(fock_eff_mo1.dim(), fock_eff_mo1.kit());
      RefSCMatrix vecs_mo1_mo2(fock_eff_mo1.dim(), fock_eff_mo1.dim(),
                               fock_eff_mo1.kit());
      fock_eff_mo1.diagonalize(vals_mo2, vecs_mo1_mo2);
      vals = vals_mo2;

      // compute the AO to new MO scf vector
      RefSCMatrix so_ao = reference_->integral()->petite_list()->sotoao();

      if (debug_ > 1) {
          vecs_mo1_mo2.t().print("vecs_mo1_mo2.t()");
          vecs_so_mo1.t().print("vecs_so_mo1.t()");
          so_ao.print("so_ao");
        }

      vecs = vecs_mo1_mo2.t() * vecs_so_mo1.t() * so_ao;
    }
  else {
      if (debug_) ExEnv::out0() << indent << "getting fock matrix" << endl;
      // get the closed shell AO fock matrices
      RefSymmSCMatrix fock_c_so = reference_->fock(0);

      // transform the AO fock matrices to the MO basis
      RefSymmSCMatrix fock_c_mo1
          = basis_matrixkit()->symmmatrix(oso_dimension());
      RefSCMatrix vecs_so_mo1 = reference_->eigenvectors();

      fock_c_mo1.assign(0.0);
      fock_c_mo1.accumulate_transform(vecs_so_mo1.t(), fock_c_so);
      fock_c_so = 0;

      if (debug_) ExEnv::out0() << indent << "diagonalizing" << endl;
      // diagonalize the fock matrix
      vals = fock_c_mo1.eigvals();

      // compute the AO to new MO scf vector
      if (debug_) ExEnv::out0() << indent << "AO to MO" << endl;
      RefSCMatrix so_ao = reference_->integral()->petite_list()->sotoao();
      vecs = vecs_so_mo1.t() * so_ao;
    }
  // fill in the occupations
  occs = matrixkit()->diagmatrix(vals.dim());
  for (i=0; i<oso_dimension()->n(); i++) {
      occs(i) = reference_->occupation(i);
    }
  // allocate storage for symmetry information
  if (!symorb_irrep_) symorb_irrep_ = new int[nbasis];
  if (!symorb_num_) symorb_num_ = new int[nbasis];
  // Check for degenerate eigenvalues.  Use unsorted eigenvalues since it
  // only matters if the degeneracies occur within a given irrep.
  BlockedDiagSCMatrix *bvals = dynamic_cast<BlockedDiagSCMatrix*>(vals.pointer());
  for (i=0; i<bvals->nblocks(); i++) {
      int done = 0;
      RefDiagSCMatrix valsi = bvals->block(i);
      for (j=1; j<valsi.n(); j++) {
          if (fabs(valsi(j)-valsi(j-1)) < 1.0e-7) {
              ExEnv::out0() << indent
                   << "NOTE: There are degenerate orbitals within an irrep."
                   << "  This will make"
                   << endl
                   << indent
                   << "      some diagnostics, such as the largest amplitude,"
                   << " nonunique."
                   << endl;
              done = 1;
              break;
            }
          if (done) break;
        }
    }
  // sort the eigenvectors and values if symmetry is not c1
  if (molecule()->point_group()->char_table().order() != 1) {
      if (debug_) ExEnv::out0() << indent << "sorting eigenvectors" << endl;
      double *evals = new double[noso];
      vals->convert(evals);
      int *indices = new int[noso];
      dquicksort(evals,indices,noso);
      delete[] evals;
      // make sure all nodes see the same indices and evals
      msg_->bcast(indices,noso);
      RefSCMatrix newvecs(vecs.rowdim(), vecs.coldim(), matrixkit());
      RefDiagSCMatrix newvals(vals.dim(), matrixkit());
      RefDiagSCMatrix newoccs(vals.dim(), matrixkit());
      for (i=0; i<noso; i++) {
          newvals(i) = vals(indices[i]);
          newoccs(i) = occs(indices[i]);
          for (j=0; j<nbasis; j++) {
              newvecs(i,j) = vecs(indices[i],j);
            }
        }
      occs = newoccs;
      vecs = newvecs;
      vals = newvals;

      // compute orbital symmetry information
      CharacterTable ct = molecule()->point_group()->char_table();
      int orbnum = 0;
      int *tmp_irrep = new int[noso];
      int *tmp_num = new int[noso];
      for (i=0; i<oso_dimension()->blocks()->nblock(); i++) {
          for (j=0; j<oso_dimension()->blocks()->size(i); j++, orbnum++) {
              tmp_irrep[orbnum] = i;
              tmp_num[orbnum] = j;
            }
        }
      for (i=0; i<noso; i++) {
          symorb_irrep_[i] = tmp_irrep[indices[i]];
          symorb_num_[i] = tmp_num[indices[i]];
        }
      delete[] tmp_irrep;
      delete[] tmp_num;

      delete[] indices;
    }
  else {
      // compute orbital symmetry information for c1
      for (i=0; i<noso; i++) {
          symorb_irrep_[i] = 0;
          symorb_num_[i] = i;
        }
    }
  // check the splitting between frozen and nonfrozen orbitals
  if (nfzc && nfzc < noso) {
      double split = vals(nfzc) - vals(nfzc-1);
      if (split < 0.2) {
          ExEnv::out0() << endl
               << indent << "WARNING: "
               << "MBPT2: gap between frozen and active occupied orbitals is "
               << split << " au" << endl << endl;
        }
    }
  if (nfzv && noso-nfzv-1 >= 0) {
      double split = vals(nbasis-nfzv) - vals(nbasis-nfzv-1);
      if (split < 0.2) {
          ExEnv::out0() << endl
               << indent << "WARNING: "
               << "MBPT2: gap between frozen and active virtual orbitals is "
               << split << " au" << endl << endl;
        }
    }
  if (debug_) ExEnv::out0() << indent << "eigen done" << endl;
}

/////////////////////////////////////////////////////////////////////////////

void
MBPT2::init_variables()
{
  nbasis = so_dimension()->n();
  noso = oso_dimension()->n();
//    if (nbasis != noso) {
//        ExEnv::outn() << "MBPT2: Noso != Nbasis: MBPT2 not checked for this case" << endl;
//        abort();
//      }
  nocc = nvir = nsocc = 0;
  for (int i=0; i<noso; i++) {
    if (reference_->occupation(i) == 2.0) nocc++;
    else if (reference_->occupation(i) == 1.0) nsocc++;
    else nvir++;
    }
}

/////////////////////////////////////////////////////////////////////////////

void
MBPT2::symmetry_changed()
{
  Wavefunction::symmetry_changed();
  reference_->symmetry_changed();
}

/////////////////////////////////////////////////////////////////////////////

int
MBPT2::nelectron()
{
  return reference_->nelectron();
}

/////////////////////////////////////////////////////////////////////////////

double
MBPT2::ref_energy()
{
  return reference_->energy();
}

/////////////////////////////////////////////////////////////////////////////

double
MBPT2::corr_energy()
{
  return energy() - reference_->energy();
}

/////////////////////////////////////////////////////////////////////////////

RefSCVector
MBPT2::ref_energy_gradient()
{
  gradient();
  return hf_gradient_;
}

/////////////////////////////////////////////////////////////////////////////

RefSCVector
MBPT2::corr_energy_gradient()
{
  gradient();
  return get_cartesian_gradient() - hf_gradient_;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:

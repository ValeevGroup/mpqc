//
// scf.cc --- implementation of the SCF abstract base class
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Edward Seidl <seidl@janed.com>
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

#include <math.h>
#include <limits.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include <util/misc/formio.h>
#include <util/state/stateio.h>
#include <util/group/mstate.h>
#include <util/misc/xmlwriter.h>
#include <util/misc/xml.h>

#include <math/scmat/local.h>
#include <math/scmat/repl.h>
#include <math/scmat/offset.h>
#include <math/scmat/blocked.h>

#include <math/optimize/diis.h>

#include <chemistry/qc/basis/petite.h>
#include <chemistry/qc/lcao/soad.h>
#include <chemistry/qc/scf/scf.h>
#include <chemistry/qc/lcao/df_runtime.h>
#ifdef MPQC_NEW_FEATURES
#  include <chemistry/qc/scf/iter_logger.h>
#endif // MPQC_NEW_FEATURES

using namespace std;
using namespace sc;

///////////////////////////////////////////////////////////////////////////
// SCF

static ClassDesc SCF_cd(
  typeid(SCF),"SCF",7,"public OneBodyWavefunction",
  0, 0, 0);

SCF::SCF(StateIn& s) :
  SavableState(s),
  OneBodyWavefunction(s)
{
  compute_guess_ = 0;

  s.get(maxiter_,"maxiter");
  if (s.version(::class_desc<SCF>()) >= 7) {
    s.get(miniter_,"miniter");
  }
  else {
    miniter_ = 0;
  }
  s.get(dens_reset_freq_);
  s.get(reset_occ_);
  s.get(local_dens_);
  if (s.version(::class_desc<SCF>()) >= 3) {
    double dstorage;
    s.get(dstorage);
    storage_ = size_t(dstorage);
  }
  else {
    unsigned int istorage;
    s.get(istorage);
    storage_ = istorage;
  }
  if (s.version(::class_desc<SCF>()) >= 2) {
    s.get(print_all_evals_);
    s.get(print_occ_evals_);
  }
  else {
    print_all_evals_ = 0;
    print_occ_evals_ = 0;
  }
  s.get(level_shift_);
  if (s.version(::class_desc<SCF>()) >= 5) {
    s.get(keep_guess_wfn_);
    guess_wfn_ << SavableState::restore_state(s);
  }
  else keep_guess_wfn_ = 0;
  if (s.version(::class_desc<SCF>()) >= 6) {
    s.get(always_use_guess_wfn_);
  }
  else always_use_guess_wfn_ = 0;

  extrap_ << SavableState::restore_state(s);
  accumdih_ << SavableState::restore_state(s);
  accumddh_ << SavableState::restore_state(s);

  scf_grp_ = basis()->matrixkit()->messagegrp();
  threadgrp_ = ThreadGrp::get_default_threadgrp();
}

SCF::SCF(const Ref<KeyVal>& keyval) :
  OneBodyWavefunction(keyval),
  compute_guess_(0),
  maxiter_(100),
  miniter_(0),
  dens_reset_freq_(10),
  reset_occ_(0),
  local_dens_(1),
  storage_(0),
  level_shift_(0)
{
  if (keyval->exists("maxiter"))
    maxiter_ = keyval->intvalue("maxiter");

  if (keyval->exists("miniter"))
    miniter_ = keyval->intvalue("miniter");

  if (keyval->exists("density_reset_frequency"))
    dens_reset_freq_ = keyval->intvalue("density_reset_frequency");

  if (keyval->exists("reset_occupations"))
    reset_occ_ = keyval->booleanvalue("reset_occupations");

  if (keyval->exists("level_shift"))
    level_shift_ = keyval->doublevalue("level_shift");

  extrap_ << keyval->describedclassvalue("extrap");
  if (extrap_.null())
    extrap_ = new DIIS;

  accumdih_ << keyval->describedclassvalue("accumdih");
  if (accumdih_.null())
    accumdih_ = new AccumHNull;

  accumddh_ << keyval->describedclassvalue("accumddh");
  if (accumddh_.null())
    accumddh_ = new AccumHNull;

  KeyValValuesize defaultmem(DEFAULT_MPQC_MEMORY);
  storage_ = keyval->sizevalue("memory",defaultmem);

  if (keyval->exists("local_density"))
    local_dens_ = keyval->booleanvalue("local_density");

  print_all_evals_ = keyval->booleanvalue("print_evals");
  print_occ_evals_ = keyval->booleanvalue("print_occupied_evals");

  scf_grp_ = basis()->matrixkit()->messagegrp();
  threadgrp_ = ThreadGrp::get_default_threadgrp();

  keep_guess_wfn_ = keyval->booleanvalue("keep_guess_wavefunction");

  always_use_guess_wfn_
    = keyval->booleanvalue("always_use_guess_wavefunction");

  // first see if guess_wavefunction is a wavefunction, then check to
  // see if it's a string.
  if (keyval->exists("guess_wavefunction")) {
    ExEnv::out0() << incindent << incindent;
    guess_wfn_ << keyval->describedclassvalue("guess_wavefunction");
    compute_guess_=1;
    if (guess_wfn_.null()) {
      compute_guess_=0;
      std::string path = keyval->stringvalue("guess_wavefunction");
      struct stat sb;
      if (!path.empty() && stat(path.c_str(), &sb)==0 && sb.st_size) {
        BcastStateInBin s(scf_grp_, path.c_str());

        // reset the default matrixkit so that the matrices in the guess
        // wavefunction will match those in this wavefunction
        Ref<SCMatrixKit> oldkit = SCMatrixKit::default_matrixkit();
        SCMatrixKit::set_default_matrixkit(basis()->matrixkit());

        guess_wfn_ << SavableState::restore_state(s);

        // go back to the original default matrixkit
        SCMatrixKit::set_default_matrixkit(oldkit);
      }
    }
    ExEnv::out0() << decindent << decindent;
  }

  // See if we have an iteration logger
  if(keyval->exists("iter_log")) {
    iter_log_ << keyval->describedclassvalue("iter_log");
  }

}

SCF::~SCF()
{
}

void
SCF::save_data_state(StateOut& s)
{
  OneBodyWavefunction::save_data_state(s);
  s.put(maxiter_);
  s.put(miniter_);
  s.put(dens_reset_freq_);
  s.put(reset_occ_);
  s.put(local_dens_);
  double dstorage = storage_;
  s.put(dstorage);
  s.put(print_all_evals_);
  s.put(print_occ_evals_);
  s.put(level_shift_);
  s.put(keep_guess_wfn_);
  SavableState::save_state(guess_wfn_.pointer(),s);
  s.put(always_use_guess_wfn_);
  SavableState::save_state(extrap_.pointer(),s);
  SavableState::save_state(accumdih_.pointer(),s);
  SavableState::save_state(accumddh_.pointer(),s);
}

#ifdef MPQC_NEW_FEATURES
boost::property_tree::ptree&
SCF::write_xml(
    boost::property_tree::ptree& parent,
    const XMLWriter& writer
)
{
  using boost::property_tree::ptree;
  ptree& my_tree = this->get_my_ptree(parent);
  if(iter_log_.nonnull()){
    writer.insert_child(my_tree, iter_log_, "iteration_log");
  }
  return OneBodyWavefunction::write_xml(parent, writer);
}
#endif // MPQC_NEW_FEATURES

RefSCMatrix
SCF::oso_eigenvectors()
{
  return oso_eigenvectors_.result();
}

RefDiagSCMatrix
SCF::eigenvalues()
{
  return eigenvalues_.result();
}

int
SCF::spin_unrestricted()
{
  return 0;
}

void
SCF::symmetry_changed()
{
  OneBodyWavefunction::symmetry_changed();
  if (guess_wfn_) {
    guess_wfn_->symmetry_changed();
  }
}

void
SCF::print(ostream&o) const
{
  OneBodyWavefunction::print(o);
  o << indent << "SCF Parameters:\n" << incindent
    << indent << "maxiter = " << maxiter_ << endl;
  if (miniter_ > 0) {
    o << indent << "miniter = " << miniter_ << endl;
  }
  o << indent << "density_reset_frequency = " << dens_reset_freq_ << endl
    << indent << scprintf("level_shift = %f\n",level_shift_)
    << decindent << endl;
}

//////////////////////////////////////////////////////////////////////////////

void
SCF::compute()
{
  local_ = (dynamic_cast<LocalSCMatrixKit*>(basis()->matrixkit().pointer()) ||
            dynamic_cast<ReplSCMatrixKit*>(basis()->matrixkit().pointer())) ? 1:0;

  const double hess_to_grad_acc = 1.0/100.0;
  if (hessian_needed())
    set_desired_gradient_accuracy(desired_hessian_accuracy()*hess_to_grad_acc);

  const double grad_to_val_acc = 1.0/100.0;
  if (gradient_needed())
    set_desired_value_accuracy(desired_gradient_accuracy()*grad_to_val_acc);

  double delta;
  if (value_needed()) {
    ExEnv::out0() << endl << indent
         << scprintf("SCF::compute: energy accuracy = %10.7e\n",
                     desired_value_accuracy())
         << endl;

    // calculate the nuclear repulsion energy
    double nucrep = nuclear_repulsion_energy();
    ExEnv::out0() << indent
                  << scprintf("nuclear repulsion energy = %20.15f", nucrep)
                  << endl << endl;

    double eelec;
    delta = compute_vector(eelec,nucrep);

    double eother = 0.0;
    if (accumddh_) eother = accumddh_->e();
    ExEnv::out0() << endl << indent
         << scprintf("total scf energy = %15.10f", eelec+eother+nucrep)
         << endl;

    set_energy(eelec+eother+nucrep);
    set_actual_value_accuracy(delta);
  }
  else {
    delta = actual_value_accuracy();
  }

  if (gradient_needed()) {
    RefSCVector gradient = matrixkit()->vector(moldim());

    ExEnv::out0() << endl << indent
         << scprintf("SCF::compute: gradient accuracy = %10.7e\n",
                     desired_gradient_accuracy())
         << endl;

    compute_gradient(gradient);
    print_natom_3(gradient,"Total Gradient:");
    set_gradient(gradient);

    set_actual_gradient_accuracy(delta/grad_to_val_acc);
  }

  if (hessian_needed()) {
    RefSymmSCMatrix hessian = matrixkit()->symmmatrix(moldim());

    ExEnv::out0() << endl << indent
         << scprintf("SCF::compute: hessian accuracy = %10.7e\n",
                     desired_hessian_accuracy())
         << endl;

    compute_hessian(hessian);
    set_hessian(hessian);

    set_actual_hessian_accuracy(delta/grad_to_val_acc/hess_to_grad_acc);
  }
}

//////////////////////////////////////////////////////////////////////////////

signed char *
SCF::init_pmax(double *pmat_data)
{
  double l2inv = 1.0/log(2.0);
  double tol = pow(2.0,-126.0);

  GaussianBasisSet& gbs = *basis().pointer();

  signed char * pmax = new signed char[i_offset(gbs.nshell())];

  int ish, jsh, ij;
  for (ish=ij=0; ish < gbs.nshell(); ish++) {
    int istart = gbs.shell_to_function(ish);
    int iend = istart + gbs(ish).nfunction();

    for (jsh=0; jsh <= ish; jsh++,ij++) {
      int jstart = gbs.shell_to_function(jsh);
      int jend = jstart + gbs(jsh).nfunction();

      double maxp=0, tmp;

      for (int i=istart; i < iend; i++) {
        int ijoff = i_offset(i) + jstart;
        for (int j=jstart; j < ((ish==jsh) ? i+1 : jend); j++,ijoff++)
          if ((tmp=fabs(pmat_data[ijoff])) > maxp)
            maxp=tmp;
      }

      if (maxp <= tol)
        maxp=tol;

      long power = long(ceil(log(maxp)*l2inv));
      if (power < SCHAR_MIN) pmax[ij] = SCHAR_MIN;
      else if (power > SCHAR_MAX) pmax[ij] = SCHAR_MAX;
      else pmax[ij] = (signed char) power;
    }
  }

  return pmax;
}

//////////////////////////////////////////////////////////////////////////////

RefSymmSCMatrix
SCF::get_local_data(const RefSymmSCMatrix& m, double*& p, Access access)
{
  RefSymmSCMatrix l = m;

  if (!dynamic_cast<LocalSymmSCMatrix*>(l.pointer())
      && !dynamic_cast<ReplSymmSCMatrix*>(l.pointer())) {
    Ref<SCMatrixKit> k = new ReplSCMatrixKit;
    l = k->symmmatrix(m.dim());
    l->convert(m);

    if (access == Accum)
      l->assign(0.0);
  } else if (scf_grp_->n() > 1 && access==Accum) {
    l = m.clone();
    l.assign(0.0);
  }

  if (dynamic_cast<ReplSymmSCMatrix*>(l.pointer()))
    p = dynamic_cast<ReplSymmSCMatrix*>(l.pointer())->get_data();
  else
    p = dynamic_cast<LocalSymmSCMatrix*>(l.pointer())->get_data();

  return l;
}

//////////////////////////////////////////////////////////////////////////////

void
SCF::initial_vector()
{
  const bool vector_is_null = oso_eigenvectors_.result_noupdate().null();
  if (vector_is_null) {
    bool soad_guess = false;
    if (guess_wfn_.null()) {
      Ref<AssignedKeyVal> akv = new AssignedKeyVal;
      akv->assign("molecule", molecule().pointer());

      // small polarized basis should be enough for guess
      // try it if at least 2 times smaller than the full basis
      // use the full basis as backup
      try {
        akv->assign("name", "Def2-SV(P)");
        Ref<GaussianBasisSet> guess_bs = new GaussianBasisSet(akv);
        if (guess_bs->nbasis() < basis()->nbasis()/2)
          akv->assign("basis", guess_bs.pointer());
        else
          akv->assign("basis", basis().pointer());
      }
      catch(...) {
        akv->assign("basis", basis().pointer());
      }

      try {
        guess_wfn_ = new SuperpositionOfAtomicDensities(Ref<KeyVal>(akv));
        //guess_wfn_ = new HCoreWfn(Ref<KeyVal>(akv));   not done because the "diagonalize-in-SO-basis" is not implemented
        soad_guess = true;
      }
      catch (...) { // failed? Oh well, will resort to core guess.
      }
    }

    // if guess_wfn_ is non-null then try to get a guess vector from it.
    // First check that the same basis is used...if not, then project the
    // guess vector into the present basis.
    if (guess_wfn_) {

      // compute guess wfn with lower accuracy than this wfn
      if (guess_wfn_->desired_value_accuracy_set_to_default())
        guess_wfn_->set_desired_value_accuracy( this->desired_value_accuracy() * this->guess_acc_ratio() );

      if (basis()->equiv(guess_wfn_->basis())
          &&orthog_method() == guess_wfn_->orthog_method()
          &&oso_dimension()->equiv(guess_wfn_->oso_dimension().pointer())) {
        ExEnv::out0() << indent
                      << "Using " << guess_wfn_->class_name()
                      << " guess wavefunction as starting vector" << endl;

        // indent output of eigenvectors() call if there is any
        ExEnv::out0() << incindent;
        SCF *g = dynamic_cast<SCF*>(guess_wfn_.pointer());
        if (!g || compute_guess_) {
          oso_eigenvectors_ = guess_wfn_->oso_eigenvectors().copy();
          eigenvalues_ = guess_wfn_->eigenvalues().copy();
          current_evals_ = eigenvalues_.result_noupdate();
        } else {
          oso_eigenvectors_ = g->oso_eigenvectors_.result_noupdate().copy();
          eigenvalues_ = g->eigenvalues_.result_noupdate().copy();
          current_evals_ = eigenvalues_.result_noupdate();
        }
        ExEnv::out0() << decindent;
      } else {
        ExEnv::out0() << indent
                      << "Projecting "
                      << guess_wfn_->class_name()
                      << " guess into the present basis set"
                      << endl;

        // indent output of projected_eigenvectors() call if there is any
        ExEnv::out0() << incindent << incindent;
        oso_eigenvectors_ = projected_eigenvectors(guess_wfn_);
        eigenvalues_ = projected_eigenvalues(guess_wfn_);
        current_evals_ = eigenvalues_.result_noupdate();
        ExEnv::out0() << decindent << decindent;
      }
      // we should only have to do this once, so free up memory used
      // for the old wavefunction, unless told otherwise
      if (!keep_guess_wfn_) guess_wfn_=0;
      // if made a SOAD guess, erase it
      if (soad_guess) guess_wfn_=0;

      ExEnv::out0() << endl;

    } else { // if all else failed
      ExEnv::out0() << indent << "Starting from core Hamiltonian guess\n"
                    << endl;
      oso_eigenvectors_ = hcore_guess(eigenvalues_.result_noupdate());
      current_evals_ = eigenvalues_.result_noupdate();
    }
  } else {
    // if vector exists do nothing
  }
} // end of initial_vector()

//////////////////////////////////////////////////////////////////////////////

void
SCF::init_mem(int nm)
{
  // if local_den_ is already 0, then that means it was set to zero by
  // the user.
  if (!local_dens_) {
    integral()->set_storage(storage_);
    return;
  }

  size_t nmem = i_offset(basis()->nbasis())*nm*sizeof(double);

  // if we're actually using local matrices, then there's no choice
  if (dynamic_cast<LocalSCMatrixKit*>(basis()->matrixkit().pointer())
      ||dynamic_cast<ReplSCMatrixKit*>(basis()->matrixkit().pointer())) {
    if (nmem > storage_)
      return;
  } else {
    if (nmem > storage_) {
      local_dens_=0;
      integral()->set_storage(storage_);
      return;
    }
  }

  integral()->set_storage(storage_-nmem);
}

/////////////////////////////////////////////////////////////////////////////

void
SCF::so_density(const RefSymmSCMatrix& d, double occ, int alp)
{
  int i,j,k;
  int me=scf_grp_->me();
  int nproc=scf_grp_->n();
  int uhf = spin_unrestricted();

  d->assign(0.0);

  RefSCMatrix oso_vector;
  if (alp || !uhf) {
    if (oso_scf_vector_)
      oso_vector = oso_scf_vector_;
  }
  else {
    if (oso_scf_vector_beta_)
      oso_vector = oso_scf_vector_beta_;
  }

  if (oso_vector.null()) {
    if (uhf) {
      if (alp)
        oso_vector = oso_alpha_eigenvectors();
      else
        oso_vector = oso_beta_eigenvectors();
    } else
      oso_vector = oso_eigenvectors();
  }

  if (debug_ > 1) oso_vector.print("ortho SO vector");

  RefSCMatrix vector = so_to_orthog_so().t() * oso_vector;
  oso_vector = 0;

  if (debug_ > 1) vector.print("SO vector");

  BlockedSCMatrix *bvec = require_dynamic_cast<BlockedSCMatrix*>(
    vector, "SCF::so_density: blocked vector");

  BlockedSymmSCMatrix *bd = require_dynamic_cast<BlockedSymmSCMatrix*>(
    d, "SCF::so_density: blocked density");

  for (int ir=0; ir < oso_dimension()->blocks()->nblock(); ir++) {
    RefSCMatrix vir = bvec->block(ir);
    RefSymmSCMatrix dir = bd->block(ir);

    if (vir.null() || vir.ncol()==0)
      continue;

    int n_orthoSO = oso_dimension()->blocks()->size(ir);
    int n_SO = so_dimension()->blocks()->size(ir);

    // figure out which columns of the scf vector we'll need
    int col0 = -1, coln = -1;
    for (i=0; i < n_orthoSO; i++) {
      double occi;
      if (!uhf)
        occi = occupation(ir, i);
      else if (alp)
        occi = alpha_occupation(ir, i);
      else
        occi = beta_occupation(ir, i);

      if (fabs(occi-occ) < 1.0e-8) {
        if (col0 == -1)
          col0 = i;
        continue;
      } else if (col0 != -1) {
        coln = i-1;
        break;
      }
    }

    if (col0 == -1)
      continue;

    if (coln == -1)
      coln = n_orthoSO-1;

    if (local_ || local_dens_) {
      RefSymmSCMatrix ldir = dir;

      RefSCMatrix occbits; // holds the occupied bits of the scf vector

      // get local copies of vector and density matrix
      if (!local_) {
        Ref<SCMatrixKit> rk = new ReplSCMatrixKit;
        RefSCMatrix lvir = rk->matrix(vir.rowdim(), vir.coldim());
        lvir->convert(vir);
        occbits = lvir->get_subblock(0, n_SO-1, col0, coln);
        lvir = 0;

        ldir = rk->symmmatrix(dir.dim());
        ldir->convert(dir);

      } else {
        occbits = vir->get_subblock(0, n_SO-1, col0, coln);
      }

      double **c;
      double *dens;

      if (dynamic_cast<LocalSCMatrix*>(occbits.pointer()))
        c = dynamic_cast<LocalSCMatrix*>(occbits.pointer())->get_rows();
      else if (dynamic_cast<ReplSCMatrix*>(occbits.pointer()))
        c = dynamic_cast<ReplSCMatrix*>(occbits.pointer())->get_rows();
      else
        abort();

      if (dynamic_cast<LocalSymmSCMatrix*>(ldir.pointer()))
        dens = dynamic_cast<LocalSymmSCMatrix*>(ldir.pointer())->get_data();
      else if (dynamic_cast<ReplSymmSCMatrix*>(ldir.pointer()))
        dens = dynamic_cast<ReplSymmSCMatrix*>(ldir.pointer())->get_data();
      else
        abort();

      int ij=0;
      for (i=0; i < n_SO; i++) {
        for (j=0; j <= i; j++, ij++) {
          if (ij%nproc != me)
            continue;

          double dv = 0;

          int kk=0;
          for (k=col0; k <= coln; k++, kk++)
            dv += c[i][kk]*c[j][kk];

          dens[ij] = dv;
        }
      }

      if (nproc > 1)
        scf_grp_->sum(dens, i_offset(n_SO));

      if (!local_) {
        dir->convert(ldir);
      }
    }

    // for now quit
    else {
      ExEnv::err0() << indent
           << "Cannot yet use anything but Local matrices"
           << endl;
      abort();
    }
  }

  if (debug_ > 0) {
    ExEnv::out0() << indent
         << "Nelectron = " << 2.0 * (d * overlap()).trace() << endl;
  }

  if (debug_ > 1) {
    d.print("SO Density");
    RefSCMatrix rd(d.dim(), d.dim(), basis_matrixkit());
    rd.assign(0.0);
    rd.accumulate(d);
    (d*overlap()*d-rd).print("SO Density idempotency error");
  }
}

double
SCF::one_body_energy()
{
  RefSymmSCMatrix dens = ao_density().copy();
  RefSymmSCMatrix hcore = dens->clone();
  hcore.assign(0.0);
  Ref<SCElementOp> hcore_op = new OneBodyIntOp(integral()->hcore());
  hcore.element_op(hcore_op);

  dens->scale_diagonal(0.5);
  SCElementScalarProduct *prod = new SCElementScalarProduct;
  prod->reference();
  Ref<SCElementOp2> op = prod;
  hcore->element_op(prod, dens);
  double e = prod->result();
  op = 0;
  prod->dereference();
  delete prod;
  return 2.0 * e;
}

void
SCF::two_body_energy(double &ec, double &ex)
{
  ExEnv::errn() << class_name() << ": two_body_energy not implemented" << endl;
}

/////////////////////////////////////////////////////////////////////////////

void
SCF::init_threads()
{
  int nthread = threadgrp_->nthread();
  size_t int_store = integral()->storage_unused()/nthread;

  // initialize the two electron integral classes
  tbis_ = new Ref<TwoBodyInt>[nthread];
  for (int i=0; i < nthread; i++) {
    tbis_[i] = integral()->electron_repulsion();
    tbis_[i]->set_integral_storage(int_store);
  }

}

void
SCF::done_threads()
{
  for (int i=0; i < threadgrp_->nthread(); i++) tbis_[i] = 0;
  delete[] tbis_;
  tbis_ = 0;
}

int *
SCF::read_occ(const Ref<KeyVal> &keyval, const char *name, int nirrep)
{
  int *occ = 0;
  if (keyval->exists(name)) {
    if (keyval->count(name) != nirrep) {
      ExEnv::err0() << indent
                   << "ERROR: SCF: have " << nirrep << " irreps but "
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

void
SCF::obsolete()
{
  OneBodyWavefunction::obsolete();
  if (guess_wfn_) guess_wfn_->obsolete();
  // do I need to obsolete the vector also here? Yes, if always_use_guess_wfn_ is set to true.
  // Otherwise, the user of this class knows the context of the call to be able to call purge(), e.g.
  // in geometry optimization vector may be reused, but in set_orthog_method it currently can't
  if (always_use_guess_wfn_) purge();
}

void
SCF::obsolete_vector() {
  oso_eigenvectors_ = RefSCMatrix(0);
}

void
SCF::purge() {
  Wavefunction::purge();
  obsolete_vector();
}

Ref<SCExtrapData>
SCF::initial_extrap_data()
{
  return 0;
}

Ref<DensityFittingInfo>
SCF::dfinfo() const {
  return 0;
}

#ifdef MPQC_NEW_FEATURES
///////////////////////////////////////////////////////////////////////////
// SCFIterationLogger

static ClassDesc SCFIterationLogger_cd(
  typeid(SCFIterationLogger),"SCFIterationLogger",1,"public XMLWritable, public DescribedClass",
  0, create<SCFIterationLogger>, 0);

SCFIterationLogger::SCFIterationLogger(const Ref<KeyVal>& keyval) :
    log_evals_(keyval->booleanvalue("log_evals", KeyValValueboolean(false))),
    log_density_(keyval->booleanvalue("log_density", KeyValValueboolean(false))),
    log_coeffs_(keyval->booleanvalue("log_coefficients", KeyValValueboolean(false)))
{

}

void
SCFIterationLogger::new_iteration(){
  SCFIterationData iteration;
  iteration.parent = this;
  iteration.number = iterations_.size() + 1;
  iterations_.push_back(iteration);
  other_iter_details_.emplace_back(0);

}

void
SCFIterationLogger::log_density(
    RefSymmSCMatrix density,
    SpinCase1 spin
) {
  if(not log_density_) return;

  switch (spin) {
  case AnySpinCase1:
    iterations_.back().density = density;
    break;
  case Alpha:
    iterations_.back().alpha_density = density;
    break;
  case Beta:
    iterations_.back().beta_density = density;
    break;
  case InvalidSpinCase1:
    throw;
  }
}

void
SCFIterationLogger::log_evals(
    RefDiagSCMatrix evals,
    SpinCase1 spin
) {
  if(not log_evals_) return;

  switch (spin) {
  case AnySpinCase1:
    iterations_.back().evals = evals;
    break;
  case Alpha:
    iterations_.back().alpha_evals = evals;
    break;
  case Beta:
    iterations_.back().beta_evals = evals;
    break;
  case InvalidSpinCase1:
    throw;
  }
}

void
SCFIterationLogger::log_coeffs(
    RefSCMatrix coeffs,
    SpinCase1 spin
) {
  if(not log_evals_) return;

  switch (spin) {
  case AnySpinCase1:
    iterations_.back().coeffs = coeffs;
    break;
  case Alpha:
    iterations_.back().alpha_coeffs = coeffs;
    break;
  case Beta:
    iterations_.back().beta_coeffs = coeffs;
    break;
  case InvalidSpinCase1:
    throw;
  }
}

#ifdef MPQC_NEW_FEATURES
boost::property_tree::ptree&
SCFIterationLogger::write_xml(
    boost::property_tree::ptree& parent,
    const XMLWriter& writer
)
{
  using boost::property_tree::ptree;
  ptree* tmp;
  if(writer.fold_in_class_name()){
    tmp = &parent;
    tmp->put("<xmlattr>.typename", "SCFIterationLogger");
  }
  else{
    ptree& tmpref = parent.add_child("SCFIterationLogger", ptree());
    tmp = &tmpref;

  }
  ptree& child = *tmp;
  child.put("density_enabled", log_density_);
  child.put("coefficients_enabled", log_coeffs_);
  child.put("evals_enabled", log_evals_);
  for(auto&& func : other_details_) {
    func(child, writer);
  }

  // Write the iterations
  ptree& iter_tree = child.add_child("iterations", ptree());
  iter_tree.put("<xmlattr>.n", iterations_.size());
  int itnum = 0;
  for(auto& itdata : iterations_){
    ptree& ichild = writer.insert_child(iter_tree, itdata, "iteration");
    for(auto&& func : other_iter_details_[itnum]) {
      func(ichild, writer);
    }
    ++itnum;
  }
  return child;
}
#endif // MPQC_NEW_FEATURES

///////////////////////////////////////////////////////////////////////////
// SCFIterationData

ptree&
SCFIterationData::write_xml(
    boost::property_tree::ptree& parent,
    const XMLWriter& writer
)
{
  ptree* tmp;
  if(writer.fold_in_class_name()){
    tmp = &parent;
    tmp->put("<xmlattr>.typename", "SCFIterationData");
  }
  else{
    ptree& tmpref = parent.add_child("SCFIterationData", ptree());
    tmp = &tmpref;
  }
  ptree& my_tree = *tmp;

  my_tree.put("<xmlattr>.number", number);

  if(evals.nonnull()){
    writer.insert_child(my_tree, evals, "evals");
  }
  if(alpha_evals.nonnull()){
    ptree& child = writer.insert_child(my_tree, alpha_evals, "evals");
    child.put("<xmlattr>.spin", "alpha");
  }
  if(beta_evals.nonnull()){
    ptree& child = writer.insert_child(my_tree, beta_evals, "evals");
    child.put("<xmlattr>.spin", "beta");
  }
  //----------------------------------------//
  if(density.nonnull()){
    ptree& child = writer.insert_child(my_tree, density, "density");
  }
  if(alpha_density.nonnull()){
    ptree& child = writer.insert_child(my_tree, alpha_density, "density");
    child.put("<xmlattr>.spin", "alpha");
  }
  if(beta_density.nonnull()){
    ptree& child = writer.insert_child(my_tree, beta_density, "density");
    child.put("<xmlattr>.spin", "beta");
  }
  //----------------------------------------//
  if(coeffs.nonnull()){
    ptree& child = writer.insert_child(my_tree, coeffs, "coefficients");
  }
  if(alpha_coeffs.nonnull()){
    ptree& child = writer.insert_child(my_tree, alpha_coeffs, "coefficients");
    child.put("<xmlattr>.spin", "alpha");
  }
  if(beta_coeffs.nonnull()){
    ptree& child = writer.insert_child(my_tree, beta_coeffs, "coefficients");
    child.put("<xmlattr>.spin", "beta");
  }

  return my_tree;
}
#endif // MPQC_NEW_FEATURES

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:

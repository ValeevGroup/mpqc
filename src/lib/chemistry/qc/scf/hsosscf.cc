//
// hsosscf.cc --- implementation of the high-spin open shell SCF class
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
#include <numeric>

#include <util/misc/regtime.h>
#include <util/misc/formio.h>
#include <util/state/stateio.h>
#include <util/misc/scexception.h>

#include <math/scmat/block.h>
#include <math/scmat/blocked.h>
#include <math/scmat/blkiter.h>
#include <math/scmat/local.h>

#include <math/optimize/scextrapmat.h>

#include <chemistry/qc/basis/petite.h>

#include <chemistry/qc/scf/scflocal.h>
#include <chemistry/qc/scf/scfops.h>
#include <chemistry/qc/scf/effh.h>
#include <chemistry/qc/scf/hsosscf.h>

#include <chemistry/qc/scf/ltbgrad.h>
#include <chemistry/qc/scf/hsoshftmpl.h>
#include <chemistry/qc/wfn/femo.h>

using namespace std;
using namespace sc;

///////////////////////////////////////////////////////////////////////////
// HSOSSCF

static ClassDesc HSOSSCF_cd(
  typeid(HSOSSCF),"HSOSSCF",3,"public SCF",
  0, 0, 0);

HSOSSCF::HSOSSCF(StateIn& s) :
  SavableState(s),
  SCF(s),
  cl_fock_(this),
  op_fock_(this),
  alpha_semican_evals_(this),
  beta_semican_evals_(this),
  alpha_semican_evecs_(this),
  beta_semican_evecs_(this)
{
  cl_fock_.result_noupdate() =
    basis_matrixkit()->symmmatrix(so_dimension());
  cl_fock_.restore_state(s);
  cl_fock_.result_noupdate().restore(s);

  op_fock_.result_noupdate() =
    basis_matrixkit()->symmmatrix(so_dimension());
  op_fock_.restore_state(s);
  op_fock_.result_noupdate().restore(s);

  s.get(user_occupations_);
  s.get(tndocc_);
  s.get(tnsocc_);
  s.get(nirrep_);
  s.get(ndocc_);
  s.get(nsocc_);
  if (s.version(::class_desc<HSOSSCF>()) >= 2) {
    s.get(initial_ndocc_);
    s.get(initial_nsocc_);
    most_recent_pg_ << SavableState::restore_state(s);
  } else {
    initial_ndocc_ = new int[nirrep_];
    memcpy(initial_ndocc_, ndocc_, sizeof(int)*nirrep_);
    initial_nsocc_ = new int[nirrep_];
    memcpy(initial_nsocc_, nsocc_, sizeof(int)*nirrep_);
  }

  if (s.version(::class_desc<HSOSSCF>()) >= 3) {
    int semican_eval; s.get(semican_eval); semicanonical_evaluated_ = (bool)semican_eval;
    if (semicanonical_evaluated_) {
      alpha_semican_evals_.result_noupdate() =
        basis_matrixkit()->diagmatrix(oso_dimension());
      alpha_semican_evals_.restore_state(s);
      alpha_semican_evals_.result_noupdate().restore(s);

      beta_semican_evals_.result_noupdate() =
        basis_matrixkit()->diagmatrix(oso_dimension());
      beta_semican_evals_.restore_state(s);
      beta_semican_evals_.result_noupdate().restore(s);

      alpha_semican_evecs_.result_noupdate() =
        basis_matrixkit()->matrix(so_dimension(),oso_dimension());
      alpha_semican_evecs_.restore_state(s);
      alpha_semican_evecs_.result_noupdate().restore(s);

      beta_semican_evecs_.result_noupdate() =
        basis_matrixkit()->matrix(so_dimension(),oso_dimension());
      beta_semican_evecs_.restore_state(s);
      beta_semican_evecs_.result_noupdate().restore(s);
    }
  }

  // now take care of memory stuff
  init_mem(4);
}

HSOSSCF::HSOSSCF(const Ref<KeyVal>& keyval) :
  SCF(keyval),
  cl_fock_(this),
  op_fock_(this),
  alpha_semican_evals_(this),
  beta_semican_evals_(this),
  alpha_semican_evecs_(this),
  beta_semican_evecs_(this),
  semicanonical_evaluated_(false)
{
  int i;

  cl_fock_.compute()=0;
  cl_fock_.computed()=0;

  op_fock_.compute()=0;
  op_fock_.computed()=0;

  alpha_semican_evals_.compute()=0;
  alpha_semican_evals_.computed()=0;
  beta_semican_evals_.compute()=0;
  beta_semican_evals_.computed()=0;
  alpha_semican_evecs_.compute()=0;
  alpha_semican_evecs_.computed()=0;
  beta_semican_evecs_.compute()=0;
  beta_semican_evecs_.computed()=0;

  // calculate the total nuclear charge
  const int Znuc=molecule()->total_Z();

  // check to see if this is to be a charged molecule
  const int charge = keyval->intvalue("total_charge", KeyValValueint(0));
  const int nelectrons = Znuc-charge;

  bool magnetic_moment_given = false;
  // first let's try to figure out how many open shells there are
  if (keyval->exists("nsocc")) {
    tnsocc_ = keyval->intvalue("nsocc");
  } else if (keyval->exists("magnetic_moment")) {
    tnsocc_ = keyval->intvalue("magnetic_moment");
    magnetic_moment_given = true;
  } else if (keyval->exists("multiplicity")) { // obsolete
    tnsocc_ = keyval->intvalue("multiplicity")-1;
    magnetic_moment_given = true;
  } else {
    // if there's an odd number of electrons, then do a doublet, otherwise
    // do a triplet
    if (nelectrons%2)
      tnsocc_=1;
    else
      tnsocc_=2;
  }

  // now do the same for the number of doubly occupied shells
  if (keyval->exists("ndocc")) {
    tndocc_ = keyval->intvalue("ndocc");
  } else {
    tndocc_ = (nelectrons-tnsocc_)/2;
    if ((nelectrons-tnsocc_)%2) {
      ExEnv::err0() << endl << indent
           << "HSOSSCF::init: Warning, there's a leftover electron.\n"
           << incindent << indent << "total_charge = " << charge << endl
           << indent << "total nuclear charge = " << Znuc << endl
           << indent << "ndocc_ = " << tndocc_ << endl
           << indent << "nsocc_ = " << tnsocc_ << endl << decindent;
    }
  }

  ExEnv::out0() << endl << indent
       << "HSOSSCF::init: total charge = " << Znuc-2*tndocc_-tnsocc_
       << endl << endl;

  nirrep_ = molecule()->point_group()->char_table().ncomp();

  ndocc_ = read_occ(keyval, "docc", nirrep_);
  nsocc_ = read_occ(keyval, "socc", nirrep_);
  if (ndocc_ && nsocc_) {
    user_occupations_=1;
    initial_ndocc_ = new int[nirrep_];
    memcpy(initial_ndocc_, ndocc_, sizeof(int)*nirrep_);
    initial_nsocc_ = new int[nirrep_];
    memcpy(initial_nsocc_, nsocc_, sizeof(int)*nirrep_);
  }
  else if ((ndocc_ && (!nsocc_)) || ((!ndocc_) && nsocc_)) {
    ExEnv::outn() << "ERROR: HSOSSCF: only one of docc and socc specified: "
                 << "give both or none" << endl;
    abort();
  }
  else {
    ndocc_=0;
    nsocc_=0;
    initial_ndocc_=0;
    initial_nsocc_=0;
    user_occupations_=0;
    // second argument: allow to change magnetic moment?
    set_occupations(0, magnetic_moment_given? false : true);
  }

  ExEnv::out0() << indent << "docc = [";
  for (i=0; i < nirrep_; i++)
    ExEnv::out0() << " " << ndocc_[i];
  ExEnv::out0() << " ]\n";

  ExEnv::out0() << indent << "socc = [";
  for (i=0; i < nirrep_; i++)
    ExEnv::out0() << " " << nsocc_[i];
  ExEnv::out0() << " ]\n";

  // check to see if this was done in SCF(keyval)
  if (!keyval->exists("maxiter"))
    maxiter_ = 100;

  if (!keyval->exists("level_shift"))
    level_shift_ = 0.25;

  // now take care of memory stuff
  init_mem(4);
}

HSOSSCF::~HSOSSCF()
{
  if (ndocc_) {
    delete[] ndocc_;
    ndocc_=0;
  }
  if (nsocc_) {
    delete[] nsocc_;
    nsocc_=0;
  }
  delete[] initial_ndocc_;
  delete[] initial_nsocc_;
}

void
HSOSSCF::save_data_state(StateOut& s)
{
  SCF::save_data_state(s);

  cl_fock_.save_data_state(s);
  cl_fock_.result_noupdate().save(s);

  op_fock_.save_data_state(s);
  op_fock_.result_noupdate().save(s);

  s.put(user_occupations_);
  s.put(tndocc_);
  s.put(tnsocc_);
  s.put(nirrep_);
  s.put(ndocc_,nirrep_);
  s.put(nsocc_,nirrep_);
  s.put(initial_ndocc_,nirrep_);
  s.put(initial_nsocc_,nirrep_);
  SavableState::save_state(most_recent_pg_.pointer(),s);

  s.put((int)semicanonical_evaluated_);
  if (semicanonical_evaluated_) {
    alpha_semican_evals_.save_data_state(s);
    alpha_semican_evals_.result_noupdate().save(s);
    beta_semican_evals_.save_data_state(s);
    beta_semican_evals_.result_noupdate().save(s);
    alpha_semican_evecs_.save_data_state(s);
    alpha_semican_evecs_.result_noupdate().save(s);
    beta_semican_evecs_.save_data_state(s);
    beta_semican_evecs_.result_noupdate().save(s);
  }
}

double
HSOSSCF::occupation(int ir, int i)
{
  if (i < ndocc_[ir]) return 2.0;
  else if (i < ndocc_[ir] + nsocc_[ir]) return 1.0;
  return 0.0;
}

double
HSOSSCF::alpha_occupation(int ir, int i)
{
  if (i < ndocc_[ir] + nsocc_[ir]) return 1.0;
  return 0.0;
}

double
HSOSSCF::beta_occupation(int ir, int i)
{
  if (i < ndocc_[ir]) return 1.0;
  return 0.0;
}

int
HSOSSCF::n_fock_matrices() const
{
  return 2;
}

RefSymmSCMatrix
HSOSSCF::fock(int n)
{
  if (n > 1) {
    ExEnv::err0() << indent
         << "HSOSSCF::fock: there are only two fock matrices, "
         << scprintf("but fock(%d) was requested\n",n);
    abort();
  }

  if (n==0)
    return cl_fock_.result();
  else
    return op_fock_.result();
}

double
HSOSSCF::magnetic_moment() const
{
  const int mm = std::accumulate(nsocc_, nsocc_+nirrep_, 0);
  return static_cast<double>(mm);
}

void
HSOSSCF::print(ostream&o) const
{
  int i;

  SCF::print(o);
  o << indent << "HSOSSCF Parameters:\n" << incindent
    << indent << "charge = " << molecule()->total_charge()
                                - 2*tndocc_ - tnsocc_ << endl
    << indent << "ndocc = " << tndocc_ << endl
    << indent << "nsocc = " << tnsocc_ << endl
    << indent << "docc = [";
  for (i=0; i < nirrep_; i++)
    o << " " << ndocc_[i];
  o << " ]" << endl;

  o << indent << "socc = [";
  for (i=0; i < nirrep_; i++)
    o << " " << nsocc_[i];
  o << " ]" << endl << decindent << endl;
}

//////////////////////////////////////////////////////////////////////////////

void
HSOSSCF::set_occupations(const RefDiagSCMatrix& ev)
{
  // set_occupations preserves magnetic_moment
  set_occupations(ev,false);
}

void
HSOSSCF::set_occupations(const RefDiagSCMatrix& ev, bool can_change_magnetic_moment)
{
  if (user_occupations_ || (initial_ndocc_ && initial_nsocc_ && ev.null())) {
    if (form_occupations(ndocc_, initial_ndocc_)
        &&form_occupations(nsocc_, initial_nsocc_)) {
      most_recent_pg_ = new PointGroup(molecule()->point_group());
      return;
    }
    delete[] ndocc_; ndocc_ = 0;
    delete[] nsocc_; nsocc_ = 0;
    ExEnv::out0() << indent
         << "HSOSSCF: WARNING: reforming occupation vectors from scratch"
         << endl;
  }

  int i,j;

  RefDiagSCMatrix evals;

  if (ev.null()) {
    initial_vector();
    evals = eigenvalues_.result_noupdate();
  }
  else
    evals = ev;

  //
  // populate orbitals
  //
  // if can change the magnetic moment, use HundsFEMOSeeker to get the FEMO with maximum magnetic moment
  Ref<FEMO> femo;
  if (can_change_magnetic_moment) {
    //                            # electron,       Etol, allow closed-shell?
    HundsFEMOSeeker femoseeker(tndocc_*2 + tnsocc_, HundsFEMOSeeker::tolerance, false,
                               evals);
    femo = femoseeker.result();
  }
  else {
    femo = new FEMO(tndocc_ + tnsocc_, tndocc_,evals);
  }
  // copy into local arrays
  int *newdocc = new int[nirrep_];
  int *newsocc = new int[nirrep_];
  tndocc_ = tnsocc_ = 0;
  for(int g=0; g<nirrep_; ++g) {
    const int na = femo->nalpha(g);
    const int nb = femo->nbeta(g);
    newdocc[g] = nb;
    newsocc[g] = na - nb;
    tndocc_ += nb;
    tnsocc_ += (na - nb);
  }

  if (!ndocc_) {
    ndocc_=newdocc;
    nsocc_=newsocc;
  } else if (most_recent_pg_.nonnull()
             && most_recent_pg_->equiv(molecule()->point_group())) {
    // test to see if newocc is different from ndocc_
    for (i=0; i < nirrep_; i++) {
      if (ndocc_[i] != newdocc[i]) {
        ExEnv::err0() << indent << "HSOSSCF::set_occupations:  WARNING!!!!\n"
             << incindent << indent
             << scprintf("occupations for irrep %d have changed\n",i+1)
             << indent
             << scprintf("ndocc was %d, changed to %d", ndocc_[i], newdocc[i])
             << endl << decindent;
      }
      if (nsocc_[i] != newsocc[i]) {
        ExEnv::err0() << indent << "HSOSSCF::set_occupations:  WARNING!!!!\n"
             << incindent << indent
             << scprintf("occupations for irrep %d have changed\n",i+1)
             << indent
             << scprintf("nsocc was %d, changed to %d", nsocc_[i], newsocc[i])
             << endl << decindent;
      }
    }

    memcpy(ndocc_,newdocc,sizeof(int)*nirrep_);
    memcpy(nsocc_,newsocc,sizeof(int)*nirrep_);
    delete[] newdocc;
    delete[] newsocc;
  }

  if (!initial_ndocc_
      || initial_pg_->equiv(molecule()->point_group())) {
    delete[] initial_ndocc_;
    initial_ndocc_ = new int[nirrep_];
    memcpy(initial_ndocc_,ndocc_,sizeof(int)*nirrep_);
  }

  if (!initial_nsocc_
      || initial_pg_->equiv(molecule()->point_group())) {
    delete[] initial_nsocc_;
    initial_nsocc_ = new int[nirrep_];
    memcpy(initial_nsocc_,nsocc_,sizeof(int)*nirrep_);
  }

  most_recent_pg_ = new PointGroup(molecule()->point_group());
}

void
HSOSSCF::symmetry_changed()
{
  SCF::symmetry_changed();
  cl_fock_.result_noupdate()=0;
  op_fock_.result_noupdate()=0;
  alpha_semican_evals_.result_noupdate()=0;
  beta_semican_evals_.result_noupdate()=0;
  alpha_semican_evecs_.result_noupdate()=0;
  beta_semican_evecs_.result_noupdate()=0;
  nirrep_ = molecule()->point_group()->char_table().ncomp();
  set_occupations(0);
}

//////////////////////////////////////////////////////////////////////////////
//
// scf things
//

void
HSOSSCF::init_vector()
{
  init_threads();

  // allocate storage for other temp matrices
  cl_dens_ = hcore_.clone();
  cl_dens_.assign(0.0);

  cl_dens_diff_ = hcore_.clone();
  cl_dens_diff_.assign(0.0);

  op_dens_ = hcore_.clone();
  op_dens_.assign(0.0);

  op_dens_diff_ = hcore_.clone();
  op_dens_diff_.assign(0.0);

  // gmat is in AO basis
  cl_gmat_ = basis()->matrixkit()->symmmatrix(basis()->basisdim());
  cl_gmat_.assign(0.0);

  op_gmat_ = cl_gmat_.clone();
  op_gmat_.assign(0.0);

  if (cl_fock_.result_noupdate().null()) {
    cl_fock_ = hcore_.clone();
    cl_fock_.result_noupdate().assign(0.0);
    op_fock_ = hcore_.clone();
    op_fock_.result_noupdate().assign(0.0);
  }

  // make sure trial vector is set up
  initial_vector();

  oso_scf_vector_ = oso_eigenvectors_.result_noupdate();
}

void
HSOSSCF::done_vector()
{
  done_threads();

  cl_gmat_ = 0;
  cl_dens_ = 0;
  cl_dens_diff_ = 0;
  op_gmat_ = 0;
  op_dens_ = 0;
  op_dens_diff_ = 0;

  oso_scf_vector_ = 0;
}

RefSymmSCMatrix
HSOSSCF::alpha_density()
{
  RefSymmSCMatrix dens1(so_dimension(), basis_matrixkit());
  RefSymmSCMatrix dens2(so_dimension(), basis_matrixkit());

  so_density(dens1, 2.0);
  so_density(dens2, 1.0);
  dens1.accumulate(dens2);
  dens2=0;

  return dens1;
}

RefSymmSCMatrix
HSOSSCF::beta_density()
{
  RefSymmSCMatrix dens(so_dimension(), basis_matrixkit());
  so_density(dens, 2.0);
  return dens;
}

void
HSOSSCF::reset_density()
{
  cl_gmat_.assign(0.0);
  cl_dens_diff_.assign(cl_dens_);

  op_gmat_.assign(0.0);
  op_dens_diff_.assign(op_dens_);
}

double
HSOSSCF::new_density()
{
  // copy current density into density diff and scale by -1.  later we'll
  // add the new density to this to get the density difference.
  cl_dens_diff_.assign(cl_dens_);
  cl_dens_diff_.scale(-1.0);

  op_dens_diff_.assign(op_dens_);
  op_dens_diff_.scale(-1.0);

  so_density(cl_dens_, 2.0);
  cl_dens_.scale(2.0);

  so_density(op_dens_, 1.0);

  cl_dens_.accumulate(op_dens_);

  cl_dens_diff_.accumulate(cl_dens_);
  op_dens_diff_.accumulate(op_dens_);

  Ref<SCElementScalarProduct> sp(new SCElementScalarProduct);
  cl_dens_diff_.element_op(sp.pointer(), cl_dens_diff_);

  double delta = sp->result();
  delta = sqrt(delta/i_offset(cl_dens_diff_.n()));

  return delta;
}

RefSymmSCMatrix
HSOSSCF::density()
{
  if (!density_.computed()) {
    RefSymmSCMatrix dens(so_dimension(), basis_matrixkit());
    RefSymmSCMatrix dens1(so_dimension(), basis_matrixkit());
    so_density(dens, 2.0);
    dens.scale(2.0);

    so_density(dens1, 1.0);
    dens.accumulate(dens1);
    dens1=0;

    density_ = dens;
    // only flag the density as computed if the calc is converged
    if (!value_needed()) density_.computed() = 1;
  }

  return density_.result_noupdate();
}

double
HSOSSCF::scf_energy()
{
  RefSymmSCMatrix t = cl_fock_.result_noupdate().copy();
  t.accumulate(hcore_);

  RefSymmSCMatrix go = op_fock_.result_noupdate().copy();
  go.scale(-1.0);
  go.accumulate(cl_fock_.result_noupdate());

  SCFEnergy *eop = new SCFEnergy;
  eop->reference();
  Ref<SCElementOp2> op = eop;
  t.element_op(op, cl_dens_);

  double cl_e = eop->result();

  eop->reset();
  go.element_op(op, op_dens_);
  double op_e = eop->result();

  op=0;
  eop->dereference();
  delete eop;

  return cl_e-op_e;
}

Ref<SCExtrapData>
HSOSSCF::extrap_data()
{
  Ref<SCExtrapData> data =
    new SymmSCMatrix2SCExtrapData(cl_fock_.result_noupdate(),
                                  op_fock_.result_noupdate());
  return data;
}

RefSymmSCMatrix
HSOSSCF::effective_fock()
{
  // use fock() instead of cl_fock_ just in case this is called from
  // someplace outside SCF::compute_vector()
  RefSymmSCMatrix mofock(oso_dimension(), basis_matrixkit());
  mofock.assign(0.0);

  RefSymmSCMatrix mofocko(oso_dimension(), basis_matrixkit());
  mofocko.assign(0.0);

  // use eigenvectors if oso_scf_vector_ is null
  if (oso_scf_vector_.null()) {
    mofock.accumulate_transform(eigenvectors(), fock(0),
                                SCMatrix::TransposeTransform);
    mofocko.accumulate_transform(eigenvectors(), fock(1),
                                 SCMatrix::TransposeTransform);
  } else {
    RefSCMatrix so_to_oso_tr = so_to_orthog_so().t();
    mofock.accumulate_transform(so_to_oso_tr * oso_scf_vector_, fock(0),
                                SCMatrix::TransposeTransform);
    mofocko.accumulate_transform(so_to_oso_tr * oso_scf_vector_, fock(1),
                                 SCMatrix::TransposeTransform);
  }

  Ref<SCElementOp2> op = new GSGeneralEffH(this);
  mofock.element_op(op, mofocko);

  return mofock;
}

/////////////////////////////////////////////////////////////////////////////

void
HSOSSCF::init_gradient()
{
  // presumably the eigenvectors have already been computed by the time
  // we get here
  oso_scf_vector_ = oso_eigenvectors_.result_noupdate();
}

void
HSOSSCF::done_gradient()
{
  cl_dens_=0;
  op_dens_=0;
  oso_scf_vector_ = 0;
}

/////////////////////////////////////////////////////////////////////////////

// MO lagrangian
//       c    o   v
//  c  |2*FC|2*FC|0|
//     -------------
//  o  |2*FC| FO |0|
//     -------------
//  v  | 0  |  0 |0|
//
RefSymmSCMatrix
HSOSSCF::lagrangian()
{
  RefSCMatrix so_to_oso_tr = so_to_orthog_so().t();

  RefSymmSCMatrix mofock(oso_dimension(), basis_matrixkit());
  mofock.assign(0.0);
  mofock.accumulate_transform(so_to_oso_tr * oso_scf_vector_,
                              cl_fock_.result_noupdate(),
                              SCMatrix::TransposeTransform);

  RefSymmSCMatrix mofocko(oso_dimension(), basis_matrixkit());
  mofocko.assign(0.0);
  mofocko.accumulate_transform(so_to_oso_tr * oso_scf_vector_,
                               op_fock_.result_noupdate(),
                               SCMatrix::TransposeTransform);

  mofock.scale(2.0);

  Ref<SCElementOp2> op = new MOLagrangian(this);
  mofock.element_op(op, mofocko);
  mofocko=0;

  // transform MO lagrangian to SO basis
  RefSymmSCMatrix so_lag(so_dimension(), basis_matrixkit());
  so_lag.assign(0.0);
  so_lag.accumulate_transform(so_to_oso_tr * oso_scf_vector_, mofock);

  // and then from SO to AO
  Ref<PetiteList> pl = integral()->petite_list();
  RefSymmSCMatrix ao_lag = pl->to_AO_basis(so_lag);

  ao_lag.scale(-1.0);

  return ao_lag;
}

RefSymmSCMatrix
HSOSSCF::gradient_density()
{
  cl_dens_ = basis_matrixkit()->symmmatrix(so_dimension());
  op_dens_ = cl_dens_.clone();

  so_density(cl_dens_, 2.0);
  cl_dens_.scale(2.0);

  so_density(op_dens_, 1.0);

  Ref<PetiteList> pl = integral()->petite_list(basis());

  cl_dens_ = pl->to_AO_basis(cl_dens_);
  op_dens_ = pl->to_AO_basis(op_dens_);

  RefSymmSCMatrix tdens = cl_dens_.copy();
  tdens.accumulate(op_dens_);

  op_dens_.scale(2.0);

  return tdens;
}

/////////////////////////////////////////////////////////////////////////////

void
HSOSSCF::init_hessian()
{
}

void
HSOSSCF::done_hessian()
{
}

/////////////////////////////////////////////////////////////////////////////

void
HSOSSCF::two_body_deriv_hf(double * tbgrad, double exchange_fraction)
{
  Ref<SCElementMaxAbs> m = new SCElementMaxAbs;
  cl_dens_.element_op(m.pointer());
  op_dens_.element_op(m.pointer());
  double pmax = m->result();
  m=0;

  // now try to figure out the matrix specialization we're dealing with.
  // if we're using Local matrices, then there's just one subblock, or
  // see if we can convert P to local matrices

  if (local_ || local_dens_) {
    // grab the data pointers from the P matrices
    double *pmat, *pmato;
    RefSymmSCMatrix ptmp = get_local_data(cl_dens_, pmat, SCF::Read);
    RefSymmSCMatrix potmp = get_local_data(op_dens_, pmato, SCF::Read);

    Ref<PetiteList> pl = integral()->petite_list();
    LocalHSOSGradContribution l(pmat,pmato);

    int i;
    int na3 = molecule()->natom()*3;
    int nthread = threadgrp_->nthread();
    double **grads = new double*[nthread];
    Ref<TwoBodyDerivInt> *tbis = new Ref<TwoBodyDerivInt>[nthread];
    for (i=0; i < nthread; i++) {
      tbis[i] = integral()->electron_repulsion_deriv();
      grads[i] = new double[na3];
      memset(grads[i], 0, sizeof(double)*na3);
    }

    LocalTBGrad<LocalHSOSGradContribution> **tblds =
      new LocalTBGrad<LocalHSOSGradContribution>*[nthread];

    for (i=0; i < nthread; i++) {
      tblds[i] = new LocalTBGrad<LocalHSOSGradContribution>(
        l, tbis[i], pl, basis(), scf_grp_, grads[i], pmax,
        desired_gradient_accuracy(), nthread, i, exchange_fraction);
      threadgrp_->add_thread(i, tblds[i]);
    }

    if (threadgrp_->start_threads() < 0
        ||threadgrp_->wait_threads() < 0) {
      ExEnv::err0() << indent
           << "HSOSSCF: error running threads" << endl;
      abort();
    }

    for (i=0; i < nthread; i++) {
      for (int j=0; j < na3; j++)
        tbgrad[j] += grads[i][j];

      delete[] grads[i];
      delete tblds[i];
      tbis[i] = 0;
    }
    delete[] tbis;

    scf_grp_->sum(tbgrad, na3);
  }

  // for now quit
  else {
    ExEnv::err0() << indent
         << "HSOSSCF::two_body_deriv: can't do gradient yet\n";
    abort();
  }
}

/////////////////////////////////////////////////////////////////////////////

RefSCMatrix
HSOSSCF::alpha_semicanonical_eigenvectors()
{
  semicanonical();
  return alpha_semican_evecs_.result_noupdate();
}

RefSCMatrix
HSOSSCF::beta_semicanonical_eigenvectors()
{
  semicanonical();
  return beta_semican_evecs_.result_noupdate();
}

RefDiagSCMatrix
HSOSSCF::alpha_semicanonical_eigenvalues()
{
  semicanonical();
  return alpha_semican_evals_.result_noupdate();
}

RefDiagSCMatrix
HSOSSCF::beta_semicanonical_eigenvalues()
{
  semicanonical();
  return beta_semican_evals_.result_noupdate();
}

#define DEBUG_SEMICANONICAL 0
void
HSOSSCF::semicanonical()
{
  // May need to update this
  if (value_needed())
    compute();
  if (semicanonical_evaluated_)
    return;

  const RefSCDimension modim = eigenvalues().dim();
  const RefSCDimension sodim = so_dimension();
  const unsigned int nirrep = molecule()->point_group()->char_table().nirrep();

  RefSCMatrix so_eigenvector = eigenvectors();
#if DEBUG_SEMICANONICAL
  eigenvalues().print("HSOSSCF orbital energies");
  so_eigenvector.print("HSOSSCF orbitals");
#endif

  // initialize dimensions for occupied and virtual blocks
  int* aoccpi = new int[nirrep];  for(unsigned int h=0; h<nirrep; ++h) aoccpi[h] = 0;
  int* boccpi = new int[nirrep];  for(unsigned int h=0; h<nirrep; ++h) boccpi[h] = 0;
  int* avirpi = new int[nirrep];  for(unsigned int h=0; h<nirrep; ++h) avirpi[h] = 0;
  int* bvirpi = new int[nirrep];  for(unsigned int h=0; h<nirrep; ++h) bvirpi[h] = 0;
  for(unsigned int h=0; h<nirrep; ++h) {
    const unsigned int n = modim->blocks()->size(h);
    for(unsigned int i=0; i<n; ++i) {
      const double o = occupation(h,i);
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

  //
  // alpha case
  //
  // diagonalize occ-occ and virt-virt fock matrices
  // Fa = Fc - Go = Fo
  RefSymmSCMatrix afock_so = fock(0).clone(); afock_so.assign(fock(0)); afock_so.scale(0.0);
  afock_so.accumulate(fock(1));
#if DEBUG_SEMICANONICAL
  RefSymmSCMatrix Fa(oso_dimension(), basis_matrixkit());
  Fa.assign(0.0);
  Fa.accumulate_transform(so_eigenvector, afock_so,
                          SCMatrix::TransposeTransform);
  Fa.print("Alpha Fock matrix");
  effective_fock().print("Effective Fock matrix");
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
  if (alpha_semican_evals_.result_noupdate().null()) {
    alpha_semican_evals_.result_noupdate() = eigenvalues().clone();  alpha_semican_evals_.result_noupdate().assign(0.0);
    alpha_semican_evecs_.result_noupdate() = so_eigenvector.clone(); alpha_semican_evecs_.result_noupdate().assign(0.0);
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
        alpha_semican_evals_.result_noupdate().set_element(mo,e);
      }
      i = 0;
      for(unsigned int mo=mostart+aoccpi[h]; mo<moend; ++mo, ++i) {
        const double e = avirevals.get_element(aviroffset + i);
        alpha_semican_evals_.result_noupdate().set_element(mo,e);
      }

      aoccoffset += aoccpi[h];
      aviroffset += avirpi[h];
    }
  }
  alpha_semican_evecs_.result_noupdate().accumulate_product(so_eigenvector,aevecs);
#if DEBUG_SEMICANONICAL
  so_eigenvector.print("Original eigenvector");
  aevecs.print("Alpha transform matrix");
  {
    RefSymmSCMatrix afock_mo(modim, basis_matrixkit()); afock_mo.assign(0.0);
    afock_mo.accumulate_transform(alpha_semican_evecs_.result_noupdate(), afock_so,
                                  SCMatrix::TransposeTransform);
    afock_mo.print("Alpha Fock matrix in the semicanonical orbitals");
  }
  alpha_semican_evals_.result_noupdate().print("Alpha semicanonical orbital energies");
  alpha_semican_evecs_.result_noupdate().print("Alpha semicanonical orbitals");
#endif
  aevecs = 0;

  //
  // beta case
  //
  // diagonalize occ-occ and virt-virt fock matrices
  // Fb = Fc + Go = 2*Fc - Fo
  RefSymmSCMatrix bfock_so = fock(0).clone(); bfock_so.assign(fock(0)); bfock_so.scale(-2.0);
  bfock_so.accumulate(fock(1)); bfock_so.scale(-1.0);
#if DEBUG_SEMICANONICAL
  RefSymmSCMatrix Fb(oso_dimension(), basis_matrixkit());
  Fb.assign(0.0);
  Fb.accumulate_transform(so_eigenvector, bfock_so,
                          SCMatrix::TransposeTransform);
  Fb.print("Beta Fock matrix");
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
  if (beta_semican_evals_.result_noupdate().null()) {
    beta_semican_evals_.result_noupdate() = eigenvalues().clone();  beta_semican_evals_.result_noupdate().assign(0.0);
    beta_semican_evecs_.result_noupdate() = so_eigenvector.clone(); beta_semican_evecs_.result_noupdate().assign(0.0);
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
        beta_semican_evals_.result_noupdate().set_element(mo,e);
      }
      i = 0;
      for(unsigned int mo=mostart+boccpi[h]; mo<moend; ++mo, ++i) {
        const double e = bvirevals.get_element(bviroffset + i);
        beta_semican_evals_.result_noupdate().set_element(mo,e);
      }

      boccoffset += boccpi[h];
      bviroffset += bvirpi[h];
    }
  }
  beta_semican_evecs_.result_noupdate().accumulate_product(so_eigenvector,bevecs);
#if DEBUG_SEMICANONICAL
  bevecs.print("Beta transform matrix");
  {
    RefSymmSCMatrix bfock_mo(modim, basis_matrixkit()); bfock_mo.assign(0.0);
    bfock_mo.accumulate_transform(beta_semican_evecs_.result_noupdate(), bfock_so,
                                  SCMatrix::TransposeTransform);
    bfock_mo.print("Beta Fock matrix in the semicanonical orbitals");
  }
  beta_semican_evals_.result_noupdate().print("Beta semicanonical orbital energies");
  beta_semican_evecs_.result_noupdate().print("Beta semicanonical orbitals");
#endif
  bevecs = 0;

  semicanonical_evaluated_ = true;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:

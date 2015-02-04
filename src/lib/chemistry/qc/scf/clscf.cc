//
// clscf.cc --- implementation of the closed shell SCF class
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

#include <util/misc/regtime.h>
#include <util/misc/formio.h>
#include <util/state/stateio.h>

#include <math/scmat/block.h>
#include <math/scmat/blocked.h>
#include <math/scmat/blkiter.h>
#include <math/scmat/local.h>

#include <math/optimize/scextrapmat.h>

#include <chemistry/qc/basis/petite.h>

#include <chemistry/qc/scf/scflocal.h>
#include <chemistry/qc/scf/scfops.h>
#include <chemistry/qc/scf/clscf.h>
#include <chemistry/qc/scf/ltbgrad.h>
#include <chemistry/qc/scf/clhftmpl.h>
#include <chemistry/qc/scf/iter_logger.h>
#include <chemistry/qc/wfn/femo.h>

using namespace std;
using namespace sc;

namespace sc {

///////////////////////////////////////////////////////////////////////////
// CLSCF

static ClassDesc CLSCF_cd(
  typeid(CLSCF),"CLSCF",2,"public SCF",
  0, 0, 0);

CLSCF::CLSCF(StateIn& s) :
  SavableState(s),
  SCF(s),
  cl_fock_(this)
{
  cl_fock_.result_noupdate() =
    basis_matrixkit()->symmmatrix(so_dimension());
  cl_fock_.restore_state(s);
  cl_fock_.result_noupdate().restore(s);

  s.get(user_occupations_);
  s.get(tndocc_);
  s.get(nirrep_);
  s.get(ndocc_);
  if (s.version(::class_desc<CLSCF>()) >= 2) {
    s.get(initial_ndocc_);
    most_recent_pg_ << SavableState::restore_state(s);
  } else {
    initial_ndocc_ = new int[nirrep_];
    memcpy(initial_ndocc_, ndocc_, sizeof(int)*nirrep_);
  }

  // now take care of memory stuff
  init_mem(2);
}

CLSCF::CLSCF(const Ref<KeyVal>& keyval) :
  SCF(keyval),
  cl_fock_(this)
{
  int i;
  int me = scf_grp_->me();

  cl_fock_.compute()=0;
  cl_fock_.computed()=0;

  // calculate the total nuclear charge
  const int Znuc=molecule()->total_Z();

  // check to see if this is to be a charged molecule
  const int charge = keyval->intvalue("total_charge", KeyValValueint(0));

  int nelectron = Znuc-charge;

  // now see if ndocc was specified
  if (keyval->exists("ndocc")) {
    tndocc_ = keyval->intvalue("ndocc");
  } else {
    tndocc_ = nelectron/2;
    if (nelectron%2 && me==0) {
      ExEnv::err0() << endl
           << indent << "CLSCF::init: Warning, there's a leftover electron.\n"
           << incindent << indent << "total_charge = " << charge << endl
           << indent << "total nuclear charge = " << Znuc << endl
           << indent << "ndocc_ = " << tndocc_ << endl << decindent;
    }
  }

  ExEnv::out0() << endl << indent
       << "CLSCF::init: total charge = " << Znuc-2*tndocc_ << endl << endl;

  nirrep_ = molecule()->point_group()->char_table().ncomp();

  ndocc_ = read_occ(keyval, "docc", nirrep_);
  if (ndocc_) {
    user_occupations_=1;
    initial_ndocc_ = new int[nirrep_];
    memcpy(initial_ndocc_, ndocc_, sizeof(int)*nirrep_);
  } else {
    initial_ndocc_=0;
    ndocc_=0;
    user_occupations_=0;
    set_occupations(0);
  }

  ExEnv::out0() << indent << "docc = [";
  for (i=0; i < nirrep_; i++)
    ExEnv::out0() << " " << ndocc_[i];
  ExEnv::out0() << " ]\n";

  ExEnv::out0() << indent << "nbasis = " << basis()->nbasis() << endl;

  // check to see if this was done in SCF(keyval)
  if (!keyval->exists("maxiter"))
    maxiter_ = 100;

  // now take care of memory stuff
  init_mem(2);
}

CLSCF::~CLSCF()
{
  if (ndocc_) {
    delete[] ndocc_;
    ndocc_=0;
  }
  delete[] initial_ndocc_;
}

void
CLSCF::save_data_state(StateOut& s)
{
  SCF::save_data_state(s);

  cl_fock_.save_data_state(s);
  cl_fock_.result_noupdate().save(s);

  s.put(user_occupations_);
  s.put(tndocc_);
  s.put(nirrep_);
  s.put(ndocc_,nirrep_);
  s.put(initial_ndocc_,nirrep_);
  SavableState::save_state(most_recent_pg_.pointer(),s);
}

double
CLSCF::occupation(int ir, int i)
{
  if (i < ndocc_[ir]) return 2.0;
  return 0.0;
}

int
CLSCF::n_fock_matrices() const
{
  return 1;
}

RefSymmSCMatrix
CLSCF::fock(int n)
{
  if (n > 0) {
    ExEnv::err0() << indent
         << "CLSCF::fock: there is only one fock matrix, "
         << scprintf("but fock(%d) was requested\n",n);
    abort();
  }

  return cl_fock_.result();
}

double
CLSCF::magnetic_moment() const
{
  return 0.0;
}

void
CLSCF::print(ostream&o) const
{
  SCF::print(o);
  o << indent << "CLSCF Parameters:\n" << incindent
    << indent << "charge = " << molecule()->total_charge()-2*tndocc_ << endl
    << indent << "ndocc = " << tndocc_ << endl
    << indent << "docc = [";
  for (int i=0; i < nirrep_; i++)
    o << " " << ndocc_[i];
  o << " ]" << endl << decindent << endl;
}

//////////////////////////////////////////////////////////////////////////////

void
CLSCF::set_occupations(const RefDiagSCMatrix& ev)
{
  if (user_occupations_ || (initial_ndocc_ && ev.null())) {
    if (form_occupations(ndocc_, initial_ndocc_)) {
      most_recent_pg_ = new PointGroup(molecule()->point_group());
      return;
    }
    ExEnv::out0() << indent
         << "CLSCF: WARNING: reforming occupation vector from scratch" << endl;
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
  Ref<FEMO> femo = new FEMO(tndocc_, tndocc_,evals);
  // copy into local arrays
  int *newocc = new int[nirrep_];
  for(int g=0; g<nirrep_; ++g) {
    const int na = femo->nalpha(g);
    newocc[g] = na;
  }

  if (!ndocc_) {
    ndocc_=newocc;
  } else if (most_recent_pg_
             && most_recent_pg_->equiv(molecule()->point_group())) {
    // test to see if newocc is different from ndocc_
    for (i=0; i < nirrep_; i++) {
      if (ndocc_[i] != newocc[i]) {
        ExEnv::err0() << indent << "CLSCF::set_occupations:  WARNING!!!!\n"
             << incindent << indent
             << scprintf("occupations for irrep %d have changed\n",i+1)
             << indent
             << scprintf("ndocc was %d, changed to %d", ndocc_[i],newocc[i])
             << endl << decindent;
      }
    }

    memcpy(ndocc_,newocc,sizeof(int)*nirrep_);
    delete[] newocc;
  }

  if (!initial_ndocc_
      || initial_pg_->equiv(molecule()->point_group())) {
    delete[] initial_ndocc_;
    initial_ndocc_ = new int[nirrep_];
    memcpy(initial_ndocc_,ndocc_,sizeof(int)*nirrep_);
  }

  most_recent_pg_ = new PointGroup(molecule()->point_group());
}

void
CLSCF::symmetry_changed()
{
  SCF::symmetry_changed();
  cl_fock_.result_noupdate()=0;
  nirrep_ = molecule()->point_group()->char_table().ncomp();
  set_occupations(0);
}

//////////////////////////////////////////////////////////////////////////////
//
// scf things
//

void
CLSCF::init_vector()
{
  init_threads();

  // initialize the two electron integral classes
  ExEnv::out0() << indent
       << "integral intermediate storage = " << integral()->storage_used()
       << " bytes" << endl;
  ExEnv::out0() << indent
       << "integral cache = " << integral()->storage_unused()
       << " bytes" << endl;

  // allocate storage for other temp matrices
  cl_dens_ = hcore_.clone();
  cl_dens_.assign(0.0);

  cl_dens_diff_ = hcore_.clone();
  cl_dens_diff_.assign(0.0);

  // gmat is in AO basis
  cl_gmat_ = basis()->matrixkit()->symmmatrix(basis()->basisdim());
  cl_gmat_.assign(0.0);

  if (cl_fock_.result_noupdate().null()) {
    cl_fock_ = hcore_.clone();
    cl_fock_.result_noupdate().assign(0.0);
  }

  // make sure trial vector is set up
  initial_vector();

  oso_scf_vector_ = oso_eigenvectors_.result_noupdate();

  if (accumddh_) accumddh_->init(this);
}

void
CLSCF::done_vector()
{
  done_threads();

  if (accumddh_) {
      accumddh_->print_summary();
      accumddh_->done();
  }

  cl_gmat_ = 0;
  cl_dens_ = 0;
  cl_dens_diff_ = 0;

  oso_scf_vector_ = 0;
}

void
CLSCF::reset_density()
{
  cl_gmat_.assign(0.0);
  cl_dens_diff_.assign(cl_dens_);
}

RefSymmSCMatrix
CLSCF::density()
{
  if (!density_.computed()) {
    RefSymmSCMatrix dens(so_dimension(), basis_matrixkit());
    so_density(dens, 2.0);
    dens.scale(2.0);

    density_ = dens;
    // only flag the density as computed if the calc is converged
    if (!value_needed()) density_.computed() = 1;
  }

  return density_.result_noupdate();
}

double
CLSCF::new_density()
{
  // copy current density into density diff and scale by -1.  later we'll
  // add the new density to this to get the density difference.
  cl_dens_diff_.assign(cl_dens_);
  cl_dens_diff_.scale(-1.0);

  so_density(cl_dens_, 2.0);
  cl_dens_.scale(2.0);
  if(!iter_log_.null()){
    iter_log_->log_density(cl_dens_.copy());
  }

  cl_dens_diff_.accumulate(cl_dens_);

  Ref<SCElementScalarProduct> sp(new SCElementScalarProduct);
  cl_dens_diff_.element_op(sp.pointer(), cl_dens_diff_);

  double delta = sp->result();
  delta = sqrt(delta/i_offset(cl_dens_diff_.n()));

  return delta;
}

double
CLSCF::scf_energy()
{
  RefSymmSCMatrix t = cl_fock_.result_noupdate().copy();
  t.accumulate(hcore_);

  SCFEnergy *eop = new SCFEnergy;
  eop->reference();
  Ref<SCElementOp2> op = eop;
  t.element_op(op,cl_dens_);
  op=0;
  eop->dereference();

  double eelec = eop->result();

  delete eop;

  return eelec;
}

Ref<SCExtrapData>
CLSCF::initial_extrap_data()
{
  Ref<SCExtrapData> data;
  // If there is an old fock matrix around, use that.
  if (cl_fock_.result_noupdate()) {
    data = new SymmSCMatrixSCExtrapData(cl_fock_.result_noupdate());
  }
  return data;
}

Ref<SCExtrapData>
CLSCF::extrap_data()
{
  Ref<SCExtrapData> data =
    new SymmSCMatrixSCExtrapData(cl_fock_.result_noupdate());
  return data;
}

RefSymmSCMatrix
CLSCF::effective_fock()
{
  // just return MO fock matrix.  use fock() instead of cl_fock_ just in
  // case this is called from someplace outside SCF::compute_vector()
  RefSymmSCMatrix mofock(oso_dimension(), basis_matrixkit());
  mofock.assign(0.0);

  if (debug_ > 1) {
    fock(0).print("CL Fock matrix in SO basis");
  }

  // use eigenvectors if scf_vector_ is null
  if (oso_scf_vector_.null())
    mofock.accumulate_transform(eigenvectors(), fock(0),
                                SCMatrix::TransposeTransform);
  else
    mofock.accumulate_transform(so_to_orthog_so().t() * oso_scf_vector_,
                                fock(0),
                                SCMatrix::TransposeTransform);

  return mofock;
}

//////////////////////////////////////////////////////////////////////////////

void
CLSCF::init_gradient()
{
  // presumably the eigenvectors have already been computed by the time
  // we get here
  oso_scf_vector_ = oso_eigenvectors_.result_noupdate();
}

void
CLSCF::done_gradient()
{
  cl_dens_=0;
  oso_scf_vector_ = 0;
}

/////////////////////////////////////////////////////////////////////////////

class CLLag : public BlockedSCElementOp {
  private:
    CLSCF *scf_;

  public:
    CLLag(CLSCF* s) : scf_(s) {}
    ~CLLag() {}

    int has_side_effects() { return 1; }

    void process(SCMatrixBlockIter& bi) {
      int ir=current_block();

      for (bi.reset(); bi; bi++) {
        double occi = scf_->occupation(ir,bi.i());

        if (occi==0.0)
          bi.set(0.0);
      }
    }
};

RefSymmSCMatrix
CLSCF::lagrangian()
{
  // the MO lagrangian is just the eigenvalues of the occupied MO's
  RefSymmSCMatrix mofock = effective_fock();

  Ref<SCElementOp> op = new CLLag(this);
  mofock.element_op(op);

  // transform MO lagrangian to SO basis
  RefSymmSCMatrix so_lag(so_dimension(), basis_matrixkit());
  so_lag.assign(0.0);
  so_lag.accumulate_transform(so_to_orthog_so().t() * oso_scf_vector_, mofock);

  // and then from SO to AO
  Ref<PetiteList> pl = integral()->petite_list();
  RefSymmSCMatrix ao_lag = pl->to_AO_basis(so_lag);
  ao_lag->scale(-2.0);

  return ao_lag;
}

RefSymmSCMatrix
CLSCF::gradient_density()
{
  cl_dens_ = basis_matrixkit()->symmmatrix(so_dimension());
  cl_dens_.assign(0.0);

  so_density(cl_dens_, 2.0);
  cl_dens_.scale(2.0);

  Ref<PetiteList> pl = integral()->petite_list(basis());

  cl_dens_ = pl->to_AO_basis(cl_dens_);

  return cl_dens_;
}

/////////////////////////////////////////////////////////////////////////////

void
CLSCF::init_hessian()
{
}

void
CLSCF::done_hessian()
{
}

/////////////////////////////////////////////////////////////////////////////

void
CLSCF::two_body_deriv_hf(double * tbgrad, double exchange_fraction)
{
  int i;
  int na3 = molecule()->natom()*3;
  int nthread = threadgrp_->nthread();

  Timer tim("setup");
  Ref<SCElementMaxAbs> m = new SCElementMaxAbs;
  cl_dens_.element_op(m.pointer());
  double pmax = m->result();
  m=0;

  double **grads = new double*[nthread];
  Ref<TwoBodyDerivInt> *tbis = new Ref<TwoBodyDerivInt>[nthread];
  for (i=0; i < nthread; i++) {
    tbis[i] = integral()->electron_repulsion_deriv();
    grads[i] = new double[na3];
    memset(grads[i], 0, sizeof(double)*na3);
  }

  Ref<PetiteList> pl = integral()->petite_list();

  tim.change("contribution");

  // now try to figure out the matrix specialization we're dealing with.
  // if we're using Local matrices, then there's just one subblock, or
  // see if we can convert P to a local matrix
  if (local_ || local_dens_) {
    double *pmat;
    RefSymmSCMatrix ptmp = get_local_data(cl_dens_, pmat, SCF::Read);

    LocalCLHFGradContribution l(pmat);
    LocalTBGrad<LocalCLHFGradContribution> **tblds =
      new LocalTBGrad<LocalCLHFGradContribution>*[nthread];

    for (i=0; i < nthread; i++) {
      tblds[i] = new LocalTBGrad<LocalCLHFGradContribution>(
        l, tbis[i], pl, basis(), scf_grp_, grads[i], pmax,
        desired_gradient_accuracy(), nthread, i, exchange_fraction);
      threadgrp_->add_thread(i, tblds[i]);
    }

    tim.enter("start thread");
    if (threadgrp_->start_threads() < 0) {
      ExEnv::err0() << indent
           << "CLSCF: error starting threads" << endl;
      abort();
    }
    tim.exit("start thread");

    tim.enter("stop thread");
    if (threadgrp_->wait_threads() < 0) {
      ExEnv::err0() << indent
           << "CLSCF: error waiting for threads" << endl;
      abort();
    }
    tim.exit("stop thread");

    for (i=0; i < nthread; i++) {
      if (i) {
        for (int j=0; j < na3; j++)
          grads[0][j] += grads[i][j];

        delete[] grads[i];
      }

      delete tblds[i];
    }

    if (scf_grp_->n() > 1)
      scf_grp_->sum(grads[0], na3);

    for (i=0; i < na3; i++)
      tbgrad[i] += grads[0][i];

    delete[] grads[0];
    delete[] tblds;
    delete[] grads;
  }

  // for now quit
  else {
    ExEnv::err0() << indent
         << "CLHF::two_body_deriv: can't do gradient yet\n";
    abort();
  }

  for (i=0; i < nthread; i++)
    tbis[i] = 0;
  delete[] tbis;

  tim.exit("contribution");
}

/////////////////////////////////////////////////////////////////////////////

}

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:

//
// ossscf.cc --- implementation of the open shell singlet SCF class
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
#include <math/scmat/local.h>

#include <math/optimize/scextrapmat.h>

#include <chemistry/qc/basis/petite.h>

#include <chemistry/qc/scf/scflocal.h>
#include <chemistry/qc/scf/scfops.h>
#include <chemistry/qc/scf/effh.h>
#include <chemistry/qc/scf/ossscf.h>

using namespace std;
using namespace sc;

///////////////////////////////////////////////////////////////////////////
// OSSSCF

static ClassDesc OSSSCF_cd(
  typeid(OSSSCF),"OSSSCF",1,"public SCF",
  0, 0, 0);

OSSSCF::OSSSCF(StateIn& s) :
  SavableState(s),
  SCF(s),
  cl_fock_(this),
  op_focka_(this),
  op_fockb_(this)
{
  cl_fock_.result_noupdate() =
    basis_matrixkit()->symmmatrix(so_dimension());
  cl_fock_.restore_state(s);
  cl_fock_.result_noupdate().restore(s);

  op_focka_.result_noupdate() =
    basis_matrixkit()->symmmatrix(so_dimension());
  op_focka_.restore_state(s);
  op_focka_.result_noupdate().restore(s);

  op_fockb_.result_noupdate() =
    basis_matrixkit()->symmmatrix(so_dimension());
  op_fockb_.restore_state(s);
  op_fockb_.result_noupdate().restore(s);

  s.get(user_occupations_);
  s.get(tndocc_);
  s.get(nirrep_);
  s.get(ndocc_);
  s.get(osa_);
  s.get(osb_);

  // now take care of memory stuff
  init_mem(6);
}

OSSSCF::OSSSCF(const Ref<KeyVal>& keyval) :
  SCF(keyval),
  cl_fock_(this),
  op_focka_(this),
  op_fockb_(this)
{
  cl_fock_.compute()=0;
  cl_fock_.computed()=0;

  op_focka_.compute()=0;
  op_focka_.computed()=0;

  op_fockb_.compute()=0;
  op_fockb_.computed()=0;

  // calculate the total nuclear charge
  const int Znuc=molecule()->total_Z();

  // check to see if this is to be a charged molecule
  const int charge = keyval->intvalue("total_charge", KeyValValueint(0));
  const int nelectrons = Znuc-charge;

  // figure out how many doubly occupied shells there are
  if (keyval->exists("ndocc")) {
    tndocc_ = keyval->intvalue("ndocc");
  } else {
    tndocc_ = (nelectrons-2)/2;
    if ((nelectrons-2)%2) {
      ExEnv::err0() << endl << indent
           << "OSSSCF::init: Warning, there's a leftover electron.\n"
           << incindent
           << indent << "total_charge = " << charge << endl
           << indent << "total nuclear charge = " << Znuc << endl
           << indent << "ndocc_ = " << tndocc_ << endl << decindent;
    }
  }

  ExEnv::out0() << endl << indent << "OSSSCF::init: total charge = "
       << Znuc-2*tndocc_-2 << endl << endl;

  nirrep_ = molecule()->point_group()->order();

  if (nirrep_==1) {
    ExEnv::err0() << indent << "OSSSCF::init: cannot do C1 symmetry\n";
    abort();
  }

  osa_=-1;
  osb_=-1;

  ndocc_ = read_occ(keyval, "docc", nirrep_);
  int *nsocc = read_occ(keyval, "socc", nirrep_);
  if (ndocc_ && nsocc) {
    user_occupations_=1;
    for (int i=0; i < nirrep_; i++) {
      int nsi = nsocc[i];
      if (nsi && osa_<0)
        osa_=i;
      else if (nsi && osb_<0)
        osb_=i;
      else if (nsi) {
        ExEnv::err0() << indent << "OSSSCF::init: too many open shells\n";
        abort();
      }
    }
    delete[] nsocc;
  }
  else if ((ndocc_ && (!nsocc)) || ((!ndocc_) && nsocc)) {
    ExEnv::outn() << "ERROR: OSSSCF: only one of docc and socc specified: "
                 << "give both or none" << endl;
    abort();
  }
  else {
    ndocc_=0;
    user_occupations_=0;
    set_occupations(0);
  }

  int i;
  ExEnv::out0() << indent << "docc = [";
  for (i=0; i < nirrep_; i++)
    ExEnv::out0() << " " << ndocc_[i];
  ExEnv::out0() << " ]\n";

  ExEnv::out0() << indent << "socc = [";
  for (i=0; i < nirrep_; i++)
    ExEnv::out0() << " " << ((i==osa_ || i==osb_) ? 1 : 0);
  ExEnv::out0() << " ]\n";

  // check to see if this was done in SCF(keyval)
  if (!keyval->exists("maxiter"))
    maxiter_ = 200;

  if (!keyval->exists("level_shift"))
    level_shift_ = 0.25;

  // now take care of memory stuff
  init_mem(6);
}

OSSSCF::~OSSSCF()
{
  if (ndocc_) {
    delete[] ndocc_;
    ndocc_=0;
  }
}

void
OSSSCF::save_data_state(StateOut& s)
{
  SCF::save_data_state(s);

  cl_fock_.save_data_state(s);
  cl_fock_.result_noupdate().save(s);

  op_focka_.save_data_state(s);
  op_focka_.result_noupdate().save(s);

  op_fockb_.save_data_state(s);
  op_fockb_.result_noupdate().save(s);

  s.put(user_occupations_);
  s.put(tndocc_);
  s.put(nirrep_);
  s.put(ndocc_,nirrep_);
  s.put(osa_);
  s.put(osb_);
}

double
OSSSCF::occupation(int ir, int i)
{
  if (i < ndocc_[ir])
    return 2.0;
  else if ((ir==osa_ || ir==osb_) && (i == ndocc_[ir]))
    return 1.0;
  return 0.0;
}

double
OSSSCF::alpha_occupation(int ir, int i)
{
  if (i < ndocc_[ir] || (ir==osa_ && i==ndocc_[ir]))
    return 1.0;
  return 0.0;
}

double
OSSSCF::beta_occupation(int ir, int i)
{
  if (i < ndocc_[ir] || (ir==osb_ && i==ndocc_[ir]))
    return 1.0;
  return 0.0;
}

int
OSSSCF::n_fock_matrices() const
{
  return 3;
}

RefSymmSCMatrix
OSSSCF::fock(int n)
{
  if (n > 2) {
    ExEnv::err0() << indent
         << "OSSSCF::fock: there are only three fock matrices, "
         << scprintf("but fock(%d) was requested\n", n);
    abort();
  }

  if (n==0)
    return cl_fock_.result();
  else if (n==1)
    return op_focka_.result();
  else
    return op_fockb_.result();
}

double
OSSSCF::magnetic_moment() const
{
  return 0;
}

void
OSSSCF::print(ostream&o) const
{
  int i;

  SCF::print(o);
  o << indent << "OSSSCF Parameters:\n" << incindent
    << indent << "ndocc = " << tndocc_ << endl
    << indent << "docc = [";
  for (i=0; i < nirrep_; i++)
    o << " " << ndocc_[i];
  o << " ]" << endl;

  o << indent << "socc = [";
  for (i=0; i < nirrep_; i++)
    o << " " << ((i==osa_ || i==osb_) ? 1 : 0);
  o << " ]" << endl << decindent << endl;
}

//////////////////////////////////////////////////////////////////////////////

void
OSSSCF::set_occupations(const RefDiagSCMatrix& ev)
{
  if (user_occupations_)
    return;

  int i,j;

  RefDiagSCMatrix evals;

  if (ev == 0) {
    initial_vector();
    evals = eigenvalues_.result_noupdate();
  }
  else
    evals = ev;

  // first convert evals to something we can deal with easily
  BlockedDiagSCMatrix *evalsb = require_dynamic_cast<BlockedDiagSCMatrix*>(evals,
                                                 "OSSSCF::set_occupations");

  double **vals = new double*[nirrep_];
  for (i=0; i < nirrep_; i++) {
    int nf=oso_dimension()->blocks()->size(i);
    if (nf) {
      vals[i] = new double[nf];
      evalsb->block(i)->convert(vals[i]);
    } else {
      vals[i] = 0;
    }
  }

  // now loop to find the tndocc_ lowest eigenvalues and populate those
  // MO's
  int *newdocc = new int[nirrep_];
  memset(newdocc,0,sizeof(int)*nirrep_);

  for (i=0; i < tndocc_; i++) {
    // find lowest eigenvalue
    int lir=0,ln=0;
    double lowest=999999999;

    for (int ir=0; ir < nirrep_; ir++) {
      int nf=oso_dimension()->blocks()->size(ir);
      if (!nf)
        continue;
      for (j=0; j < nf; j++) {
        if (vals[ir][j] < lowest) {
          lowest=vals[ir][j];
          lir=ir;
          ln=j;
        }
      }
    }
    vals[lir][ln]=999999999;
    newdocc[lir]++;
  }

  int osa=-1, osb=-1;

  for (i=0; i < 2; i++) {
    // find lowest eigenvalue
    int lir=0,ln=0;
    double lowest=999999999;

    for (int ir=0; ir < nirrep_; ir++) {
      int nf=oso_dimension()->blocks()->size(ir);
      if (!nf)
        continue;
      for (j=0; j < nf; j++) {
        if (vals[ir][j] < lowest) {
          lowest=vals[ir][j];
          lir=ir;
          ln=j;
        }
      }
    }
    vals[lir][ln]=999999999;

    if (!i) {
      osa=lir;
    } else {
      if (lir==osa) {
        i--;
        continue;
      }
      osb=lir;
    }
  }

  // get rid of vals
  for (i=0; i < nirrep_; i++)
    if (vals[i])
      delete[] vals[i];
  delete[] vals;

  if (!ndocc_) {
    ndocc_=newdocc;
    osa_=osa;
    osb_=osb;
  } else {
    // test to see if newocc is different from ndocc_
    for (i=0; i < nirrep_; i++) {
      if (ndocc_[i] != newdocc[i]) {
        ExEnv::err0() << indent << "OSSSCF::set_occupations:  WARNING!!!!\n"
             << incindent << indent
             << scprintf("occupations for irrep %d have changed\n", i+1)
             << indent
             << scprintf("ndocc was %d, changed to %d", ndocc_[i], newdocc[i])
             << endl << decindent;
      }
      if ((osa != osa_ && osa != osb_) || (osb != osb_ && osb != osa_)) {
        ExEnv::err0() << indent << "OSSSCF::set_occupations:  WARNING!!!!\n"
             << incindent << indent
             << "open shell occupations have changed"
             << endl << decindent;
        osa_=osa;
        osb_=osb;
        reset_density();
      }
    }

    memcpy(ndocc_,newdocc,sizeof(int)*nirrep_);

    delete[] newdocc;
  }
}

void
OSSSCF::symmetry_changed()
{
  SCF::symmetry_changed();
  cl_fock_.result_noupdate()=0;
  op_focka_.result_noupdate()=0;
  op_fockb_.result_noupdate()=0;
  nirrep_ = molecule()->point_group()->char_table().ncomp();
  set_occupations(0);
}

//////////////////////////////////////////////////////////////////////////////
//
// scf things
//

void
OSSSCF::init_vector()
{
  init_threads();

  // allocate storage for other temp matrices
  cl_dens_ = hcore_.clone();
  cl_dens_.assign(0.0);

  cl_dens_diff_ = hcore_.clone();
  cl_dens_diff_.assign(0.0);

  op_densa_ = hcore_.clone();
  op_densa_.assign(0.0);

  op_densa_diff_ = hcore_.clone();
  op_densa_diff_.assign(0.0);

  op_densb_ = hcore_.clone();
  op_densb_.assign(0.0);

  op_densb_diff_ = hcore_.clone();
  op_densb_diff_.assign(0.0);

  // gmat is in AO basis
  cl_gmat_ = basis()->matrixkit()->symmmatrix(basis()->basisdim());
  cl_gmat_.assign(0.0);

  op_gmata_ = cl_gmat_.clone();
  op_gmata_.assign(0.0);

  op_gmatb_ = cl_gmat_.clone();
  op_gmatb_.assign(0.0);

  // test to see if we need a guess vector.
  if (cl_fock_.result_noupdate() == 0) {
    cl_fock_ = hcore_.clone();
    cl_fock_.result_noupdate().assign(0.0);
    op_focka_ = hcore_.clone();
    op_focka_.result_noupdate().assign(0.0);
    op_fockb_ = hcore_.clone();
    op_fockb_.result_noupdate().assign(0.0);
  }

  // make sure trial vector is set up
  initial_vector();

  oso_scf_vector_ = oso_eigenvectors_.result_noupdate();
}

void
OSSSCF::done_vector()
{
  done_threads();

  cl_gmat_ = 0;
  cl_dens_ = 0;
  cl_dens_diff_ = 0;
  op_gmata_ = 0;
  op_densa_ = 0;
  op_densa_diff_ = 0;
  op_gmatb_ = 0;
  op_densb_ = 0;
  op_densb_diff_ = 0;

  oso_scf_vector_ = 0;
}

RefSymmSCMatrix
OSSSCF::density()
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

RefSymmSCMatrix
OSSSCF::alpha_density()
{
  RefSymmSCMatrix dens1(so_dimension(), basis_matrixkit());
  RefSymmSCMatrix dens2(so_dimension(), basis_matrixkit());

  so_density(dens1, 2.0);
  so_density(dens2, 1.0);
  dynamic_cast<BlockedSymmSCMatrix*>(dens2.pointer())->block(osb_)->assign(0.0);

  dens1.accumulate(dens2);
  dens2=0;

  return dens1;
}

RefSymmSCMatrix
OSSSCF::beta_density()
{
  RefSymmSCMatrix dens1(so_dimension(), basis_matrixkit());
  RefSymmSCMatrix dens2(so_dimension(), basis_matrixkit());

  so_density(dens1, 2.0);
  so_density(dens2, 1.0);
  dynamic_cast<BlockedSymmSCMatrix*>(dens2.pointer())->block(osa_)->assign(0.0);

  dens1.accumulate(dens2);
  dens2=0;

  return dens1;
}

void
OSSSCF::reset_density()
{
  cl_gmat_.assign(0.0);
  cl_dens_diff_.assign(cl_dens_);

  op_gmata_.assign(0.0);
  op_densa_diff_.assign(op_densa_);

  op_gmatb_.assign(0.0);
  op_densb_diff_.assign(op_densb_);
}

double
OSSSCF::new_density()
{
  // copy current density into density diff and scale by -1.  later we'll
  // add the new density to this to get the density difference.
  cl_dens_diff_.assign(cl_dens_);
  cl_dens_diff_.scale(-1.0);

  op_densa_diff_.assign(op_densa_);
  op_densa_diff_.scale(-1.0);

  op_densb_diff_.assign(op_densb_);
  op_densb_diff_.scale(-1.0);

  so_density(cl_dens_, 2.0);
  cl_dens_.scale(2.0);

  so_density(op_densa_, 1.0);

  cl_dens_.accumulate(op_densa_);

  op_densb_.assign(op_densa_);
  dynamic_cast<BlockedSymmSCMatrix*>(op_densa_.pointer())->block(osb_)->assign(0.0);
  dynamic_cast<BlockedSymmSCMatrix*>(op_densb_.pointer())->block(osa_)->assign(0.0);

  cl_dens_diff_.accumulate(cl_dens_);
  op_densa_diff_.accumulate(op_densa_);
  op_densb_diff_.accumulate(op_densb_);

  Ref<SCElementScalarProduct> sp(new SCElementScalarProduct);
  cl_dens_diff_.element_op(sp.pointer(), cl_dens_diff_);

  double delta = sp->result();
  delta = sqrt(delta/i_offset(cl_dens_diff_.n()));

  return delta;
}

double
OSSSCF::scf_energy()
{
  RefSymmSCMatrix t = cl_fock_.result_noupdate().copy();
  t.accumulate(hcore_);

  RefSymmSCMatrix ga = op_focka_.result_noupdate().copy();
  ga.scale(-1.0);
  ga.accumulate(cl_fock_.result_noupdate());

  RefSymmSCMatrix gb = op_fockb_.result_noupdate().copy();
  gb.scale(-1.0);
  gb.accumulate(cl_fock_.result_noupdate());

  SCFEnergy *eop = new SCFEnergy;
  eop->reference();
  Ref<SCElementOp2> op = eop;
  t.element_op(op, cl_dens_);

  double cl_e = eop->result();

  eop->reset();
  ga.element_op(op, op_densa_);
  double opa_e = eop->result();

  eop->reset();
  gb.element_op(op, op_densb_);
  double opb_e = eop->result();

  op=0;
  eop->dereference();
  delete eop;

  return cl_e-opa_e-opb_e;
}

////////////////////////////////////////////////////////////////////////////

Ref<SCExtrapData>
OSSSCF::extrap_data()
{
  RefSymmSCMatrix *m = new RefSymmSCMatrix[3];
  m[0] = cl_fock_.result_noupdate();
  m[1] = op_focka_.result_noupdate();
  m[2] = op_fockb_.result_noupdate();

  Ref<SCExtrapData> data = new SymmSCMatrixNSCExtrapData(3, m);
  delete[] m;

  return data;
}

RefSymmSCMatrix
OSSSCF::effective_fock()
{
  // use fock() instead of cl_fock_ just in case this is called from
  // someplace outside SCF::compute_vector()
  RefSymmSCMatrix mofock(oso_dimension(), basis_matrixkit());
  mofock.assign(0.0);

  RefSymmSCMatrix mofocka(oso_dimension(), basis_matrixkit());
  mofocka.assign(0.0);

  RefSymmSCMatrix mofockb(oso_dimension(), basis_matrixkit());
  mofockb.assign(0.0);

  // use eigenvectors if oso_scf_vector_ is null
  RefSCMatrix vec;
  if (oso_scf_vector_ == 0) {
    vec = eigenvectors();
  } else {
    vec = so_to_orthog_so().t() * oso_scf_vector_;
  }

  mofock.accumulate_transform(vec, fock(0),
                              SCMatrix::TransposeTransform);
  mofocka.accumulate_transform(vec, fock(1),
                               SCMatrix::TransposeTransform);
  mofockb.accumulate_transform(vec, fock(2),
                               SCMatrix::TransposeTransform);

  dynamic_cast<BlockedSymmSCMatrix*>(mofocka.pointer())->block(osb_)->assign(0.0);
  dynamic_cast<BlockedSymmSCMatrix*>(mofockb.pointer())->block(osa_)->assign(0.0);

  mofocka.accumulate(mofockb);
  mofockb=0;

  Ref<SCElementOp2> op = new GSGeneralEffH(this);
  mofock.element_op(op, mofocka);

  return mofock;
}

/////////////////////////////////////////////////////////////////////////////

void
OSSSCF::init_gradient()
{
  // presumably the eigenvectors have already been computed by the time
  // we get here
  oso_scf_vector_ = oso_eigenvectors_.result_noupdate();
}

void
OSSSCF::done_gradient()
{
  cl_dens_=0;
  op_densa_=0;
  op_densb_=0;
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
OSSSCF::lagrangian()
{
  RefSCMatrix vec = so_to_orthog_so().t() * oso_scf_vector_;

  RefSymmSCMatrix mofock(oso_dimension(), basis_matrixkit());
  mofock.assign(0.0);
  mofock.accumulate_transform(vec, cl_fock_.result_noupdate(),
                              SCMatrix::TransposeTransform);

  RefSymmSCMatrix mofocka(oso_dimension(), basis_matrixkit());
  mofocka.assign(0.0);
  mofocka.accumulate_transform(vec, op_focka_.result_noupdate(),
                               SCMatrix::TransposeTransform);

  RefSymmSCMatrix mofockb(oso_dimension(), basis_matrixkit());
  mofockb.assign(0.0);
  mofockb.accumulate_transform(vec, op_fockb_.result_noupdate(),
                               SCMatrix::TransposeTransform);

  dynamic_cast<BlockedSymmSCMatrix*>(mofocka.pointer())->block(osb_)->assign(0.0);
  dynamic_cast<BlockedSymmSCMatrix*>(mofockb.pointer())->block(osa_)->assign(0.0);

  mofocka.accumulate(mofockb);
  mofockb=0;

  mofock.scale(2.0);

  Ref<SCElementOp2> op = new MOLagrangian(this);
  mofock.element_op(op, mofocka);
  mofocka=0;

  // transform MO lagrangian to SO basis
  RefSymmSCMatrix so_lag(so_dimension(), basis_matrixkit());
  so_lag.assign(0.0);
  so_lag.accumulate_transform(vec, mofock);

  // and then from SO to AO
  Ref<PetiteList> pl = integral()->petite_list();
  RefSymmSCMatrix ao_lag = pl->to_AO_basis(so_lag);

  ao_lag.scale(-1.0);

  return ao_lag;
}

RefSymmSCMatrix
OSSSCF::gradient_density()
{
  cl_dens_ = basis_matrixkit()->symmmatrix(so_dimension());
  op_densa_ = cl_dens_.clone();
  op_densb_ = cl_dens_.clone();

  so_density(cl_dens_, 2.0);
  cl_dens_.scale(2.0);

  so_density(op_densa_, 1.0);
  op_densb_.assign(op_densa_);

  dynamic_cast<BlockedSymmSCMatrix*>(op_densa_.pointer())->block(osb_)->assign(0.0);
  dynamic_cast<BlockedSymmSCMatrix*>(op_densb_.pointer())->block(osa_)->assign(0.0);

  Ref<PetiteList> pl = integral()->petite_list(basis());

  cl_dens_ = pl->to_AO_basis(cl_dens_);
  op_densa_ = pl->to_AO_basis(op_densa_);
  op_densb_ = pl->to_AO_basis(op_densb_);

  RefSymmSCMatrix tdens = cl_dens_.copy();
  tdens.accumulate(op_densa_);
  tdens.accumulate(op_densb_);

  op_densa_.scale(2.0);
  op_densb_.scale(2.0);

  return tdens;
}

/////////////////////////////////////////////////////////////////////////////

void
OSSSCF::init_hessian()
{
}

void
OSSSCF::done_hessian()
{
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:

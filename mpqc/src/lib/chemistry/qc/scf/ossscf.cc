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

#ifdef __GNUC__
#pragma implementation
#pragma implementation "osscont.h"
#endif

#include <math.h>

#include <util/misc/timer.h>
#include <util/misc/formio.h>

#include <math/scmat/block.h>
#include <math/scmat/blocked.h>
#include <math/scmat/local.h>

#include <math/optimize/scextrapmat.h>

#include <chemistry/qc/basis/petite.h>

#include <chemistry/qc/scf/scflocal.h>
#include <chemistry/qc/scf/scfden.h>
#include <chemistry/qc/scf/effh.h>

#include <chemistry/qc/scf/ossscf.h>
#include <chemistry/qc/scf/osscont.h>
#include <chemistry/qc/scf/lgbuild.h>
#include <chemistry/qc/scf/ltbgrad.h>

///////////////////////////////////////////////////////////////////////////

#ifdef __GNUC__
template class GBuild<LocalOSSContribution>;
template class LocalGBuild<LocalOSSContribution>;

template class TBGrad<LocalOSSGradContribution>;
template class LocalTBGrad<LocalOSSGradContribution>;
#endif

///////////////////////////////////////////////////////////////////////////
// OSSSCF

#define CLASSNAME OSSSCF
#define HAVE_STATEIN_CTOR
#define HAVE_KEYVAL_CTOR
#define PARENTS public SCF
#include <util/class/classi.h>
void *
OSSSCF::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SCF::_castdown(cd);
  return do_castdowns(casts,cd);
}

OSSSCF::OSSSCF(StateIn& s) :
  SCF(s),
  cl_fock_(this),
  op_focka_(this),
  op_fockb_(this)
  maybe_SavableState(s)
{
  cl_fock_.result_noupdate() =
    basis_matrixkit()->symmmatrix(basis_dimension());
  cl_fock_.restore_state(s);
  cl_fock_.result_noupdate().restore(s);
  
  op_focka_.result_noupdate() =
    basis_matrixkit()->symmmatrix(basis_dimension());
  op_focka_.restore_state(s);
  op_focka_.result_noupdate().restore(s);
  
  op_fockb_.result_noupdate() =
    basis_matrixkit()->symmmatrix(basis_dimension());
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

OSSSCF::OSSSCF(const RefKeyVal& keyval) :
  SCF(keyval),
  cl_fock_(this),
  op_focka_(this),
  op_fockb_(this)
{
  int me = scf_grp_->me();
  
  cl_fock_.compute()=0;
  cl_fock_.computed()=0;
  
  op_focka_.compute()=0;
  op_focka_.computed()=0;
  
  op_fockb_.compute()=0;
  op_fockb_.computed()=0;
  
  // calculate the total nuclear charge
  int Znuc=0;
  PointBag_double *z = molecule()->charges();
  
  for (Pix p=z->first(); p; z->next(p)) Znuc += (int) z->get(p);

  // check to see if this is to be a charged molecule
  int charge = keyval->intvalue("total_charge");
  int nelectrons = Znuc-charge;

  // figure out how many doubly occupied shells there are
  if (keyval->exists("ndocc")) {
    tndocc_ = keyval->intvalue("ndocc");
  } else {
    tndocc_ = (nelectrons-2)/2;
    if ((nelectrons-2)%2) {
      cerr << node0 << endl << indent
           << "OSSSCF::init: Warning, there's a leftover electron.\n"
           << incindent
           << indent << "total_charge = " << charge << endl
           << indent << "total nuclear charge = " << Znuc << endl
           << indent << "ndocc_ = " << tndocc_ << endl << decindent;
    }
  }

  cout << node0 << endl << indent << "OSSSCF::init: total charge = "
       << Znuc-2*tndocc_-2 << endl << endl;

  nirrep_ = molecule()->point_group().char_table().ncomp();

  if (nirrep_==1) {
    cerr << node0 << indent << "OSSSCF::init: cannot do C1 symmetry\n";
    abort();
  }

  osa_=-1;
  osb_=-1;
  
  if (keyval->exists("docc") && keyval->exists("socc")) {
    ndocc_ = new int[nirrep_];
    user_occupations_=1;
    for (int i=0; i < nirrep_; i++) {
      ndocc_[i] = keyval->intvalue("docc",i);
      int nsi = keyval->intvalue("socc",i);
      if (nsi && osa_<0)
        osa_==i;
      else if (nsi && osb_<0)
        osb_==i;
      else if (nsi) {
        cerr << node0 << indent << "OSSSCF::init: too many open shells\n";
        abort();
      }
    }
  } else {
    ndocc_=0;
    user_occupations_=0;
    set_occupations(0);
  }

  int i;
  cout << node0 << indent << "docc = [";
  for (i=0; i < nirrep_; i++)
    cout << node0 << " " << ndocc_[i];
  cout << node0 << " ]\n";

  cout << node0 << indent << "socc = [";
  for (i=0; i < nirrep_; i++)
    cout << node0 << " " << (i==osa_ || i==osb_) ? 1 : 0;
  cout << node0 << " ]\n";

  // check to see if this was done in SCF(keyval)
  if (!keyval->exists("maxiter"))
    maxiter_ = 200;

  if (!keyval->exists("level_shift"))
    level_shift_ = 1.0;

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
  else
    return 0.0;
}

int
OSSSCF::n_fock_matrices() const
{
  return 2;
}

RefSymmSCMatrix
OSSSCF::fock(int n)
{
  if (n > 2) {
    cerr << node0 << indent
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

int
OSSSCF::value_implemented()
{
  return 1;
}

int
OSSSCF::gradient_implemented()
{
  return 1;
}

int
OSSSCF::hessian_implemented()
{
  return 0;
}

void
OSSSCF::print(ostream&o)
{
  int i;
  
  SCF::print(o);
  o << node0 << indent << "OSSSCF Parameters:\n" << incindent
    << indent << "ndocc = " << tndocc_ << endl
    << indent << "docc = [";
  for (i=0; i < nirrep_; i++)
    o << node0 << " " << ndocc_[i];
  o << node0 << " ]" << endl;

  o << node0 << indent << "socc = [";
  for (i=0; i < nirrep_; i++)
    o << node0 << " " << (i==osa_ || i==osb_) ? 1 : 0;
  o << node0 << " ]" << endl << decindent << endl;
}

//////////////////////////////////////////////////////////////////////////////

void
OSSSCF::set_occupations(const RefDiagSCMatrix& ev)
{
  if (user_occupations_)
    return;
  
  int i,j;
  
  RefDiagSCMatrix evals;
  
  if (ev.null()) {
    initial_vector(0);
    evals = eigenvalues_.result_noupdate();
  }
  else
    evals = ev;

  // first convert evals to something we can deal with easily
  BlockedDiagSCMatrix *evalsb = BlockedDiagSCMatrix::require_castdown(evals,
                                                 "OSSSCF::set_occupations");
  
  RefPetiteList pl = integral()->petite_list(basis());
  
  double **vals = new double*[nirrep_];
  for (i=0; i < nirrep_; i++) {
    int nf=pl->nfunction(i);
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
    int lir,ln;
    double lowest=999999999;

    for (int ir=0; ir < nirrep_; ir++) {
      int nf=pl->nfunction(ir);
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
    int lir,ln;
    double lowest=999999999;

    for (int ir=0; ir < nirrep_; ir++) {
      int nf=pl->nfunction(ir);
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
        cerr << node0 << indent << "OSSSCF::set_occupations:  WARNING!!!!\n"
             << incindent << indent
             << scprintf("occupations for irrep %d have changed\n", i+1)
             << indent
             << scprintf("ndocc was %d, changed to %d", ndocc_[i], newdocc[i])
             << endl << decindent;
      }
      if ((osa != osa_ && osa != osb_) || (osb != osb_ && osb != osa_)) {
        cerr << node0 << indent << "OSSSCF::set_occupations:  WARNING!!!!\n"
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

//////////////////////////////////////////////////////////////////////////////
//
// scf things
//

void
OSSSCF::init_vector()
{
  // initialize the two electron integral classes
  tbi_ = integral()->electron_repulsion();
  tbi_->set_integral_storage(integral()->storage_unused());

  // calculate the core Hamiltonian
  cl_hcore_ = core_hamiltonian();
  
  // allocate storage for other temp matrices
  cl_dens_ = cl_hcore_.clone();
  cl_dens_.assign(0.0);
  
  cl_dens_diff_ = cl_hcore_.clone();
  cl_dens_diff_.assign(0.0);

  op_densa_ = cl_hcore_.clone();
  op_densa_.assign(0.0);
  
  op_densa_diff_ = cl_hcore_.clone();
  op_densa_diff_.assign(0.0);

  op_densb_ = cl_hcore_.clone();
  op_densb_.assign(0.0);
  
  op_densb_diff_ = cl_hcore_.clone();
  op_densb_diff_.assign(0.0);

  // gmat is in AO basis
  cl_gmat_ = basis()->matrixkit()->symmmatrix(basis()->basisdim());
  cl_gmat_.assign(0.0);

  op_gmata_ = cl_gmat_.clone();
  op_gmata_.assign(0.0);

  op_gmatb_ = cl_gmat_.clone();
  op_gmatb_.assign(0.0);

  // test to see if we need a guess vector.
  if (cl_fock_.result_noupdate().null()) {
    cl_fock_ = cl_hcore_.clone();
    cl_fock_.result_noupdate().assign(0.0);
    op_focka_ = cl_hcore_.clone();
    op_focka_.result_noupdate().assign(0.0);
    op_fockb_ = cl_hcore_.clone();
    op_fockb_.result_noupdate().assign(0.0);
  }

  // set up trial vector
  initial_vector(1);

  scf_vector_ = eigenvectors_.result_noupdate();
}

void
OSSSCF::done_vector()
{
  tbi_=0;
  
  cl_hcore_ = 0;
  cl_gmat_ = 0;
  cl_dens_ = 0;
  cl_dens_diff_ = 0;
  op_gmata_ = 0;
  op_densa_ = 0;
  op_densa_diff_ = 0;
  op_gmatb_ = 0;
  op_densb_ = 0;
  op_densb_diff_ = 0;

  scf_vector_ = 0;
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

  cl_dens_.assign(0.0);
  RefSCElementOp op = new SCFDensity(this, scf_vector_, 2.0);
  cl_dens_.element_op(op);
  cl_dens_.scale(2.0);

  op_densa_.assign(0.0);
  op = new SCFDensity(this, scf_vector_, 1.0);
  op_densa_.element_op(op);

  cl_dens_.accumulate(op_densa_);

  op_densb_.assign(op_densa_);
  BlockedSymmSCMatrix::castdown(op_densa_.pointer())->block(osb_)->assign(0.0);
  BlockedSymmSCMatrix::castdown(op_densb_.pointer())->block(osa_)->assign(0.0);
  
  cl_dens_diff_.accumulate(cl_dens_);
  op_densa_diff_.accumulate(op_densa_);
  op_densb_diff_.accumulate(op_densb_);

  RefSCElementScalarProduct sp(new SCElementScalarProduct);
  cl_dens_diff_.element_op(sp, cl_dens_diff_);
  
  double delta = sp->result();
  delta = sqrt(delta/i_offset(cl_dens_diff_.n()));

  return delta;
}

double
OSSSCF::scf_energy()
{
  RefSymmSCMatrix t = cl_fock_.result_noupdate().copy();
  t.accumulate(cl_hcore_);

  RefSymmSCMatrix ga = op_focka_.result_noupdate().copy();
  ga.scale(-1.0);
  ga.accumulate(cl_fock_.result_noupdate());
  
  RefSymmSCMatrix gb = op_fockb_.result_noupdate().copy();
  gb.scale(-1.0);
  gb.accumulate(cl_fock_.result_noupdate());
  
  SCFEnergy *eop = new SCFEnergy;
  eop->reference();
  RefSCElementOp2 op = eop;
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
    
RefSCExtrapData
OSSSCF::extrap_data()
{
  RefSymmSCMatrix *m = new RefSymmSCMatrix[3];
  m[0] = cl_fock_.result_noupdate();
  m[1] = op_focka_.result_noupdate();
  m[2] = op_fockb_.result_noupdate();
  
  RefSCExtrapData data = new SymmSCMatrixNSCExtrapData(3, m);
  delete[] m;
  
  return data;
}

RefSymmSCMatrix
OSSSCF::effective_fock()
{
  // use fock() instead of cl_fock_ just in case this is called from
  // someplace outside SCF::compute_vector()
  RefSymmSCMatrix mofock = fock(0).clone();
  mofock.assign(0.0);

  RefSymmSCMatrix mofocka = fock(1).clone();
  mofocka.assign(0.0);
  
  RefSymmSCMatrix mofockb = fock(2).clone();
  mofockb.assign(0.0);

  // use eigenvectors if scf_vector_ is null
  if (scf_vector_.null()) {
    mofock.accumulate_transform(eigenvectors().t(), fock(0));
    mofocka.accumulate_transform(eigenvectors().t(), fock(1));
    mofockb.accumulate_transform(eigenvectors().t(), fock(2));
  } else {
    mofock.accumulate_transform(scf_vector_.t(), fock(0));
    mofocka.accumulate_transform(scf_vector_.t(), fock(1));
    mofockb.accumulate_transform(scf_vector_.t(), fock(2));
  }
  
  BlockedSymmSCMatrix::castdown(mofocka.pointer())->block(osb_)->assign(0.0);
  BlockedSymmSCMatrix::castdown(mofockb.pointer())->block(osa_)->assign(0.0);
  
  mofocka.accumulate(mofockb);
  mofockb=0;
  
  RefSCElementOp2 op = new GSGeneralEffH(this);
  mofock.element_op(op, mofocka);

  return mofock;
}

//////////////////////////////////////////////////////////////////////////////

void
OSSSCF::ao_fock()
{
  RefPetiteList pl = integral()->petite_list(basis());
  
  // calculate G.  First transform cl_dens_diff_ to the AO basis, then
  // scale the off-diagonal elements by 2.0
  RefSymmSCMatrix dd = cl_dens_diff_;
  cl_dens_diff_ = pl->to_AO_basis(dd);
  cl_dens_diff_->scale(2.0);
  cl_dens_diff_->scale_diagonal(0.5);

  RefSymmSCMatrix dda = op_densa_diff_;
  op_densa_diff_ = pl->to_AO_basis(dda);
  op_densa_diff_->scale(2.0);
  op_densa_diff_->scale_diagonal(0.5);
  
  RefSymmSCMatrix ddb = op_densb_diff_;
  op_densb_diff_ = pl->to_AO_basis(ddb);
  op_densb_diff_->scale(2.0);
  op_densb_diff_->scale_diagonal(0.5);
  
  // now try to figure out the matrix specialization we're dealing with
  // if we're using Local matrices, then there's just one subblock, or
  // see if we can convert G and P to local matrices
  if (local_ || local_dens_) {
    // grab the data pointers from the G and P matrices
    double *gmat, *gmata, *gmatb, *pmat, *pmata, *pmatb;
    RefSymmSCMatrix gtmp = get_local_data(cl_gmat_, gmat, SCF::Accum);
    RefSymmSCMatrix ptmp = get_local_data(cl_dens_diff_, pmat, SCF::Read);
    RefSymmSCMatrix gatmp = get_local_data(op_gmata_, gmata, SCF::Accum);
    RefSymmSCMatrix patmp = get_local_data(op_densa_diff_, pmata, SCF::Read);
    RefSymmSCMatrix gbtmp = get_local_data(op_gmatb_, gmatb, SCF::Accum);
    RefSymmSCMatrix pbtmp = get_local_data(op_densb_diff_, pmatb, SCF::Read);
    
    signed char * pmax = init_pmax(pmat);
  
    LocalOSSContribution lclc(gmat, pmat, gmata, pmata, gmatb, pmatb);
    LocalGBuild<LocalOSSContribution>
      gb(lclc, tbi_, integral(), basis(), scf_grp_, pmax);
    gb.build_gmat(desired_value_accuracy()/100.0);

    delete[] pmax;

    // if we're running on multiple processors, then sum the G matrices
    if (scf_grp_->n() > 1) {
      scf_grp_->sum(gmat, i_offset(basis()->nbasis()));
      scf_grp_->sum(gmata, i_offset(basis()->nbasis()));
      scf_grp_->sum(gmatb, i_offset(basis()->nbasis()));
    }
    
    // if we're running on multiple processors, or we don't have local
    // matrices, then accumulate gtmp back into G
    if (!local_ || scf_grp_->n() > 1) {
      cl_gmat_->convert_accumulate(gtmp);
      op_gmata_->convert_accumulate(gatmp);
      op_gmatb_->convert_accumulate(gbtmp);
    }
  }

  // for now quit
  else {
    cerr << node0 << indent << "Cannot yet use anything but Local matrices\n";
    abort();
  }
  
  // get rid of AO delta P
  cl_dens_diff_ = dd;
  dd = cl_dens_diff_.clone();

  op_densa_diff_ = dda;
  dda = op_densa_diff_.clone();

  op_densb_diff_ = ddb;
  ddb = op_densb_diff_.clone();

  // now symmetrize the skeleton G matrix, placing the result in dd
  RefSymmSCMatrix skel_gmat = cl_gmat_.copy();
  skel_gmat.scale(1.0/(double)pl->order());
  pl->symmetrize(skel_gmat,dd);

  skel_gmat = op_gmata_.copy();
  skel_gmat.scale(1.0/(double)pl->order());
  pl->symmetrize(skel_gmat,dda);
  
  skel_gmat = op_gmatb_.copy();
  skel_gmat.scale(1.0/(double)pl->order());
  pl->symmetrize(skel_gmat,ddb);
  
  // F = H+G
  cl_fock_.result_noupdate().assign(cl_hcore_);
  cl_fock_.result_noupdate().accumulate(dd);

  // Fa = H+G-Ga
  op_focka_.result_noupdate().assign(cl_fock_.result_noupdate());
  dda.scale(-1.0);
  op_focka_.result_noupdate().accumulate(dda);

  // Fb = H+G-Gb
  op_fockb_.result_noupdate().assign(cl_fock_.result_noupdate());
  ddb.scale(-1.0);
  op_fockb_.result_noupdate().accumulate(ddb);

  cl_fock_.computed()=1;
  op_focka_.computed()=1;
  op_fockb_.computed()=1;
}

/////////////////////////////////////////////////////////////////////////////

void
OSSSCF::init_gradient()
{
  // presumably the eigenvectors have already been computed by the time
  // we get here
  scf_vector_ = eigenvectors_.result_noupdate();
}

void
OSSSCF::done_gradient()
{
  cl_dens_=0;
  op_densa_=0;
  op_densb_=0;
  scf_vector_ = 0;
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
  RefSymmSCMatrix mofock = cl_fock_.result_noupdate().clone();
  mofock.assign(0.0);
  mofock.accumulate_transform(scf_vector_.t(), cl_fock_.result_noupdate());

  RefSymmSCMatrix mofocka = mofock.clone();
  mofocka.assign(0.0);
  mofocka.accumulate_transform(scf_vector_.t(), op_focka_.result_noupdate());
  
  RefSymmSCMatrix mofockb = mofock.clone();
  mofockb.assign(0.0);
  mofockb.accumulate_transform(scf_vector_.t(), op_fockb_.result_noupdate());
  
  BlockedSymmSCMatrix::castdown(mofocka.pointer())->block(osb_)->assign(0.0);
  BlockedSymmSCMatrix::castdown(mofockb.pointer())->block(osa_)->assign(0.0);
  
  mofocka.accumulate(mofockb);
  mofockb=0;
  
  mofock.scale(2.0);
  
  RefSCElementOp2 op = new MOLagrangian(this);
  mofock.element_op(op, mofocka);
  mofocka=0;

  // transform MO lagrangian to SO basis
  RefSymmSCMatrix so_lag(basis_dimension(), basis_matrixkit());
  so_lag.assign(0.0);
  so_lag.accumulate_transform(scf_vector_, mofock);
  
  // and then from SO to AO
  RefPetiteList pl = integral()->petite_list();
  RefSymmSCMatrix ao_lag = pl->to_AO_basis(so_lag);

  ao_lag.scale(-1.0);

  return ao_lag;
}

RefSymmSCMatrix
OSSSCF::gradient_density()
{
  cl_dens_ = basis_matrixkit()->symmmatrix(basis_dimension());
  op_densa_ = cl_dens_.clone();
  op_densb_ = cl_dens_.clone();
  
  cl_dens_.assign(0.0);
  op_densa_.assign(0.0);
  op_densb_.assign(0.0);
  
  RefSCElementOp op = new SCFDensity(this, scf_vector_, 2.0);
  cl_dens_.element_op(op);
  cl_dens_.scale(2.0);
  
  op = new SCFDensity(this, scf_vector_, 1.0);
  op_densa_.element_op(op);
  op_densb_.assign(op_densa_);

  BlockedSymmSCMatrix::castdown(op_densa_.pointer())->block(osb_)->assign(0.0);
  BlockedSymmSCMatrix::castdown(op_densb_.pointer())->block(osa_)->assign(0.0);
  
  RefPetiteList pl = integral()->petite_list(basis());
  
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
OSSSCF::two_body_deriv(double * tbgrad)
{
  RefSCElementMaxAbs m = new SCElementMaxAbs;
  cl_dens_.element_op(m);
  double pmax = m->result();
  m=0;

  // now try to figure out the matrix specialization we're dealing with.
  // if we're using Local matrices, then there's just one subblock, or
  // see if we can convert P to local matrices

  if (local_ || local_dens_) {
    // grab the data pointers from the P matrices
    double *pmat, *pmata, *pmatb;
    RefSymmSCMatrix ptmp = get_local_data(cl_dens_, pmat, SCF::Read);
    RefSymmSCMatrix patmp = get_local_data(op_densa_, pmata, SCF::Read);
    RefSymmSCMatrix pbtmp = get_local_data(op_densb_, pmatb, SCF::Read);
  
    LocalOSSGradContribution l(pmat,pmata,pmatb);
    LocalTBGrad<LocalOSSGradContribution> tb(l, integral(), basis(), scf_grp_);
    tb.build_tbgrad(tbgrad, pmax, desired_gradient_accuracy());
  }

  // for now quit
  else {
    cerr << node0 << indent
         << "OSSSCF::two_body_deriv: can't do gradient yet\n";
    abort();
  }
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

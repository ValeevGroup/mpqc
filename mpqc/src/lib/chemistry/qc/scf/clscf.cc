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

#ifdef __GNUC__
#pragma implementation
#pragma implementation "clcont.h"
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
#include <chemistry/qc/scf/clscf.h>
#include <chemistry/qc/scf/clcont.h>
#include <chemistry/qc/scf/lgbuild.h>
#include <chemistry/qc/scf/lbgbuild.h>
#include <chemistry/qc/scf/ltbgrad.h>

///////////////////////////////////////////////////////////////////////////

#ifdef __GNUC__
template class GBuild<LocalCLContribution>;
template class LocalGBuild<LocalCLContribution>;
template class LocalLBGBuild<LocalCLContribution>;

template class TBGrad<LocalCLGradContribution>;
template class LocalTBGrad<LocalCLGradContribution>;
#endif

///////////////////////////////////////////////////////////////////////////
// CLSCF

#define CLASSNAME CLSCF
#define HAVE_STATEIN_CTOR
#define HAVE_KEYVAL_CTOR
#define PARENTS public SCF
#include <util/class/classi.h>
void *
CLSCF::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SCF::_castdown(cd);
  return do_castdowns(casts,cd);
}

CLSCF::CLSCF(StateIn& s) :
  SCF(s),
  cl_fock_(this)
  maybe_SavableState(s)
{
  cl_fock_.result_noupdate() =
    basis_matrixkit()->symmmatrix(basis_dimension());
  cl_fock_.restore_state(s);
  cl_fock_.result_noupdate().restore(s);
  
  s.get(user_occupations_);
  s.get(tndocc_);
  s.get(nirrep_);
  s.get(ndocc_);

  // now take care of memory stuff
  init_mem(2);
}

CLSCF::CLSCF(const RefKeyVal& keyval) :
  SCF(keyval),
  cl_fock_(this)
{
  int i;
  int me = scf_grp_->me();
  
  cl_fock_.compute()=0;
  cl_fock_.computed()=0;
  
  // calculate the total nuclear charge
  int Znuc=0;
  PointBag_double *z = molecule()->charges();
  
  for (Pix p=z->first(); p; z->next(p)) Znuc += (int) z->get(p);

  // check to see if this is to be a charged molecule
  int charge = keyval->intvalue("total_charge");
  
  // now see if ndocc was specified
  if (keyval->exists("ndocc")) {
    tndocc_ = keyval->intvalue("ndocc");
  } else {
    tndocc_ = (Znuc-charge)/2;
    if ((Znuc-charge)%2 && me==0) {
      cerr << node0 << endl
           << indent << "CLSCF::init: Warning, there's a leftover electron.\n"
           << incindent << indent << "total_charge = " << charge << endl
           << indent << "total nuclear charge = " << Znuc << endl
           << indent << "ndocc_ = " << tndocc_ << endl << decindent;
    }
  }

  cout << node0 << endl << indent
       << "CLSCF::init: total charge = " << Znuc-2*tndocc_ << endl << endl;

  nirrep_ = molecule()->point_group().char_table().ncomp();

  if (keyval->exists("docc")) {
    ndocc_ = new int[nirrep_];
    user_occupations_=1;
    for (i=0; i < nirrep_; i++) {
      ndocc_[i]=0;

      if (keyval->exists("docc",i))
        ndocc_[i] = keyval->intvalue("docc",i);
    }
  } else {
    ndocc_=0;
    user_occupations_=0;
    set_occupations(0);
  }

  cout << node0 << indent << "docc = [";
  for (i=0; i < nirrep_; i++)
    cout << node0 << " " << ndocc_[i];
  cout << node0 << " ]\n";

  cout << node0 << indent << "nbasis = " << basis()->nbasis() << endl;

  // check to see if this was done in SCF(keyval)
  if (!keyval->exists("maxiter"))
    maxiter_ = 40;

  // now take care of memory stuff
  init_mem(2);
}

CLSCF::~CLSCF()
{
  if (ndocc_) {
    delete[] ndocc_;
    ndocc_=0;
  }
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
    cerr << node0 << indent
         << "CLSCF::fock: there is only one fock matrix, "
         << scprintf("but fock(%d) was requested\n",n);
    abort();
  }

  return cl_fock_.result();
}

int
CLSCF::value_implemented()
{
  return 1;
}

int
CLSCF::gradient_implemented()
{
  return 1;
}

int
CLSCF::hessian_implemented()
{
  return 0;
}

void
CLSCF::print(ostream&o)
{
  SCF::print(o);
  o << node0 << indent << "CLSCF Parameters:\n" << incindent
    << indent << "ndocc = " << tndocc_ << endl
    << indent << "docc = [";
  for (int i=0; i < nirrep_; i++)
    o << node0 << " " << ndocc_[i];
  o << node0 << " ]" << endl << decindent << endl;
}

//////////////////////////////////////////////////////////////////////////////

void
CLSCF::set_occupations(const RefDiagSCMatrix& ev)
{
  if (user_occupations_)
    return;
  
  if (nirrep_==1) {
    if (!ndocc_) {
      ndocc_=new int[1];
      ndocc_[0]=tndocc_;
    }
    return;
  }
  
  int i,j;
  
  RefDiagSCMatrix evals;
  
  if (ev.null())
    evals = core_hamiltonian().eigvals();
  else
    evals = ev;

  // first convert evals to something we can deal with easily
  BlockedDiagSCMatrix *evalsb = BlockedDiagSCMatrix::require_castdown(evals,
                                                 "CLSCF::set_occupations");
  
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
  int *newocc = new int[nirrep_];
  memset(newocc,0,sizeof(int)*nirrep_);

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
    newocc[lir]++;
  }

  // get rid of vals
  for (i=0; i < nirrep_; i++)
    if (vals[i])
      delete[] vals[i];
  delete[] vals;

  if (!ndocc_) {
    ndocc_=newocc;
  } else {
    // test to see if newocc is different from ndocc_
    for (i=0; i < nirrep_; i++) {
      if (ndocc_[i] != newocc[i]) {
        cerr << node0 << indent << "CLSCF::set_occupations:  WARNING!!!!\n"
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
}

//////////////////////////////////////////////////////////////////////////////
//
// scf things
//

void
CLSCF::init_vector()
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

  // gmat is in AO basis
  cl_gmat_ = basis()->matrixkit()->symmmatrix(basis()->basisdim());
  cl_gmat_.assign(0.0);

  // test to see if we need a guess vector.
  if (eigenvectors_.result_noupdate().null()) {
    initial_vector();
    
    cl_fock_ = cl_hcore_.clone();
    cl_fock_.result_noupdate().assign(0.0);
  }

  scf_vector_ = eigenvectors_.result_noupdate();
}

void
CLSCF::done_vector()
{
  tbi_=0;
  
  cl_hcore_ = 0;
  cl_gmat_ = 0;
  cl_dens_ = 0;
  cl_dens_diff_ = 0;

  scf_vector_ = 0;
}

void
CLSCF::reset_density()
{
  cl_gmat_.assign(0.0);
  cl_dens_diff_.assign(cl_dens_);
}

double
CLSCF::new_density()
{
  // copy current density into density diff and scale by -1.  later we'll
  // add the new density to this to get the density difference.
  cl_dens_diff_.assign(cl_dens_);
  cl_dens_diff_.scale(-1.0);
  
  cl_dens_.assign(0.0);
  RefSCElementOp op = new SCFDensity(this, scf_vector_, 2.0);
  cl_dens_.element_op(op);
  cl_dens_.scale(2.0);

  cl_dens_diff_.accumulate(cl_dens_);
  
  RefSCElementScalarProduct sp(new SCElementScalarProduct);
  cl_dens_diff_.element_op(sp, cl_dens_diff_);
  
  double delta = sp->result();
  delta = sqrt(delta/i_offset(cl_dens_diff_.n()));

  return delta;
}

double
CLSCF::scf_energy()
{
  RefSymmSCMatrix t = cl_fock_.result_noupdate().copy();
  t.accumulate(cl_hcore_);

  SCFEnergy *eop = new SCFEnergy;
  eop->reference();
  RefSCElementOp2 op = eop;
  t.element_op(op,cl_dens_);
  op=0;
  eop->dereference();

  double eelec = eop->result();

  delete eop;
  
  return eelec;
}

RefSCExtrapData
CLSCF::extrap_data()
{
  RefSCExtrapData data =
    new SymmSCMatrixSCExtrapData(cl_fock_.result_noupdate());
  return data;
}

RefSymmSCMatrix
CLSCF::effective_fock()
{
  // just return MO fock matrix.  use fock() instead of cl_fock_ just in
  // case this is called from someplace outside SCF::compute_vector()
  RefSymmSCMatrix mofock = fock(0).clone();
  mofock.assign(0.0);

  // use eigenvectors if scf_vector_ is null
  if (scf_vector_.null())
    mofock.accumulate_transform(eigenvectors().t(), fock(0));
  else
    mofock.accumulate_transform(scf_vector_.t(), fock(0));

  return mofock;
}

//////////////////////////////////////////////////////////////////////////////

void
CLSCF::ao_fock()
{
  tim_enter("petite");
  RefPetiteList pl = integral()->petite_list(basis());
  tim_exit("petite");
  
  // calculate G.  First transform cl_dens_diff_ to the AO basis, then
  // scale the off-diagonal elements by 2.0
  tim_enter("setup");
  RefSymmSCMatrix dd = cl_dens_diff_;
  cl_dens_diff_ = pl->to_AO_basis(dd);
  cl_dens_diff_->scale(2.0);
  cl_dens_diff_->scale_diagonal(0.5);
  tim_exit("setup");

  // now try to figure out the matrix specialization we're dealing with
  // if we're using Local matrices, then there's just one subblock, or
  // see if we can convert G and P to local matrices

  if (local_ || local_dens_) {
    // grab the data pointers from the G and P matrices
    double *gmat, *pmat;
    tim_enter("local data");
    RefSymmSCMatrix gtmp = get_local_data(cl_gmat_, gmat, SCF::Accum);
    RefSymmSCMatrix ptmp = get_local_data(cl_dens_diff_, pmat, SCF::Read);
    tim_exit("local data");

    tim_enter("init pmax");
    signed char * pmax = init_pmax(pmat);
    tim_exit("init pmax");
  
    LocalCLContribution lclc(gmat, pmat);
    LocalGBuild<LocalCLContribution>
      gb(lclc, tbi_, integral(), basis(), scf_grp_, pmax);
    gb.build_gmat(desired_value_accuracy()/100.0);

    delete[] pmax;

    // if we're running on multiple processors, then sum the G matrix
    tim_enter("sum");
    if (scf_grp_->n() > 1)
      scf_grp_->sum(gmat, i_offset(basis()->nbasis()));
    tim_exit("sum");

    // if we're running on multiple processors, or we don't have local
    // matrices, then accumulate gtmp back into G
    tim_enter("accum");
    if (!local_ || scf_grp_->n() > 1)
      cl_gmat_->convert_accumulate(gtmp);
    tim_exit("accum");
  }

  // for now quit
  else {
    cerr << node0 << indent << "Cannot yet use anything but Local matrices\n";
    abort();
  }
  
  tim_enter("symm");
  // get rid of AO delta P
  cl_dens_diff_ = dd;
  dd = cl_dens_diff_.clone();

  // now symmetrize the skeleton G matrix, placing the result in dd
  RefSymmSCMatrix skel_gmat = cl_gmat_.copy();
  skel_gmat.scale(1.0/(double)pl->order());
  pl->symmetrize(skel_gmat,dd);
  tim_exit("symm");
  
  // F = H+G
  cl_fock_.result_noupdate().assign(cl_hcore_);
  cl_fock_.result_noupdate().accumulate(dd);
  cl_fock_.computed()=1;
}

/////////////////////////////////////////////////////////////////////////////

void
CLSCF::init_gradient()
{
  // presumably the eigenvectors have already been computed by the time
  // we get here
  scf_vector_ = eigenvectors_.result_noupdate();
}

void
CLSCF::done_gradient()
{
  cl_dens_=0;
  scf_vector_ = 0;
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



  RefSCElementOp op = new CLLag(this);
  mofock.element_op(op);
  
  // transform MO lagrangian to SO basis
  RefSymmSCMatrix so_lag(basis_dimension(), basis_matrixkit());
  so_lag.assign(0.0);
  so_lag.accumulate_transform(scf_vector_, mofock);
  
  // and then from SO to AO
  RefPetiteList pl = integral()->petite_list();
  RefSymmSCMatrix ao_lag = pl->to_AO_basis(so_lag);
  ao_lag->scale(-2.0);

  return ao_lag;
}

RefSymmSCMatrix
CLSCF::gradient_density()
{
  cl_dens_ = basis_matrixkit()->symmmatrix(basis_dimension());
  cl_dens_.assign(0.0);
  
  RefSCElementOp op = new SCFDensity(this, scf_vector_, 2.0);
  cl_dens_.element_op(op);
  cl_dens_.scale(2.0);
  
  RefPetiteList pl = integral()->petite_list(basis());
  
  cl_dens_ = pl->to_AO_basis(cl_dens_);
  cl_dens_.scale(1.0);

  return cl_dens_;
}

/////////////////////////////////////////////////////////////////////////////

void
CLSCF::two_body_deriv(double * tbgrad)
{
  RefSCElementMaxAbs m = new SCElementMaxAbs;
  cl_dens_.element_op(m);
  double pmax = m->result();
  m=0;

  // now try to figure out the matrix specialization we're dealing with.
  // if we're using Local matrices, then there's just one subblock, or
  // see if we can convert P to a local matrix
  if (local_ || local_dens_) {
    double *pmat;
    RefSymmSCMatrix ptmp = get_local_data(cl_dens_, pmat, SCF::Read);

    LocalCLGradContribution l(pmat);
    LocalTBGrad<LocalCLGradContribution> tb(l, integral(), basis(), scf_grp_);
    tb.build_tbgrad(tbgrad, pmax, desired_gradient_accuracy());
  }

  // for now quit
  else {
    cerr << node0 << indent
         << "CLSCF::two_body_deriv: can't do gradient yet\n";
    abort();
  }
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

// Local Variables:
// mode: c++
// eval: (c-set-style "ETS")
// End:

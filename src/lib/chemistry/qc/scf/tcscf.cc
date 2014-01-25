//
// tcscf.cc --- implementation of the two-configuration SCF class
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
#include <chemistry/qc/scf/effh.h>

#include <chemistry/qc/scf/tcscf.h>

using namespace std;
using namespace sc;

///////////////////////////////////////////////////////////////////////////
// TCSCF

static ClassDesc TCSCF_cd(
  typeid(TCSCF),"TCSCF",1,"public SCF",
  0, 0, 0);

TCSCF::TCSCF(StateIn& s) :
  SavableState(s),
  SCF(s),
  focka_(this),
  fockb_(this),
  ka_(this),
  kb_(this)
{
  focka_.result_noupdate() = basis_matrixkit()->symmmatrix(so_dimension());
  focka_.restore_state(s);
  focka_.result_noupdate().restore(s);

  fockb_.result_noupdate() = basis_matrixkit()->symmmatrix(so_dimension());
  fockb_.restore_state(s);
  fockb_.result_noupdate().restore(s);

  ka_.result_noupdate() = basis_matrixkit()->symmmatrix(so_dimension());
  ka_.restore_state(s);
  ka_.result_noupdate().restore(s);

  kb_.result_noupdate() = basis_matrixkit()->symmmatrix(so_dimension());
  kb_.restore_state(s);
  kb_.result_noupdate().restore(s);

  s.get(user_occupations_);
  s.get(tndocc_);
  s.get(nirrep_);
  s.get(ndocc_);
  s.get(osa_);
  s.get(osb_);
  s.get(occa_);
  s.get(occb_);
  s.get(ci1_);
  s.get(ci2_);

  // now take care of memory stuff
  init_mem(8);
}

TCSCF::TCSCF(const Ref<KeyVal>& keyval) :
  SCF(keyval),
  focka_(this),
  fockb_(this),
  ka_(this),
  kb_(this)
{
  focka_.compute()=0;
  focka_.computed()=0;

  fockb_.compute()=0;
  fockb_.computed()=0;

  ka_.compute()=0;
  ka_.computed()=0;

  kb_.compute()=0;
  kb_.computed()=0;

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
           << "TCSCF::init: Warning, there's a leftover electron.\n"
           << incindent
           << indent << "total_charge = " << charge << endl
           << indent << "total nuclear charge = " << Znuc << endl
           << indent << "ndocc_ = " << tndocc_ << endl << decindent;
    }
  }

  ExEnv::out0() << endl << indent << "TCSCF::init: total charge = "
       << Znuc-2*tndocc_-2 << endl << endl;

  nirrep_ = molecule()->point_group()->char_table().ncomp();

  if (nirrep_==1) {
    ExEnv::err0() << indent << "TCSCF::init: cannot do C1 symmetry\n";
    abort();
  }

  occa_=occb_=1.0;
  ci1_=ci2_ = 0.5*sqrt(2.0);

  if (keyval->exists("ci1")) {
    ci1_ = keyval->doublevalue("ci1");
    ci2_ = sqrt(1.0 - ci1_*ci1_);
    occa_ = 2.0*ci1_*ci1_;
    occb_ = 2.0*ci2_*ci2_;
  }

  if (keyval->exists("occa")) {
    occa_ = keyval->doublevalue("occa");
    ci1_ = sqrt(occa_/2.0);
    ci2_ = sqrt(1.0 - ci1_*ci1_);
    occb_ = 2.0*ci2_*ci2_;
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
        ExEnv::err0() << indent << "TCSCF::init: too many open shells\n";
        abort();
      }
    }
    delete[] nsocc;
  }
  else if ((ndocc_ && (!nsocc)) || ((!ndocc_) && nsocc)) {
    ExEnv::outn() << "ERROR: TCSCF: only one of docc and socc specified: "
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
  init_mem(8);
}

TCSCF::~TCSCF()
{
  if (ndocc_) {
    delete[] ndocc_;
    ndocc_=0;
  }
}

void
TCSCF::save_data_state(StateOut& s)
{
  SCF::save_data_state(s);

  focka_.save_data_state(s);
  focka_.result_noupdate().save(s);

  fockb_.save_data_state(s);
  fockb_.result_noupdate().save(s);

  ka_.save_data_state(s);
  ka_.result_noupdate().save(s);

  kb_.save_data_state(s);
  kb_.result_noupdate().save(s);

  s.put(user_occupations_);
  s.put(tndocc_);
  s.put(nirrep_);
  s.put(ndocc_,nirrep_);
  s.put(osa_);
  s.put(osb_);
  s.put(occa_);
  s.put(occb_);
  s.put(ci1_);
  s.put(ci2_);
}

double
TCSCF::occupation(int ir, int i)
{
  if (i < ndocc_[ir])
    return 2.0;
  else if (ir==osa_ && i==ndocc_[ir])
    return occa_;
  else if (ir==osb_ && i==ndocc_[ir])
    return occb_;
  else
    return 0.0;
}

double
TCSCF::alpha_occupation(int ir, int i)
{
  if (i < ndocc_[ir])
    return 1.0;
  else if (ir==osa_ && i==ndocc_[ir])
    return 0.5*occa_;
  return 0.0;
}

double
TCSCF::beta_occupation(int ir, int i)
{
  if (i < ndocc_[ir])
    return 1.0;
  else if (ir==osb_ && i==ndocc_[ir])
    return 0.5*occb_;
  return 0.0;
}

int
TCSCF::n_fock_matrices() const
{
  return 4;
}

RefSymmSCMatrix
TCSCF::fock(int n)
{
  if (n > 3) {
    ExEnv::err0() << indent
         << "TCSCF::fock: there are only four fock matrices, "
         << scprintf("but fock(%d) was requested\n", n);
    abort();
  }

  if (n==0)
    return focka_.result();
  else if (n==1)
    return fockb_.result();
  else if (n==2)
    return ka_.result();
  else
    return kb_.result();
}

double
TCSCF::magnetic_moment() const
{
  return 0;
}

void
TCSCF::print(ostream&o) const
{
  int i;

  SCF::print(o);

  o << indent << "TCSCF Parameters:\n" << incindent
    << indent << "ndocc = " << tndocc_ << endl
    << indent << scprintf("occa = %f", occa_) << endl
    << indent << scprintf("occb = %f", occb_) << endl
    << indent << scprintf("ci1 = %9.6f", ci1_) << endl
    << indent << scprintf("ci2 = %9.6f", ci2_) << endl
    << indent << "docc = [";
  for (i=0; i < nirrep_; i++)
    o << " " << ndocc_[i];
  o << " ]" << endl
    << indent << "socc = [";
  for (i=0; i < nirrep_; i++)
    o << " " << ((i==osa_ || i==osb_) ? 1 : 0);
  o << " ]" << endl << decindent << endl;
}

//////////////////////////////////////////////////////////////////////////////

void
TCSCF::set_occupations(const RefDiagSCMatrix& ev)
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
                                                 "TCSCF::set_occupations");

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

   if (osa > osb) {
     int tmp=osa;
     osa=osb;
     osb=tmp;
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
        ExEnv::err0() << indent << "TCSCF::set_occupations:  WARNING!!!!\n"
             << incindent << indent
             << scprintf("occupations for irrep %d have changed\n", i+1)
             << indent
             << scprintf("ndocc was %d, changed to %d", ndocc_[i], newdocc[i])
             << endl << decindent;
      }
      if (((osa != osa_ && osa != osb_) || (osb != osb_ && osb != osa_))) {
        ExEnv::err0() << indent << "TCSCF::set_occupations:  WARNING!!!!\n"
             << incindent << indent << "open shell occupations have changed"
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
TCSCF::symmetry_changed()
{
  SCF::symmetry_changed();
  focka_.result_noupdate()=0;
  fockb_.result_noupdate()=0;
  ka_.result_noupdate()=0;
  kb_.result_noupdate()=0;
  nirrep_ = molecule()->point_group()->char_table().ncomp();
  set_occupations(0);
}

//////////////////////////////////////////////////////////////////////////////
//
// scf things
//

void
TCSCF::init_vector()
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
  ao_gmata_ = basis()->matrixkit()->symmmatrix(basis()->basisdim());
  ao_gmata_.assign(0.0);

  ao_gmatb_ = ao_gmata_.clone();
  ao_gmatb_.assign(0.0);

  ao_ka_ = ao_gmata_.clone();
  ao_ka_.assign(0.0);

  ao_kb_ = ao_gmata_.clone();
  ao_kb_.assign(0.0);

  // test to see if we need a guess vector
  if (focka_.result_noupdate() == 0) {
    focka_ = hcore_.clone();
    focka_.result_noupdate().assign(0.0);
    fockb_ = hcore_.clone();
    fockb_.result_noupdate().assign(0.0);
    ka_ = hcore_.clone();
    ka_.result_noupdate().assign(0.0);
    kb_ = hcore_.clone();
    kb_.result_noupdate().assign(0.0);
  }

  // make sure trial vector is set up
  initial_vector();

  oso_scf_vector_ = oso_eigenvectors_.result_noupdate();
}

void
TCSCF::done_vector()
{
  done_threads();

  cl_dens_ = 0;
  cl_dens_diff_ = 0;
  op_densa_ = 0;
  op_densa_diff_ = 0;
  op_densb_ = 0;
  op_densb_diff_ = 0;

  ao_gmata_ = 0;
  ao_gmatb_ = 0;
  ao_ka_ = 0;
  ao_kb_ = 0;

  oso_scf_vector_ = 0;
}

////////////////////////////////////////////////////////////////////////////

RefSymmSCMatrix
TCSCF::density()
{
  if (!density_.computed()) {
    RefSymmSCMatrix dens(so_dimension(), basis_matrixkit());
    RefSymmSCMatrix dens1(so_dimension(), basis_matrixkit());

    so_density(dens, 2.0);
    dens.scale(2.0);

    so_density(dens1, occa_);
    dens1.scale(occa_);
    dens.accumulate(dens1);

    so_density(dens1, occb_);
    dens1.scale(occb_);
    dens.accumulate(dens1);

    dens1=0;

    density_ = dens;
    // only flag the density as computed if the calc is converged
    if (!value_needed()) density_.computed() = 1;
  }

  return density_.result_noupdate();
}

RefSymmSCMatrix
TCSCF::alpha_density()
{
  RefSymmSCMatrix dens(so_dimension(), basis_matrixkit());
  RefSymmSCMatrix dens1(so_dimension(), basis_matrixkit());

  so_density(dens, 2.0);
  so_density(dens1, occa_);
  dens.accumulate(dens1);

  dens.scale(2.0);
  return dens;
}

RefSymmSCMatrix
TCSCF::beta_density()
{
  RefSymmSCMatrix dens(so_dimension(), basis_matrixkit());
  RefSymmSCMatrix dens1(so_dimension(), basis_matrixkit());

  so_density(dens, 2.0);
  so_density(dens1, occb_);
  dens.accumulate(dens1);

  dens.scale(2.0);
  return dens;
}

void
TCSCF::reset_density()
{
  cl_dens_diff_.assign(cl_dens_);

  ao_gmata_.assign(0.0);
  op_densa_diff_.assign(op_densa_);

  ao_gmatb_.assign(0.0);
  op_densb_diff_.assign(op_densb_);

  ao_ka_.assign(0.0);
  ao_kb_.assign(0.0);
}

double
TCSCF::new_density()
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

  so_density(op_densa_, occa_);
  dynamic_cast<BlockedSymmSCMatrix*>(op_densa_.pointer())->block(osb_)->assign(0.0);
  op_densa_.scale(2.0);

  so_density(op_densb_, occb_);
  dynamic_cast<BlockedSymmSCMatrix*>(op_densb_.pointer())->block(osa_)->assign(0.0);
  op_densb_.scale(2.0);

  cl_dens_diff_.accumulate(cl_dens_);
  op_densa_diff_.accumulate(op_densa_);
  op_densb_diff_.accumulate(op_densb_);

  RefSymmSCMatrix del = cl_dens_diff_.copy();
  del.accumulate(op_densa_diff_);
  del.accumulate(op_densb_diff_);

  Ref<SCElementScalarProduct> sp(new SCElementScalarProduct);
  del.element_op(sp.pointer(), del);

  double delta = sp->result();
  delta = sqrt(delta/i_offset(cl_dens_diff_.n()));

  return delta;
}

double
TCSCF::scf_energy()
{
  // first calculate the elements of the CI matrix
  SCFEnergy *eop = new SCFEnergy;
  eop->reference();
  Ref<SCElementOp2> op = eop;

  RefSymmSCMatrix t = focka_.result_noupdate().copy();
  t.accumulate(hcore_);

  RefSymmSCMatrix d = cl_dens_.copy();
  d.accumulate(op_densa_);

  t.element_op(op, d);
  double h11 = eop->result();

  t.assign(fockb_.result_noupdate().copy());
  t.accumulate(hcore_);

  d.assign(cl_dens_);
  d.accumulate(op_densb_);

  eop->reset();
  t.element_op(op, d);
  double h22 = eop->result();

  //t = ka_.result_noupdate();
  //eop->reset();
  //t.element_op(op, op_densb_);
  //double h21 = eop->result();

  t = kb_.result_noupdate();
  eop->reset();
  t.element_op(op, op_densa_);
  double h12 = eop->result();

  op=0;
  eop->dereference();
  delete eop;

  // now diagonalize the CI matrix to get the coefficients
  RefSCDimension l2 = new SCDimension(2);
  Ref<SCMatrixKit> lkit = new LocalSCMatrixKit;
  RefSymmSCMatrix h = lkit->symmmatrix(l2);
  RefSCMatrix hv = lkit->matrix(l2,l2);
  RefDiagSCMatrix hl = lkit->diagmatrix(l2);

  h.set_element(0,0,h11);
  h.set_element(1,1,h22);
  h.set_element(1,0,h12);
  h.diagonalize(hl,hv);

  ci1_ = hv.get_element(0,0);
  ci2_ = hv.get_element(1,0);
  double c1c2 = ci1_*ci2_;

  ExEnv::out0() << indent
               << scprintf("c1 = %10.7f c2 = %10.7f", ci1_, ci2_)
               << endl;

  occa_ = 2*ci1_*ci1_;
  occb_ = 2*ci2_*ci2_;

  double eelec = 0.5*occa_*h11 + 0.5*occb_*h22 + 2.0*c1c2*h12;

  return eelec;
}

Ref<SCExtrapData>
TCSCF::extrap_data()
{
  RefSymmSCMatrix *m = new RefSymmSCMatrix[4];
  m[0] = focka_.result_noupdate();
  m[1] = fockb_.result_noupdate();
  m[2] = ka_.result_noupdate();
  m[3] = kb_.result_noupdate();

  Ref<SCExtrapData> data = new SymmSCMatrixNSCExtrapData(4, m);
  delete[] m;

  return data;
}

RefSymmSCMatrix
TCSCF::effective_fock()
{
  // use fock() instead of cl_fock_ just in case this is called from
  // someplace outside SCF::compute_vector()
  RefSymmSCMatrix mofocka(oso_dimension(), basis_matrixkit());
  mofocka.assign(0.0);

  RefSymmSCMatrix mofockb(oso_dimension(), basis_matrixkit());
  mofockb.assign(0.0);

  RefSymmSCMatrix moka = mofocka.clone();
  moka.assign(0.0);

  RefSymmSCMatrix mokb = mofocka.clone();
  mokb.assign(0.0);

  // use eigenvectors if oso_scf_vector_ is null
  RefSCMatrix vec;
  if (oso_scf_vector_ == 0) {
    vec = eigenvectors();
  } else {
    vec = so_to_orthog_so().t() * oso_scf_vector_;
  }
  mofocka.accumulate_transform(vec, fock(0),
                               SCMatrix::TransposeTransform);
  mofockb.accumulate_transform(vec, fock(1),
                               SCMatrix::TransposeTransform);
  moka.accumulate_transform(vec, fock(2),
                            SCMatrix::TransposeTransform);
  mokb.accumulate_transform(vec, fock(3),
                            SCMatrix::TransposeTransform);

  mofocka.scale(ci1_*ci1_);
  mofockb.scale(ci2_*ci2_);
  moka.scale(ci1_*ci2_);
  mokb.scale(ci1_*ci2_);

  RefSymmSCMatrix mofock = mofocka.copy();
  mofock.accumulate(mofockb);

  BlockedSymmSCMatrix *F = dynamic_cast<BlockedSymmSCMatrix*>(mofock.pointer());
  BlockedSymmSCMatrix *Fa = dynamic_cast<BlockedSymmSCMatrix*>(mofocka.pointer());
  BlockedSymmSCMatrix *Fb = dynamic_cast<BlockedSymmSCMatrix*>(mofockb.pointer());
  BlockedSymmSCMatrix *Ka = dynamic_cast<BlockedSymmSCMatrix*>(moka.pointer());
  BlockedSymmSCMatrix *Kb = dynamic_cast<BlockedSymmSCMatrix*>(mokb.pointer());

  double scalea = (fabs(ci1_) < fabs(ci2_)) ? 1.0/(ci1_*ci1_ + 0.05) : 1.0;
  double scaleb = (fabs(ci2_) < fabs(ci1_)) ? 1.0/(ci2_*ci2_ + 0.05) : 1.0;

  for (int b=0; b < Fa->nblocks(); b++) {
    if (b==osa_) {
      RefSymmSCMatrix f = F->block(b);
      RefSymmSCMatrix fa = Fa->block(b);
      RefSymmSCMatrix fb = Fb->block(b);
      RefSymmSCMatrix kb = Kb->block(b);

      int i,j;

      i=ndocc_[b];
      for (j=0; j < ndocc_[b]; j++)
        f->set_element(i,j,
                       scaleb*(fb->get_element(i,j)-kb->get_element(i,j)));

      j=ndocc_[b];
      for (i=ndocc_[b]+1; i < f->n(); i++)
        f->set_element(i,j,
                       scalea*(fa->get_element(i,j)+kb->get_element(i,j)));

    } else if (b==osb_) {
      RefSymmSCMatrix f = F->block(b);
      RefSymmSCMatrix fa = Fa->block(b);
      RefSymmSCMatrix fb = Fb->block(b);
      RefSymmSCMatrix ka = Ka->block(b);

      int i,j;

      i=ndocc_[b];
      for (j=0; j < ndocc_[b]; j++)
        f->set_element(i,j,
                       scalea*(fa->get_element(i,j)-ka->get_element(i,j)));

      j=ndocc_[b];
      for (i=ndocc_[b]+1; i < f->n(); i++)
        f->set_element(i,j,
                       scaleb*(fb->get_element(i,j)+ka->get_element(i,j)));
    }
  }

  return mofock;
}

/////////////////////////////////////////////////////////////////////////////

void
TCSCF::init_gradient()
{
  // presumably the eigenvectors have already been computed by the time
  // we get here
  oso_scf_vector_ = oso_eigenvectors_.result_noupdate();
}

void
TCSCF::done_gradient()
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
TCSCF::lagrangian()
{
  RefSCMatrix vec = so_to_orthog_so().t() * oso_scf_vector_;

  RefSymmSCMatrix mofocka = focka_.result_noupdate().clone();
  mofocka.assign(0.0);
  mofocka.accumulate_transform(vec, focka_.result_noupdate(),
                               SCMatrix::TransposeTransform);
  mofocka.scale(ci1_*ci1_);

  RefSymmSCMatrix mofockb = mofocka.clone();
  mofockb.assign(0.0);
  mofockb.accumulate_transform(vec, fockb_.result_noupdate(),
                               SCMatrix::TransposeTransform);
  mofockb.scale(ci2_*ci2_);

  // FOa = c1^2*Fa + c1c2*Kb
  RefSymmSCMatrix moka = mofocka.clone();
  moka.assign(0.0);
  moka.accumulate_transform(vec, kb_.result_noupdate(),
                            SCMatrix::TransposeTransform);
  moka.scale(ci1_*ci2_);
  moka.accumulate(mofocka);

  // FOb = c1^2*Fb + c1c2*Ka
  RefSymmSCMatrix mokb = mofocka.clone();
  mokb.assign(0.0);
  mokb.accumulate_transform(vec, ka_.result_noupdate(),
                            SCMatrix::TransposeTransform);
  mokb.scale(ci1_*ci2_);
  mokb.accumulate(mofockb);

  dynamic_cast<BlockedSymmSCMatrix*>(moka.pointer())->block(osb_)->assign(0.0);
  dynamic_cast<BlockedSymmSCMatrix*>(mokb.pointer())->block(osa_)->assign(0.0);

  moka.accumulate(mokb);
  mokb=0;

  // FC = c1^2*Fa + c2^2*Fb
  mofocka.accumulate(mofockb);
  mofockb=0;

  Ref<SCElementOp2> op = new MOLagrangian(this);
  mofocka.element_op(op, moka);
  moka=0;
  mofocka.scale(2.0);

  // transform MO lagrangian to SO basis
  RefSymmSCMatrix so_lag(so_dimension(), basis_matrixkit());
  so_lag.assign(0.0);
  so_lag.accumulate_transform(vec, mofocka);

  // and then from SO to AO
  Ref<PetiteList> pl = integral()->petite_list();
  RefSymmSCMatrix ao_lag = pl->to_AO_basis(so_lag);

  ao_lag.scale(-1.0);

  return ao_lag;
}

RefSymmSCMatrix
TCSCF::gradient_density()
{
  cl_dens_ = basis_matrixkit()->symmmatrix(so_dimension());
  op_densa_ = cl_dens_.clone();
  op_densb_ = cl_dens_.clone();

  so_density(cl_dens_, 2.0);
  cl_dens_.scale(2.0);

  so_density(op_densa_, occa_);
  op_densa_.scale(occa_);

  so_density(op_densb_, occb_);
  op_densb_.scale(occb_);

  dynamic_cast<BlockedSymmSCMatrix*>(op_densa_.pointer())->block(osb_)->assign(0.0);
  dynamic_cast<BlockedSymmSCMatrix*>(op_densb_.pointer())->block(osa_)->assign(0.0);

  Ref<PetiteList> pl = integral()->petite_list(basis());

  cl_dens_ = pl->to_AO_basis(cl_dens_);
  op_densa_ = pl->to_AO_basis(op_densa_);
  op_densb_ = pl->to_AO_basis(op_densb_);

  RefSymmSCMatrix tdens = cl_dens_.copy();
  tdens.accumulate(op_densa_);
  tdens.accumulate(op_densb_);

  op_densa_.scale(2.0/occa_);
  op_densb_.scale(2.0/occb_);

  return tdens;
}

/////////////////////////////////////////////////////////////////////////////

void
TCSCF::init_hessian()
{
}

void
TCSCF::done_hessian()
{
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:

//
// hsosks.cc --- implementation of restricted open shell Kohn-Sham SCF
// derived from clks.cc
//
// Copyright (C) 1997 Limit Point Systems, Inc.
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
#include <util/misc/scexception.h>

#include <math/optimize/scextrapmat.h>

#include <chemistry/qc/basis/petite.h>

#include <chemistry/qc/dft/hsosks.h>
#include <chemistry/qc/scf/lgbuild.h>
#include <chemistry/qc/scf/ltbgrad.h>
#include <chemistry/qc/scf/effh.h>
#include <chemistry/qc/scf/scfops.h>

#include <chemistry/qc/dft/hsoskstmpl.h>

using namespace std;
using namespace sc;

///////////////////////////////////////////////////////////////////////////
// HSOSKS

static ClassDesc HSOSKS_cd(
  typeid(HSOSKS),"HSOSKS",1,"public HSOSSCF",
  0, create<HSOSKS>, create<HSOSKS>);

HSOSKS::HSOSKS(StateIn& s) :
  SavableState(s),
  HSOSSCF(s)
{
  exc_=0;
  integrator_ << SavableState::restore_state(s);
  functional_ << SavableState::restore_state(s);
  vxc_a_ = basis_matrixkit()->symmmatrix(so_dimension());
  vxc_a_.restore(s);
  vxc_b_ = basis_matrixkit()->symmmatrix(so_dimension());
  vxc_b_.restore(s);
}

HSOSKS::HSOSKS(const Ref<KeyVal>& keyval) :
  HSOSSCF(keyval)
{
  exc_=0;
  integrator_ << keyval->describedclassvalue("integrator");
  if (integrator_.null()) integrator_ = new RadialAngularIntegrator();

  functional_ << keyval->describedclassvalue("functional");
  if (functional_.null()) {
    ExEnv::outn() << "ERROR: " << class_name() << ": no \"functional\" given" << endl;
    abort();
  }
}

HSOSKS::~HSOSKS()
{
}

void
HSOSKS::save_data_state(StateOut& s)
{
  HSOSSCF::save_data_state(s);
  SavableState::save_state(integrator_.pointer(),s);
  SavableState::save_state(functional_.pointer(),s);
  vxc_a_.save(s);
  vxc_b_.save(s);
}

int
HSOSKS::value_implemented() const
{
  return 1;
}

bool
HSOSKS::analytic_gradient_implemented() const
{
  return true;
}

void
HSOSKS::print(ostream&o) const
{
  o << indent
    << "High-Spin (Spin-Restricted)  Open-Shell Kohn-Sham (HSOSKS):" << endl;
  o << incindent;
  HSOSSCF::print(o);
  o << indent << "Functional:" << endl;
  o << incindent;
  functional_->print(o);
  o << decindent;
  o << indent << "Integrator:" << endl;
  o << incindent;
  integrator_->print(o);
  o << decindent;
  o << decindent;
}

double
HSOSKS::scf_energy()
{
  double ehf = HSOSSCF::scf_energy();

  return ehf+exc_;
}

RefSymmSCMatrix
HSOSKS::effective_fock()
{
  RefSymmSCMatrix mofock(oso_dimension(), basis_matrixkit());
  mofock.assign(0.0);

  RefSymmSCMatrix mofocko(oso_dimension(), basis_matrixkit());
  mofocko.assign(0.0);

  // use eigenvectors if oso_scf_vector_ is null
  if (oso_scf_vector_.null()) {
    mofock.accumulate_transform(eigenvectors(), fock(0)+cl_vxc(),
                                SCMatrix::TransposeTransform);
    mofocko.accumulate_transform(eigenvectors(), fock(1)+op_vxc(),
                                 SCMatrix::TransposeTransform);
  } else {
    RefSCMatrix so_to_oso_tr = so_to_orthog_so().t();
    mofock.accumulate_transform(so_to_oso_tr * oso_scf_vector_,
                                fock(0)+cl_vxc(),
                                SCMatrix::TransposeTransform);
    mofocko.accumulate_transform(so_to_oso_tr * oso_scf_vector_,
                                 fock(1)+op_vxc(),
                                 SCMatrix::TransposeTransform);
  }

  Ref<SCElementOp2> op = new GSGeneralEffH(this);
  mofock.element_op(op, mofocko);

  return mofock;
}

RefSymmSCMatrix
HSOSKS::lagrangian()
{
  RefSCMatrix so_to_oso_tr = so_to_orthog_so().t();

  RefSymmSCMatrix mofock(oso_dimension(), basis_matrixkit());
  mofock.assign(0.0);
  mofock.accumulate_transform(so_to_oso_tr * oso_scf_vector_,
                              cl_fock_.result_noupdate()+cl_vxc(),
                              SCMatrix::TransposeTransform);

  RefSymmSCMatrix mofocko(oso_dimension(), basis_matrixkit());
  mofocko.assign(0.0);
  mofocko.accumulate_transform(so_to_oso_tr * oso_scf_vector_,
                               op_fock_.result_noupdate()+op_vxc(),
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

Ref<SCExtrapData>
HSOSKS::extrap_data()
{
  Ref<SCExtrapData> data =
    new SymmSCMatrix4SCExtrapData(cl_fock_.result_noupdate(),
                                  op_fock_.result_noupdate(),
                                  vxc_a_, vxc_b_);
  return data;
}

//////////////////////////////////////////////////////////////////////////////

void
HSOSKS::ao_fock(double accuracy)
{
  Ref<PetiteList> pl = integral()->petite_list(basis());
  
  // calculate G.  First transform cl_dens_diff_ to the AO basis, then
  // scale the off-diagonal elements by 2.0
  RefSymmSCMatrix dd = cl_dens_diff_;
  cl_dens_diff_ = pl->to_AO_basis(dd);
  cl_dens_diff_->scale(2.0);
  cl_dens_diff_->scale_diagonal(0.5);

  RefSymmSCMatrix ddo = op_dens_diff_;
  op_dens_diff_ = pl->to_AO_basis(ddo);
  op_dens_diff_->scale(2.0);
  op_dens_diff_->scale_diagonal(0.5);
  
  // now try to figure out the matrix specialization we're dealing with
  // if we're using Local matrices, then there's just one subblock, or
  // see if we can convert G and P to local matrices
  if (local_ || local_dens_) {
    double *gmat, *gmato, *pmat, *pmato;
    
    // grab the data pointers from the G and P matrices
    RefSymmSCMatrix gtmp = get_local_data(cl_gmat_, gmat, SCF::Accum);
    RefSymmSCMatrix ptmp = get_local_data(cl_dens_diff_, pmat, SCF::Read);
    RefSymmSCMatrix gotmp = get_local_data(op_gmat_, gmato, SCF::Accum);
    RefSymmSCMatrix potmp = get_local_data(op_dens_diff_, pmato, SCF::Read);

    signed char * pmax = init_pmax(pmat);
  
//      LocalHSOSKSContribution lclc(gmat, pmat, gmato, pmato, functional_->a0());
//      LocalGBuild<LocalHSOSKSContribution>
//        gb(lclc, tbi_, pl, basis(), scf_grp_, pmax, desired_value_accuracy()/100.0);
//      gb.run();
    int i;
    int nthread = threadgrp_->nthread();
    LocalGBuild<LocalHSOSKSContribution> **gblds =
      new LocalGBuild<LocalHSOSKSContribution>*[nthread];
    LocalHSOSKSContribution **conts = new LocalHSOSKSContribution*[nthread];
    
    double **gmats = new double*[nthread];
    gmats[0] = gmat;
    double **gmatos = new double*[nthread];
    gmatos[0] = gmato;
    
    Ref<GaussianBasisSet> bs = basis();
    int ntri = i_offset(bs->nbasis());

    double gmat_accuracy = accuracy;
    if (min_orthog_res() < 1.0) { gmat_accuracy *= min_orthog_res(); }

    for (i=0; i < nthread; i++) {
      if (i) {
        gmats[i] = new double[ntri];
        memset(gmats[i], 0, sizeof(double)*ntri);
        gmatos[i] = new double[ntri];
        memset(gmatos[i], 0, sizeof(double)*ntri);
      }
      conts[i] = new LocalHSOSKSContribution(gmats[i], pmat, gmatos[i], pmato,
                                             functional_->a0());
      gblds[i] = new LocalGBuild<LocalHSOSKSContribution>(*conts[i], tbis_[i],
        pl, bs, scf_grp_, pmax, gmat_accuracy, nthread, i
        );

      threadgrp_->add_thread(i, gblds[i]);
    }

    Timer tim("start thread");
    if (threadgrp_->start_threads() < 0) {
      ExEnv::err0() << indent
           << "HSOSKS: error starting threads" << endl;
      abort();
    }
    tim.exit("start thread");

    tim.enter("stop thread");
    if (threadgrp_->wait_threads() < 0) {
      ExEnv::err0() << indent
           << "HSOSKS: error waiting for threads" << endl;
      abort();
    }
    tim.exit("stop thread");
      
    double tnint=0;
    for (i=0; i < nthread; i++) {
      tnint += gblds[i]->tnint;

      if (i) {
        for (int j=0; j < ntri; j++) {
          gmat[j] += gmats[i][j];
          gmato[j] += gmatos[i][j];
        }
        delete[] gmats[i];
        delete[] gmatos[i];
      }

      delete gblds[i];
      delete conts[i];
    }

    delete[] gmats;
    delete[] gmatos;
    delete[] gblds;
    delete[] conts;

    delete[] pmax;

    scf_grp_->sum(&tnint, 1, 0, 0);
    ExEnv::out0() << indent << scprintf("%20.0f integrals\n", tnint);
    
    // if we're running on multiple processors, then sum the G matrices
    if (scf_grp_->n() > 1) {
      scf_grp_->sum(gmat, i_offset(basis()->nbasis()));
      scf_grp_->sum(gmato, i_offset(basis()->nbasis()));
    }
    
    // if we're running on multiple processors, or we don't have local
    // matrices, then accumulate gtmp back into G
    if (!local_ || scf_grp_->n() > 1) {
      cl_gmat_->convert_accumulate(gtmp);
      op_gmat_->convert_accumulate(gotmp);
    }
  }

  // for now quit
  else {
    ExEnv::err0() << indent << "Cannot yet use anything but Local matrices\n";
    abort();
  }

  RefSymmSCMatrix dens_a = alpha_ao_density();
  RefSymmSCMatrix dens_b = beta_ao_density();
  integrator_->set_compute_potential_integrals(1);
  integrator_->set_accuracy(accuracy);
  integrator_->integrate(functional_, dens_a, dens_b);
  exc_ = integrator_->value();
  vxc_a_ = dens_a.clone();
  vxc_a_->assign((double*)integrator_->alpha_vmat());
  vxc_a_ = pl->to_SO_basis(vxc_a_);
  vxc_b_ = dens_b.clone();
  vxc_b_->assign((double*)integrator_->beta_vmat());
  vxc_b_ = pl->to_SO_basis(vxc_b_);
  
  // get rid of AO delta P
  cl_dens_diff_ = dd;
  dd = cl_dens_diff_.clone();

  op_dens_diff_ = ddo;
  ddo = op_dens_diff_.clone();

  // now symmetrize the skeleton G matrix, placing the result in dd
  RefSymmSCMatrix skel_gmat = cl_gmat_.copy();
  skel_gmat.scale(1.0/(double)pl->order());
  pl->symmetrize(skel_gmat,dd);

  skel_gmat = op_gmat_.copy();
  skel_gmat.scale(1.0/(double)pl->order());
  pl->symmetrize(skel_gmat,ddo);
  
  // F = H+G
  cl_fock_.result_noupdate().assign(hcore_);
  cl_fock_.result_noupdate().accumulate(dd);

  // Fo = H+G-Go
  op_fock_.result_noupdate().assign(cl_fock_.result_noupdate());
  ddo.scale(-1.0);
  op_fock_.result_noupdate().accumulate(ddo);
  ddo=0;

  dd.assign(0.0);
  accumddh_->accum(dd);
  cl_fock_.result_noupdate().accumulate(dd);
  op_fock_.result_noupdate().accumulate(dd);
  dd=0;

  cl_fock_.computed()=1;
  op_fock_.computed()=1;
}

/////////////////////////////////////////////////////////////////////////////

void
HSOSKS::two_body_energy(double &ec, double &ex)
{
  Timer tim("hsosks e2");
  ec = 0.0;
  ex = 0.0;
  if (local_ || local_dens_) {
    // grab the data pointers from the G and P matrices
    double *dpmat;
    double *spmat;
    tim.enter("local data");
    RefSymmSCMatrix ddens = beta_ao_density();
    RefSymmSCMatrix sdens = alpha_ao_density() - ddens;
    ddens->scale(2.0);
    ddens->accumulate(sdens);
    ddens->scale(2.0);
    ddens->scale_diagonal(0.5);
    sdens->scale(2.0);
    sdens->scale_diagonal(0.5);
    RefSymmSCMatrix dptmp = get_local_data(ddens, dpmat, SCF::Read);
    RefSymmSCMatrix sptmp = get_local_data(sdens, spmat, SCF::Read);
    tim.exit("local data");

    // initialize the two electron integral classes
    Ref<TwoBodyInt> tbi = integral()->electron_repulsion();
    tbi->set_integral_storage(0);

    signed char * pmax = init_pmax(dpmat);
  
    LocalHSOSKSEnergyContribution lclc(dpmat, spmat, functional_->a0());
    Ref<PetiteList> pl = integral()->petite_list();
    LocalGBuild<LocalHSOSKSEnergyContribution>
      gb(lclc, tbi, pl, basis(), scf_grp_, pmax,
         desired_value_accuracy()/100.0);
    gb.run();

    delete[] pmax;

    ec = lclc.ec;
    ex = lclc.ex;
  }
  else {
    ExEnv::err0() << indent << "Cannot yet use anything but Local matrices\n";
    abort();
  }
  tim.exit("hsoshf e2");
}

/////////////////////////////////////////////////////////////////////////////

void
HSOSKS::two_body_deriv(double * tbgrad)
{
  Timer tim("grad");

  int natom3 = 3*molecule()->natom();

  tim.enter("two-body");
  double *hfgrad = new double[natom3];
  memset(hfgrad,0,sizeof(double)*natom3);
  two_body_deriv_hf(hfgrad,functional_->a0());
  //print_natom_3(hfgrad, "Two-body contribution to DFT gradient");
  tim.exit("two-body");

  double *dftgrad = new double[natom3];
  memset(dftgrad,0,sizeof(double)*natom3);
  RefSymmSCMatrix dens_a = alpha_ao_density();
  RefSymmSCMatrix dens_b = beta_ao_density();
  integrator_->init(this);
  integrator_->set_compute_potential_integrals(0);
  integrator_->set_accuracy(desired_gradient_accuracy());
  integrator_->integrate(functional_, dens_a, dens_b, dftgrad);
  // must unset the wavefunction so we don't have a circular list that
  // will not be freed with the reference counting memory manager
  integrator_->done();
  //print_natom_3(dftgrad, "E-X contribution to DFT gradient");

  scf_grp_->sum(dftgrad, natom3);

  for (int i=0; i<natom3; i++) tbgrad[i] += dftgrad[i] + hfgrad[i];
  delete[] dftgrad;
  delete[] hfgrad;

  tim.exit("grad");
}

RefSymmSCMatrix
HSOSKS::cl_vxc()
{
  RefSymmSCMatrix r = vxc_a_+vxc_b_;
  r.scale(0.5);
  return r;
}

RefSymmSCMatrix
HSOSKS::op_vxc()
{
  RefSymmSCMatrix r = vxc_a_.copy();
  return r;
}

/////////////////////////////////////////////////////////////////////////////

void
HSOSKS::init_vector()
{
  integrator_->init(this);
  HSOSSCF::init_vector();
}

void
HSOSKS::done_vector()
{
  integrator_->done();
  HSOSSCF::done_vector();
}

void
HSOSKS::semicanonical()
{
  throw FeatureNotImplemented("Semicanonical orbitals cannot yet be computed for HSOSKS wave function",__FILE__,__LINE__);
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:

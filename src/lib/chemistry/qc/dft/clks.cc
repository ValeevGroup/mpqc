//
// clks.cc --- implementation of the closed shell Kohn-Sham SCF class
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

#include <math/optimize/scextrapmat.h>

#include <chemistry/qc/basis/petite.h>

#include <chemistry/qc/dft/clks.h>
#include <chemistry/qc/scf/lgbuild.h>
#include <chemistry/qc/scf/ltbgrad.h>

#include <chemistry/qc/dft/clkstmpl.h>

using namespace std;
using namespace sc;

///////////////////////////////////////////////////////////////////////////
// CLKS

static ClassDesc CLKS_cd(
  typeid(CLKS),"CLKS",1,"public CLSCF",
  0, create<CLKS>, create<CLKS>);

CLKS::CLKS(StateIn& s) :
  SavableState(s),
  CLSCF(s)
{
  exc_=0;
  integrator_ << SavableState::restore_state(s);
  functional_ << SavableState::restore_state(s);
  vxc_ = basis_matrixkit()->symmmatrix(so_dimension());
  vxc_.restore(s);
}

CLKS::CLKS(const Ref<KeyVal>& keyval) :
  CLSCF(keyval)
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

CLKS::~CLKS()
{
}

void
CLKS::save_data_state(StateOut& s)
{
  CLSCF::save_data_state(s);
  SavableState::save_state(integrator_.pointer(),s);
  SavableState::save_state(functional_.pointer(),s);
  vxc_.save(s);
}

int
CLKS::value_implemented() const
{
  return 1;
}

bool
CLKS::analytic_gradient_implemented() const
{
  return true;
}

void
CLKS::print(ostream&o) const
{
  o << indent << "Closed Shell Kohn-Sham (CLKS):" << endl;
  o << incindent;
  CLSCF::print(o);
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

RefSymmSCMatrix
CLKS::density()
{
  RefSymmSCMatrix dens(so_dimension(), basis_matrixkit());
  so_density(dens, 2.0);
  dens.scale(2.0);
  return dens;
}

double
CLKS::scf_energy()
{
  double ehf = CLSCF::scf_energy();
  return ehf+exc_;
}

RefSymmSCMatrix
CLKS::effective_fock()
{
  RefSymmSCMatrix fa = fock(0) + vxc_;

  RefSymmSCMatrix mofock(oso_dimension(), basis_matrixkit());
  mofock.assign(0.0);

  // use eigenvectors if scf_vector_ is null
  if (oso_scf_vector_.null())
    mofock.accumulate_transform(eigenvectors(), fa,
                                SCMatrix::TransposeTransform);
  else
    mofock.accumulate_transform(so_to_orthog_so().t() * oso_scf_vector_, fa,
                                SCMatrix::TransposeTransform);

  return mofock;
}

Ref<SCExtrapData>
CLKS::initial_extrap_data()
{
  Ref<SCExtrapData> data;
  // If there is an old fock matrix and v_xc around, use that.
  if (cl_fock_.result_noupdate().nonnull()
      && vxc_.nonnull()) {
    data = new SymmSCMatrix2SCExtrapData(cl_fock_.result_noupdate(),
                                         vxc_);
  }
  return data;
}

Ref<SCExtrapData>
CLKS::extrap_data()
{
  Ref<SCExtrapData> data =
    new SymmSCMatrix2SCExtrapData(cl_fock_.result_noupdate(), vxc_);
  return data;
}

//////////////////////////////////////////////////////////////////////////////

void
CLKS::ao_fock(double accuracy)
{
  Ref<PetiteList> pl = integral()->petite_list(basis());
  
  // calculate G.  First transform cl_dens_diff_ to the AO basis, then
  // scale the off-diagonal elements by 2.0
  Timer tim("setup");
  RefSymmSCMatrix dd = cl_dens_diff_;
  cl_dens_diff_ = pl->to_AO_basis(dd);
  cl_dens_diff_->scale(2.0);
  cl_dens_diff_->scale_diagonal(0.5);
  tim.exit("setup");

  // now try to figure out the matrix specialization we're dealing with
  // if we're using Local matrices, then there's just one subblock, or
  // see if we can convert G and P to local matrices

  if (local_ || local_dens_) {
    // grab the data pointers from the G and P matrices
    double *gmat, *pmat;
    tim.enter("local data");
    RefSymmSCMatrix gtmp = get_local_data(cl_gmat_, gmat, SCF::Accum);
    RefSymmSCMatrix ptmp = get_local_data(cl_dens_diff_, pmat, SCF::Read);
    tim.exit("local data");

    tim.enter("init pmax");
    signed char * pmax = init_pmax(pmat);
    tim.exit("init pmax");
  
//      LocalCLKSContribution lclc(gmat, pmat, functional_->a0());
//      LocalGBuild<LocalCLKSContribution>
//        gb(lclc, tbi_, pl, basis(), scf_grp_, pmax, desired_value_accuracy()/100.0);
//      gb.run();
    int i;
    int nthread = threadgrp_->nthread();
    LocalGBuild<LocalCLKSContribution> **gblds =
      new LocalGBuild<LocalCLKSContribution>*[nthread];
    LocalCLKSContribution **conts = new LocalCLKSContribution*[nthread];
    
    double **gmats = new double*[nthread];
    gmats[0] = gmat;
    
    Ref<GaussianBasisSet> bs = basis();
    int ntri = i_offset(bs->nbasis());

    double gmat_accuracy = accuracy;
    if (min_orthog_res() < 1.0) { gmat_accuracy *= min_orthog_res(); }

    for (i=0; i < nthread; i++) {
      if (i) {
        gmats[i] = new double[ntri];
        memset(gmats[i], 0, sizeof(double)*ntri);
      }
      conts[i] = new LocalCLKSContribution(gmats[i], pmat, functional_->a0());
      gblds[i] = new LocalGBuild<LocalCLKSContribution>(*conts[i], tbis_[i],
        pl, bs, scf_grp_, pmax, gmat_accuracy, nthread, i
        );

      threadgrp_->add_thread(i, gblds[i]);
    }

    tim.enter("start thread");
    if (threadgrp_->start_threads() < 0) {
      ExEnv::err0() << indent
           << "CLKS: error starting threads" << endl;
      abort();
    }
    tim.exit("start thread");

    tim.enter("stop thread");
    if (threadgrp_->wait_threads() < 0) {
      ExEnv::err0() << indent
           << "CLKS: error waiting for threads" << endl;
      abort();
    }
    tim.exit("stop thread");

    double tnint=0;
    for (i=0; i < nthread; i++) {
      tnint += gblds[i]->tnint;
      
      if (i) {
        for (int j=0; j < ntri; j++)
          gmat[j] += gmats[i][j];
        delete[] gmats[i];
      }
      delete gblds[i];
      delete conts[i];
    }

    delete[] gmats;
    delete[] gblds;
    delete[] conts;

    delete[] pmax;

    scf_grp_->sum(&tnint, 1, 0, 0);
    ExEnv::out0() << indent << scprintf("%20.0f integrals\n", tnint);
    
    // if we're running on multiple processors, then sum the G matrix
    tim.enter("sum");
    if (scf_grp_->n() > 1)
      scf_grp_->sum(gmat, i_offset(basis()->nbasis()));
    tim.exit("sum");

    // if we're running on multiple processors, or we don't have local
    // matrices, then accumulate gtmp back into G
    tim.enter("accum");
    if (!local_ || scf_grp_->n() > 1)
      cl_gmat_->convert_accumulate(gtmp);
    tim.exit("accum");
  }

  // for now quit
  else {
    ExEnv::out0() << indent << "Cannot yet use anything but Local matrices\n";
    abort();
  }
  
  cl_dens_diff_ = pl->to_AO_basis(cl_dens_);
  cl_dens_diff_.scale(0.5);
  integrator_->set_compute_potential_integrals(1);
  integrator_->set_accuracy(accuracy);
  integrator_->integrate(functional_, cl_dens_diff_, cl_dens_diff_);
  exc_ = integrator_->value();
  RefSymmSCMatrix vxa = cl_gmat_.clone();
  vxa->assign((double*)integrator_->alpha_vmat());
  vxa = pl->to_SO_basis(vxa);
  vxc_ = vxa;

  tim.enter("symm");
  // get rid of AO delta P
  cl_dens_diff_ = dd;
  dd = cl_dens_diff_.clone();

  // now symmetrize the skeleton G matrix, placing the result in dd
  RefSymmSCMatrix skel_gmat = cl_gmat_.copy();
  skel_gmat.scale(1.0/(double)pl->order());
  pl->symmetrize(skel_gmat,dd);
  tim.exit("symm");
  
  
  // F = H+G
  cl_fock_.result_noupdate().assign(hcore_);
  cl_fock_.result_noupdate().accumulate(dd);
  accumddh_->accum(cl_fock_.result_noupdate());
  cl_fock_.computed()=1;
}

/////////////////////////////////////////////////////////////////////////////

void
CLKS::two_body_energy(double &ec, double &ex)
{
  Timer tim("clks e2");
  ec = 0.0;
  ex = 0.0;

  if (local_ || local_dens_) {
    // grab the data pointers from the G and P matrices
    double *pmat;
    tim.enter("local data");
    RefSymmSCMatrix dens = ao_density();
    dens->scale(2.0);
    dens->scale_diagonal(0.5);
    RefSymmSCMatrix ptmp = get_local_data(dens, pmat, SCF::Read);
    tim.exit("local data");

    // initialize the two electron integral classes
    Ref<TwoBodyInt> tbi = integral()->electron_repulsion();
    tbi->set_integral_storage(0);

    tim.enter("init pmax");
    signed char * pmax = init_pmax(pmat);
    tim.exit("init pmax");
  
    LocalCLKSEnergyContribution lclc(pmat, functional_->a0());
    Ref<PetiteList> pl = integral()->petite_list();
    LocalGBuild<LocalCLKSEnergyContribution>
      gb(lclc, tbi, pl, basis(), scf_grp_, pmax,
         desired_value_accuracy()/100.0);
    gb.run();

    delete[] pmax;

    ec = lclc.ec;
    ex = lclc.ex;
  }
  else {
    ExEnv::out0() << indent << "Cannot yet use anything but Local matrices\n";
    abort();
  }
  tim.exit("clks e2");
}

/////////////////////////////////////////////////////////////////////////////

void
CLKS::two_body_deriv(double * tbgrad)
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
  Ref<PetiteList> pl = integral()->petite_list(basis());
  RefSymmSCMatrix aodens = gradient_density();
  aodens.scale(0.5);
  integrator_->set_compute_potential_integrals(0);
  integrator_->init(this);
  integrator_->set_accuracy(desired_gradient_accuracy());
  integrator_->integrate(functional_, aodens, aodens, dftgrad);
  integrator_->done();
  //print_natom_3(dftgrad, "E-X contribution to DFT gradient");

  scf_grp_->sum(dftgrad, natom3);

  for (int i=0; i<natom3; i++) tbgrad[i] += dftgrad[i] + hfgrad[i];
  delete[] dftgrad;
  delete[] hfgrad;

  tim.exit("grad");
}

/////////////////////////////////////////////////////////////////////////////

void
CLKS::init_vector()
{
  integrator_->init(this);
  CLSCF::init_vector();
}

void
CLKS::done_vector()
{
  integrator_->done();
  CLSCF::done_vector();
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:

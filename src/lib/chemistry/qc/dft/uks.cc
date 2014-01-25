//
// uks.cc --- implementation of the unrestricted Hartree-Fock class
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
#include <algorithm>

#include <util/misc/regtime.h>
#include <util/misc/formio.h>
#include <util/state/stateio.h>

#include <math/optimize/scextrapmat.h>

#include <chemistry/qc/basis/petite.h>

#include <chemistry/qc/dft/uks.h>
#include <chemistry/qc/scf/lgbuild.h>
#include <chemistry/qc/scf/ltbgrad.h>

#include <chemistry/qc/dft/ukstmpl.h>

using namespace std;
using namespace sc;

///////////////////////////////////////////////////////////////////////////
// UKS

static ClassDesc UKS_cd(
  typeid(UKS),"UKS",1,"public UnrestrictedSCF",
  0, create<UKS>, create<UKS>);

UKS::UKS(StateIn& s) :
  SavableState(s),
  UnrestrictedSCF(s)
{
  exc_=0;
  integrator_ << SavableState::restore_state(s);
  functional_ << SavableState::restore_state(s);
  vaxc_ = basis_matrixkit()->symmmatrix(so_dimension());
  vaxc_.restore(s);
  vbxc_ = basis_matrixkit()->symmmatrix(so_dimension());
  vbxc_.restore(s);
}

UKS::UKS(const Ref<KeyVal>& keyval) :
  UnrestrictedSCF(keyval)
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

UKS::~UKS()
{
}

void
UKS::save_data_state(StateOut& s)
{
  UnrestrictedSCF::save_data_state(s);
  SavableState::save_state(integrator_.pointer(),s);
  SavableState::save_state(functional_.pointer(),s);
  vaxc_.save(s);
  vbxc_.save(s);
}

int
UKS::value_implemented() const
{
  return 1;
}

bool
UKS::analytic_gradient_implemented() const
{
  return true;
}

double
UKS::scf_energy()
{
  RefSymmSCMatrix mva = vaxc_.copy();
  mva.scale(-1.0);
  focka_.result_noupdate().accumulate(mva);
  RefSymmSCMatrix mvb = vbxc_.copy();
  mvb.scale(-1.0);
  fockb_.result_noupdate().accumulate(mvb);
  double ehf = UnrestrictedSCF::scf_energy();
  focka_.result_noupdate().accumulate(vaxc_);
  fockb_.result_noupdate().accumulate(vbxc_);
  return ehf + exc_;
}

Ref<SCExtrapData>
UKS::extrap_data()
{
  RefSymmSCMatrix *m = new RefSymmSCMatrix[4];
  m[0] = focka_.result_noupdate();
  m[1] = fockb_.result_noupdate();
  m[2] = vaxc_;
  m[3] = vbxc_;
  
  Ref<SCExtrapData> data = new SymmSCMatrixNSCExtrapData(4, m);
  delete[] m;
  
  return data;
}

void
UKS::print(ostream&o) const
{
  o << indent << "(Spin-)Unrestricted Kohn-Sham (UKS):" << endl;
  o << incindent;
  UnrestrictedSCF::print(o);
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

//////////////////////////////////////////////////////////////////////////////

void
UKS::two_body_energy(double &ec, double &ex)
{
  Timer tim("uks e2");
  ec = 0.0;
  ex = 0.0;
  if (local_ || local_dens_) {
    // grab the data pointers from the G and P matrices
    double *apmat;
    double *bpmat;
    tim.enter("local data");
    RefSymmSCMatrix adens = alpha_ao_density();
    RefSymmSCMatrix bdens = beta_ao_density();
    adens->scale(2.0);
    adens->scale_diagonal(0.5);
    bdens->scale(2.0);
    bdens->scale_diagonal(0.5);
    RefSymmSCMatrix aptmp = get_local_data(adens, apmat, SCF::Read);
    RefSymmSCMatrix bptmp = get_local_data(bdens, bpmat, SCF::Read);
    tim.exit("local data");

    // initialize the two electron integral classes
    Ref<TwoBodyInt> tbi = integral()->electron_repulsion();
    tbi->set_integral_storage(0);

    Ref<GaussianBasisSet> bs = basis();
    const int ntri = i_offset(bs->nbasis());

    // use total density to compute bounds
    std::vector<double> pmat(ntri);
    std::transform(apmat, apmat+ntri, bpmat, pmat.begin(),
                   plus<double>());
    signed char * pmax = init_pmax(&(pmat[0]));
  
    LocalUKSEnergyContribution lclc(apmat, bpmat, 0);
    Ref<PetiteList> pl = integral()->petite_list();
    LocalGBuild<LocalUKSEnergyContribution>
      gb(lclc, tbi_, pl, basis(), scf_grp_, pmax, desired_value_accuracy()/100.0);
    gb.run();

    delete[] pmax;

    ec = lclc.ec;
    ex = lclc.ex;
  }
  else {
    ExEnv::out0() << indent << "Cannot yet use anything but Local matrices\n";
    abort();
  }
  tim.exit("uks e2");
}

//////////////////////////////////////////////////////////////////////////////

void
UKS::ao_fock(double accuracy)
{
  Ref<PetiteList> pl = integral()->petite_list(basis());
  
  // calculate G.  First transform diff_densa_ to the AO basis, then
  // scale the off-diagonal elements by 2.0
  RefSymmSCMatrix dda = diff_densa_;
  diff_densa_ = pl->to_AO_basis(dda);
  diff_densa_->scale(2.0);
  diff_densa_->scale_diagonal(0.5);

  RefSymmSCMatrix ddb = diff_densb_;
  diff_densb_ = pl->to_AO_basis(ddb);
  diff_densb_->scale(2.0);
  diff_densb_->scale_diagonal(0.5);

  // now try to figure out the matrix specialization we're dealing with
  // if we're using Local matrices, then there's just one subblock, or
  // see if we can convert G and P to local matrices
  if (local_ || local_dens_) {
    double *gmat, *gmato, *pmat, *pmato;
    
    // grab the data pointers from the G and P matrices
    RefSymmSCMatrix gtmp = get_local_data(gmata_, gmat, SCF::Accum);
    RefSymmSCMatrix ptmp = get_local_data(diff_densa_, pmat, SCF::Read);
    RefSymmSCMatrix gotmp = get_local_data(gmatb_, gmato, SCF::Accum);
    RefSymmSCMatrix potmp = get_local_data(diff_densb_, pmato, SCF::Read);

    Ref<GaussianBasisSet> bs = basis();
    const int ntri = i_offset(bs->nbasis());

    // use total density to compute bounds
    std::vector<double> pmat_total(ntri);
    std::transform(pmat, pmat+ntri, pmato, pmat_total.begin(),
                   plus<double>());
    signed char * pmax = init_pmax(&(pmat_total[0]));
  
//      LocalUKSContribution lclc(gmat, pmat, gmato, pmato, functional_->a0());
//      LocalGBuild<LocalUKSContribution>
//        gb(lclc, tbi_, pl, basis(), scf_grp_, pmax, desired_value_accuracy()/100.0);
//      gb.run();
    int i;
    int nthread = threadgrp_->nthread();
    LocalGBuild<LocalUKSContribution> **gblds =
      new LocalGBuild<LocalUKSContribution>*[nthread];
    LocalUKSContribution **conts = new LocalUKSContribution*[nthread];
    
    double **gmats = new double*[nthread];
    gmats[0] = gmat;
    double **gmatos = new double*[nthread];
    gmatos[0] = gmato;
    
    double gmat_accuracy = accuracy;
    if (min_orthog_res() < 1.0) { gmat_accuracy *= min_orthog_res(); }

    for (i=0; i < nthread; i++) {
      if (i) {
        gmats[i] = new double[ntri];
        memset(gmats[i], 0, sizeof(double)*ntri);
        gmatos[i] = new double[ntri];
        memset(gmatos[i], 0, sizeof(double)*ntri);
      }
      conts[i] = new LocalUKSContribution(gmats[i], pmat, gmatos[i], pmato,
                                          functional_->a0());
      gblds[i] = new LocalGBuild<LocalUKSContribution>(*conts[i], tbis_[i],
        pl, bs, scf_grp_, pmax, gmat_accuracy, nthread, i
        );

      threadgrp_->add_thread(i, gblds[i]);
    }

    Timer tim("start thread");
    if (threadgrp_->start_threads() < 0) {
      ExEnv::err0() << indent
           << "UKS: error starting threads" << endl;
      abort();
    }
    tim.exit("start thread");

    tim.enter("stop thread");
    if (threadgrp_->wait_threads() < 0) {
      ExEnv::err0() << indent
           << "UKS: error waiting for threads" << endl;
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
      gmata_->convert_accumulate(gtmp);
      gmatb_->convert_accumulate(gotmp);
    }
  }

  // for now quit
  else {
    ExEnv::out0() << indent << "Cannot yet use anything but Local matrices\n";
    abort();
  }
  
  diff_densa_ = pl->to_AO_basis(densa_);
  diff_densb_ = pl->to_AO_basis(densb_);
  integrator_->set_compute_potential_integrals(1);
  integrator_->set_accuracy(accuracy);
  integrator_->integrate(functional_, diff_densa_, diff_densb_);
  exc_ = integrator_->value();
  RefSymmSCMatrix vxa = gmata_.clone();
  RefSymmSCMatrix vxb = gmatb_.clone();
  vxa->assign((double*)integrator_->alpha_vmat());
  vxb->assign((double*)integrator_->beta_vmat());
  vxa = pl->to_SO_basis(vxa);
  vxb = pl->to_SO_basis(vxb);
  vaxc_ = vxa;
  vbxc_ = vxb;

  // get rid of AO delta P
  diff_densa_ = dda;
  dda = diff_densa_.clone();

  diff_densb_ = ddb;
  ddb = diff_densb_.clone();

  // now symmetrize the skeleton G matrix, placing the result in dda
  RefSymmSCMatrix skel_gmat = gmata_.copy();
  skel_gmat.scale(1.0/(double)pl->order());
  pl->symmetrize(skel_gmat,dda);

  skel_gmat = gmatb_.copy();
  skel_gmat.scale(1.0/(double)pl->order());
  pl->symmetrize(skel_gmat,ddb);
  
  // Fa = H+Ga
  focka_.result_noupdate().assign(hcore_);
  focka_.result_noupdate().accumulate(dda);
  focka_.result_noupdate().accumulate(vaxc_);

  // Fb = H+Gb
  fockb_.result_noupdate().assign(hcore_);
  fockb_.result_noupdate().accumulate(ddb);
  fockb_.result_noupdate().accumulate(vbxc_);

  dda.assign(0.0);
  accumddh_->accum(dda);
  focka_.result_noupdate().accumulate(dda);
  fockb_.result_noupdate().accumulate(dda);

  focka_.computed()=1;
  fockb_.computed()=1;
}

/////////////////////////////////////////////////////////////////////////////

void
UKS::two_body_deriv(double * tbgrad)
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
  RefSymmSCMatrix ao_dens_a = alpha_ao_density();
  RefSymmSCMatrix ao_dens_b = beta_ao_density();
  integrator_->init(this);
  integrator_->set_compute_potential_integrals(0);
  integrator_->set_accuracy(desired_gradient_accuracy());
  integrator_->integrate(functional_, ao_dens_a, ao_dens_b, dftgrad);
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
UKS::init_vector()
{
  integrator_->init(this);
  UnrestrictedSCF::init_vector();
}

void
UKS::done_vector()
{
  integrator_->done();
  UnrestrictedSCF::done_vector();
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:

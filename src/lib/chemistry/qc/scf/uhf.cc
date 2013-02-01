//
// uhf.cc --- implementation of the unrestricted Hartree-Fock class
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

#include <chemistry/qc/basis/petite.h>

#include <chemistry/qc/scf/uhf.h>
#include <chemistry/qc/scf/lgbuild.h>
#include <chemistry/qc/scf/uhftmpl.h>

using namespace std;
using namespace sc;

///////////////////////////////////////////////////////////////////////////
// UHF

static ClassDesc UHF_cd(
  typeid(UHF),"UHF",1,"public UnrestrictedSCF",
  0, create<UHF>, create<UHF>);

UHF::UHF(StateIn& s) :
  SavableState(s),
  UnrestrictedSCF(s)
{
}

UHF::UHF(const Ref<KeyVal>& keyval) :
  UnrestrictedSCF(keyval)
{
}

UHF::~UHF()
{
}

void
UHF::save_data_state(StateOut& s)
{
  UnrestrictedSCF::save_data_state(s);
}

int
UHF::value_implemented() const
{
  return 1;
}

bool
UHF::analytic_gradient_implemented() const
{
  return true;
}

void
UHF::print(ostream&o) const
{
  o << indent << "(Spin-)Unrestricted Hartree-Fock (UHF):" << endl;
  o << incindent;
  UnrestrictedSCF::print(o);
  o << decindent;
}

//////////////////////////////////////////////////////////////////////////////

void
UHF::two_body_energy(double &ec, double &ex)
{
  Timer tim("uhf e2");
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
  
    LocalUHFEnergyContribution lclc(apmat, bpmat);
    Ref<PetiteList> pl = integral()->petite_list();
    LocalGBuild<LocalUHFEnergyContribution>
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
  tim.exit("uhf e2");
}

//////////////////////////////////////////////////////////////////////////////

void
UHF::ao_fock(double accuracy)
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
  
//      LocalUHFContribution lclc(gmat, pmat, gmato, pmato);
//      LocalGBuild<LocalUHFContribution>
//        gb(lclc, tbi_, pl, basis(), scf_grp_, pmax,
//           desired_value_accuracy()/100.0);
//      gb.run();
    int i;
    int nthread = threadgrp_->nthread();
    LocalGBuild<LocalUHFContribution> **gblds =
      new LocalGBuild<LocalUHFContribution>*[nthread];
    LocalUHFContribution **conts = new LocalUHFContribution*[nthread];
    
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
      conts[i] = new LocalUHFContribution(gmats[i], pmat, gmatos[i], pmato);
      gblds[i] = new LocalGBuild<LocalUHFContribution>(*conts[i], tbis_[i],
        pl, bs, scf_grp_, pmax, gmat_accuracy, nthread, i
        );

      threadgrp_->add_thread(i, gblds[i]);
    }

    Timer tim("start thread");
    if (threadgrp_->start_threads() < 0) {
      ExEnv::err0() << indent
           << "UHF: error starting threads" << endl;
      abort();
    }
    tim.exit("start thread");

    tim.enter("stop thread");
    if (threadgrp_->wait_threads() < 0) {
      ExEnv::err0() << indent
           << "UHF: error waiting for threads" << endl;
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
    ExEnv::err0() << indent << "Cannot yet use anything but Local matrices\n";
    abort();
  }
  
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

  // Fb = H+Gb
  fockb_.result_noupdate().assign(hcore_);
  fockb_.result_noupdate().accumulate(ddb);

  dda.assign(0.0);
  accumddh_->accum(dda);
  focka_.result_noupdate().accumulate(dda);
  fockb_.result_noupdate().accumulate(dda);

  focka_.computed()=1;
  fockb_.computed()=1;
}

/////////////////////////////////////////////////////////////////////////////

void
UHF::two_body_deriv(double * tbgrad)
{
  two_body_deriv_hf(tbgrad, 1.0);
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:

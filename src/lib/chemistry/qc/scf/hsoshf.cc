//
// hsoshf.cc --- implementation of the high-spin open shell Hartree-Fock SCF class
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

#include <chemistry/qc/basis/petite.h>

#include <chemistry/qc/scf/hsoshf.h>
#include <chemistry/qc/scf/lgbuild.h>
#include <chemistry/qc/scf/hsoshftmpl.h>

using namespace std;
using namespace sc;

///////////////////////////////////////////////////////////////////////////
// HSOSHF

static ClassDesc HSOSHF_cd(
  typeid(HSOSHF),"HSOSHF",1,"public HSOSSCF",
  0, create<HSOSHF>, create<HSOSHF>);

HSOSHF::HSOSHF(StateIn& s) :
  SavableState(s),
  HSOSSCF(s)
{
}

HSOSHF::HSOSHF(const Ref<KeyVal>& keyval) :
  HSOSSCF(keyval)
{
}

HSOSHF::~HSOSHF()
{
}

void
HSOSHF::save_data_state(StateOut& s)
{
  HSOSSCF::save_data_state(s);
}

int
HSOSHF::value_implemented() const
{
  return 1;
}

bool
HSOSHF::analytic_gradient_implemented() const
{
  return true;
}

void
HSOSHF::print(ostream&o) const
{
  o << indent
    << "High-Spin (Spin-Restricted) Open-Shell Hartree-Fock (HSOSHF):" << endl;
  o << incindent;
  HSOSSCF::print(o);
  o << decindent;
}

//////////////////////////////////////////////////////////////////////////////

void
HSOSHF::two_body_energy(double &ec, double &ex)
{
  Timer tim("hsoshf e2");
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
  
    LocalHSOSEnergyContribution lclc(dpmat, spmat);
    Ref<PetiteList> pl = integral()->petite_list();
    LocalGBuild<LocalHSOSEnergyContribution>
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

//////////////////////////////////////////////////////////////////////////////

void
HSOSHF::ao_fock(double accuracy)
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
  
    int i;
    int nthread = threadgrp_->nthread();
    LocalGBuild<LocalHSOSContribution> **gblds =
      new LocalGBuild<LocalHSOSContribution>*[nthread];
    LocalHSOSContribution **conts = new LocalHSOSContribution*[nthread];
    
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
      conts[i] = new LocalHSOSContribution(gmats[i], pmat, gmatos[i], pmato);
      gblds[i] = new LocalGBuild<LocalHSOSContribution>(*conts[i], tbis_[i],
        pl, bs, scf_grp_, pmax, gmat_accuracy, nthread, i
        );

      threadgrp_->add_thread(i, gblds[i]);
    }

    Timer tim("start thread");
    if (threadgrp_->start_threads() < 0) {
      ExEnv::err0() << indent
           << "HSOSHF: error starting threads" << endl;
      abort();
    }
    tim.exit("start thread");

    tim.enter("stop thread");
    if (threadgrp_->wait_threads() < 0) {
      ExEnv::err0() << indent
           << "HSOSHF: error waiting for threads" << endl;
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
HSOSHF::two_body_deriv(double * tbgrad)
{
  two_body_deriv_hf(tbgrad, 1.0);
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:

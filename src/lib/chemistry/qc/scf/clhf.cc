//
// clhf.cc --- implementation of the closed shell Hartree-Fock SCF class
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

#include <chemistry/qc/scf/clhf.h>
#include <chemistry/qc/scf/lgbuild.h>
#include <chemistry/qc/scf/clhftmpl.h>

using namespace std;
using namespace sc;

///////////////////////////////////////////////////////////////////////////
// CLHF

static ClassDesc CLHF_cd(
  typeid(CLHF),"CLHF",1,"public CLSCF",
  0, create<CLHF>, create<CLHF>);

CLHF::CLHF(StateIn& s) :
  SavableState(s),
  CLSCF(s)
{
}

CLHF::CLHF(const Ref<KeyVal>& keyval) :
  CLSCF(keyval)
{
}

CLHF::~CLHF()
{
}

void
CLHF::save_data_state(StateOut& s)
{
  CLSCF::save_data_state(s);
}

int
CLHF::value_implemented() const
{
  return 1;
}

bool
CLHF::analytic_gradient_implemented() const
{
  return true;
}

void
CLHF::print(ostream&o) const
{
  o << indent << "Closed Shell Hartree-Fock (CLHF):" << endl;
  o << incindent;
  CLSCF::print(o);
  o << decindent;
}

//////////////////////////////////////////////////////////////////////////////

void
CLHF::ao_fock(double accuracy)
{
  int i;
  int nthread = threadgrp_->nthread();

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

  if (debug_>1) {
    cl_gmat_.print("cl_gmat before build");
    cl_dens_diff_.print("cl_dens_diff before build");
  }

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
  
    tim.enter("ao_gmat");
    LocalGBuild<LocalCLHFContribution> **gblds =
      new LocalGBuild<LocalCLHFContribution>*[nthread];
    LocalCLHFContribution **conts = new LocalCLHFContribution*[nthread];
    
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
      conts[i] = new LocalCLHFContribution(gmats[i], pmat);
      gblds[i] = new LocalGBuild<LocalCLHFContribution>(*conts[i], tbis_[i],
               pl, bs, scf_grp_, pmax, gmat_accuracy, nthread, i
        );

      threadgrp_->add_thread(i, gblds[i]);
    }

    tim.enter("start thread");
    if (threadgrp_->start_threads() < 0) {
      ExEnv::err0() << indent
           << "CLHF: error starting threads" << endl;
      abort();
    }
    tim.exit("start thread");

    tim.enter("stop thread");
    if (threadgrp_->wait_threads() < 0) {
      ExEnv::err0() << indent
           << "CLHF: error waiting for threads" << endl;
      abort();
    }
    tim.exit("stop thread");
      
    double tnint=0;
    double tinttime=0.0;
    for (i=0; i < nthread; i++) {
      tnint += gblds[i]->tnint;
      tinttime += gblds[i]->tinttime;

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
    scf_grp_->sum(&tinttime, 1, 0, 0);
    ExEnv::out0() << indent << scprintf("%20.0f integrals", tnint);
    if (tinttime != 0.0)
      ExEnv::out0() << scprintf("  ( %7.3lf sec)", tinttime);
    ExEnv::out0() << std::endl;


    tim.exit("ao_gmat");

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
    ExEnv::err0() << indent << "Cannot yet use anything but Local matrices\n";
    abort();
  }
  
  tim.enter("symm");
  // get rid of AO delta P
  cl_dens_diff_ = dd;
  dd = cl_dens_diff_.clone();

  // now symmetrize the skeleton G matrix, placing the result in dd
  RefSymmSCMatrix skel_gmat = cl_gmat_.copy();
  skel_gmat.scale(1.0/(double)pl->order());
  if (debug_>1) {
    skel_gmat.print("skel_gmat before symmetrize");
  }
  pl->symmetrize(skel_gmat,dd);
  if (debug_>1) {
    dd.print("dd after symmetrize");
  }
  tim.exit("symm");

  // F = H+G
  cl_fock_.result_noupdate().assign(hcore_);
  cl_fock_.result_noupdate().accumulate(dd);
  accumddh_->accum(cl_fock_.result_noupdate());
  cl_fock_.computed()=1;
}

/////////////////////////////////////////////////////////////////////////////

void
CLHF::two_body_energy(double &ec, double &ex)
{
  Timer tim("clhf e2");
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
  
    LocalCLHFEnergyContribution lclc(pmat);
    Ref<PetiteList> pl = integral()->petite_list();
    LocalGBuild<LocalCLHFEnergyContribution>
      gb(lclc, tbi, pl, basis(), scf_grp_, pmax,
         1.e-20/*desired_value_accuracy()/100.0*/);
    gb.run();

    delete[] pmax;

    ec = lclc.ec;
    ex = lclc.ex;
  }
  else {
    ExEnv::err0() << indent << "Cannot yet use anything but Local matrices\n";
    abort();
  }
  tim.exit("clhf e2");
}

/////////////////////////////////////////////////////////////////////////////

void
CLHF::two_body_deriv(double * tbgrad)
{
  two_body_deriv_hf(tbgrad, 1.0);
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:

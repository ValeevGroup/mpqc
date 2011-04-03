//
// osshf.cc --- implementation of the open shell singlet Hartree-Fock SCF class
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

#include <chemistry/qc/basis/petite.h>

#include <chemistry/qc/scf/osshf.h>
#include <chemistry/qc/scf/lgbuild.h>
#include <chemistry/qc/scf/ltbgrad.h>

#include <chemistry/qc/scf/osshftmpl.h>

using namespace std;
using namespace sc;

///////////////////////////////////////////////////////////////////////////
// OSSHF

static ClassDesc OSSHF_cd(
  typeid(OSSHF),"OSSHF",1,"public OSSSCF",
  0, create<OSSHF>, create<OSSHF>);

OSSHF::OSSHF(StateIn& s) :
  SavableState(s),
  OSSSCF(s)
{
}

OSSHF::OSSHF(const Ref<KeyVal>& keyval) :
  OSSSCF(keyval)
{
}

OSSHF::~OSSHF()
{
}

void
OSSHF::save_data_state(StateOut& s)
{
  OSSSCF::save_data_state(s);
}

int
OSSHF::value_implemented() const
{
  return 1;
}

int
OSSHF::gradient_implemented() const
{
  return 1;
}

void
OSSHF::print(ostream&o) const
{
  OSSSCF::print(o);
}

//////////////////////////////////////////////////////////////////////////////

void
OSSHF::ao_fock(double accuracy)
{
  Ref<PetiteList> pl = integral()->petite_list(basis());
  
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
  
//      LocalOSSContribution lclc(gmat, pmat, gmata, pmata, gmatb, pmatb);
//      LocalGBuild<LocalOSSContribution>
//        gb(lclc, tbi_, pl, basis(), scf_grp_, pmax,
//           desired_value_accuracy()/100.0);
//      gb.run();
    int nthread = threadgrp_->nthread();
    LocalGBuild<LocalOSSContribution> **gblds =
      new LocalGBuild<LocalOSSContribution>*[nthread];
    LocalOSSContribution **conts = new LocalOSSContribution*[nthread];
    
    double **gmatas = new double*[nthread];
    gmatas[0] = gmata;
    double **gmatbs = new double*[nthread];
    gmatbs[0] = gmatb;
    double **gmats = new double*[nthread];
    gmats[0] = gmat;
    
    Ref<GaussianBasisSet> bs = basis();
    int ntri = i_offset(bs->nbasis());

    double gmat_accuracy = accuracy;
    if (min_orthog_res() < 1.0) { gmat_accuracy *= min_orthog_res(); }

    int i;
    for (i=0; i < nthread; i++) {
      if (i) {
        gmatas[i] = new double[ntri];
        memset(gmatas[i], 0, sizeof(double)*ntri);
        gmatbs[i] = new double[ntri];
        memset(gmatbs[i], 0, sizeof(double)*ntri);
        gmats[i] = new double[ntri];
        memset(gmats[i], 0, sizeof(double)*ntri);
      }
      conts[i] = new LocalOSSContribution(gmats[i], pmat,
                                          gmatas[i], pmata, gmatbs[i], pmatb);
      gblds[i] = new LocalGBuild<LocalOSSContribution>(*conts[i], tbis_[i],
        pl, bs, scf_grp_, pmax, gmat_accuracy, nthread, i
        );

      threadgrp_->add_thread(i, gblds[i]);
    }

    Timer tim("start thread");
    if (threadgrp_->start_threads() < 0) {
      ExEnv::err0() << indent
           << "OSSHF: error starting threads" << endl;
      abort();
    }
    tim.exit("start thread");

    tim.enter("stop thread");
    if (threadgrp_->wait_threads() < 0) {
      ExEnv::err0() << indent
           << "OSSHF: error waiting for threads" << endl;
      abort();
    }
    tim.exit("stop thread");
      
    double tnint=0;
    for (i=0; i < nthread; i++) {
      tnint += gblds[i]->tnint;

      if (i) {
        for (int j=0; j < ntri; j++) {
          gmata[j] += gmatas[i][j];
          gmatb[j] += gmatbs[i][j];
          gmat[j]  += gmats[i][j];
        }
        delete[] gmatas[i];
        delete[] gmatbs[i];
        delete[] gmats[i];
      }

      delete gblds[i];
      delete conts[i];
    }

    delete[] gmatas;
    delete[] gmatbs;
    delete[] gmats;
    delete[] gblds;
    delete[] conts;

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
    ExEnv::err0() << indent << "Cannot yet use anything but Local matrices\n";
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
  cl_fock_.result_noupdate().assign(hcore_);
  cl_fock_.result_noupdate().accumulate(dd);

  // Fa = H+G-Ga
  op_focka_.result_noupdate().assign(cl_fock_.result_noupdate());
  dda.scale(-1.0);
  op_focka_.result_noupdate().accumulate(dda);

  // Fb = H+G-Gb
  op_fockb_.result_noupdate().assign(cl_fock_.result_noupdate());
  ddb.scale(-1.0);
  op_fockb_.result_noupdate().accumulate(ddb);

  dd.assign(0.0);
  accumddh_->accum(dd);
  cl_fock_.result_noupdate().accumulate(dd);
  op_focka_.result_noupdate().accumulate(dd);
  op_fockb_.result_noupdate().accumulate(dd);

  cl_fock_.computed()=1;
  op_focka_.computed()=1;
  op_fockb_.computed()=1;
}

//////////////////////////////////////////////////////////////////////////////

void
OSSHF::two_body_energy(double& ec, double& ex)
{
  Timer tim("oshf e2");
  ec = 0.0;
  ex = 0.0;

  if (local_ || local_dens_) {
    Ref<PetiteList> pl = integral()->petite_list(basis());
    
    // grab the data pointers from the G and P matrices
    double *dpmat;
    double *sapmat;
    double *sbpmat;
    tim.enter("local data");
    RefSymmSCMatrix adens = alpha_density();
    RefSymmSCMatrix bdens = beta_density();
    RefSymmSCMatrix ddens = adens+bdens;

    // 2C+a+b - 2(c+b) = a-b
    RefSymmSCMatrix sdensa = bdens.copy();
    sdensa.scale(-2.0);
    sdensa.accumulate(ddens);
    dynamic_cast<BlockedSymmSCMatrix*>(sdensa.pointer())->block(osb_)->assign(0.0);

    // 2C+a+b - 2(c+a) = b-a
    RefSymmSCMatrix sdensb = adens.copy();
    sdensb.scale(-2.0);
    sdensb.accumulate(ddens);
    dynamic_cast<BlockedSymmSCMatrix*>(sdensb.pointer())->block(osa_)->assign(0.0);

    adens=0;
    bdens=0;

    ddens = pl->to_AO_basis(ddens);
    sdensa = pl->to_AO_basis(sdensa);
    sdensb = pl->to_AO_basis(sdensb);
    
    ddens->scale(2.0);
    ddens->scale_diagonal(0.5);
    sdensa->scale(2.0);
    sdensa->scale_diagonal(0.5);
    sdensb->scale(2.0);
    sdensb->scale_diagonal(0.5);

    RefSymmSCMatrix dptmp = get_local_data(ddens, dpmat, SCF::Read);
    RefSymmSCMatrix saptmp = get_local_data(sdensa, sapmat, SCF::Read);
    RefSymmSCMatrix sbptmp = get_local_data(sdensb, sbpmat, SCF::Read);
    tim.exit("local data");

    // initialize the two electron integral classes
    Ref<TwoBodyInt> tbi = integral()->electron_repulsion();
    tbi->set_integral_storage(0);

    signed char * pmax = init_pmax(dpmat);
  
    LocalOSSEnergyContribution lclc(dpmat, sapmat, sbpmat);
    LocalGBuild<LocalOSSEnergyContribution>
      gb(lclc, tbi, pl, basis(), scf_grp_, pmax,
         desired_value_accuracy()/100.0);
    gb.run();

    delete[] pmax;

    ec = lclc.ec;
    ex = lclc.ex;
  }

  // for now quit
  else {
    ExEnv::err0() << indent << "Cannot yet use anything but Local matrices\n";
    abort();
  }
  tim.exit("oshf e2");
}

/////////////////////////////////////////////////////////////////////////////

void
OSSHF::two_body_deriv(double * tbgrad)
{
  Ref<SCElementMaxAbs> m = new SCElementMaxAbs;
  cl_dens_.element_op(m.pointer());
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
    Ref<TwoBodyDerivInt> tbi = integral()->electron_repulsion_deriv();
    Ref<PetiteList> pl = integral()->petite_list();
    LocalTBGrad<LocalOSSGradContribution> tb(l, tbi, pl, basis(), scf_grp_,
                                             tbgrad, pmax, desired_gradient_accuracy());
    tb.run();
    scf_grp_->sum(tbgrad,3 * basis()->molecule()->natom());
  }

  // for now quit
  else {
    ExEnv::err0() << indent
         << "OSSHF::two_body_deriv: can't do gradient yet\n";
    abort();
  }
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:

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

#ifdef __GNUC__
#pragma implementation
#endif

#include <math.h>

#include <util/misc/timer.h>
#include <util/misc/formio.h>

#include <chemistry/qc/basis/petite.h>

#include <chemistry/qc/scf/osshf.h>
#include <chemistry/qc/scf/lgbuild.h>
#include <chemistry/qc/scf/ltbgrad.h>

#include <chemistry/qc/scf/osshftmpl.h>

///////////////////////////////////////////////////////////////////////////
// OSSHF

#define CLASSNAME OSSHF
#define HAVE_STATEIN_CTOR
#define HAVE_KEYVAL_CTOR
#define PARENTS public OSSSCF
#include <util/class/classi.h>
void *
OSSHF::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = OSSSCF::_castdown(cd);
  return do_castdowns(casts,cd);
}

OSSHF::OSSHF(StateIn& s) :
  OSSSCF(s)
  maybe_SavableState(s)
{
}

OSSHF::OSSHF(const RefKeyVal& keyval) :
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
OSSHF::value_implemented()
{
  return 1;
}

int
OSSHF::gradient_implemented()
{
  return 1;
}

int
OSSHF::hessian_implemented()
{
  return 0;
}

void
OSSHF::print(ostream&o)
{
  OSSSCF::print(o);
}

//////////////////////////////////////////////////////////////////////////////

void
OSSHF::ao_fock()
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
  tim_enter("oshf e2");
  ec = 0.0;
  ex = 0.0;

  if (local_ || local_dens_) {
    RefPetiteList pl = integral()->petite_list(basis());
    
    // grab the data pointers from the G and P matrices
    double *dpmat;
    double *sapmat;
    double *sbpmat;
    tim_enter("local data");
    RefSymmSCMatrix adens = alpha_density();
    RefSymmSCMatrix bdens = beta_density();
    RefSymmSCMatrix ddens = adens+bdens;

    // 2C+a+b - 2(c+b) = a-b
    RefSymmSCMatrix sdensa = bdens.copy();
    sdensa.scale(-2.0);
    sdensa.accumulate(ddens);
    BlockedSymmSCMatrix::castdown(sdensa.pointer())->block(osb_)->assign(0.0);

    // 2C+a+b - 2(c+a) = b-a
    RefSymmSCMatrix sdensb = adens.copy();
    sdensb.scale(-2.0);
    sdensb.accumulate(ddens);
    BlockedSymmSCMatrix::castdown(sdensb.pointer())->block(osa_)->assign(0.0);

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
    tim_exit("local data");

    // initialize the two electron integral classes
    tbi_ = integral()->electron_repulsion();
    tbi_->set_integral_storage(0);

    signed char * pmax = init_pmax(dpmat);
  
    LocalOSSEnergyContribution lclc(dpmat, sapmat, sbpmat);
    LocalGBuild<LocalOSSEnergyContribution>
      gb(lclc, tbi_, integral(), basis(), scf_grp_, pmax);
    gb.build_gmat(desired_value_accuracy()/100.0);

    delete[] pmax;

    tbi_ = 0;

    ec = lclc.ec;
    ex = lclc.ex;
  }

  // for now quit
  else {
    cerr << node0 << indent << "Cannot yet use anything but Local matrices\n";
    abort();
  }
  tim_exit("oshf e2");
}

/////////////////////////////////////////////////////////////////////////////

void
OSSHF::two_body_deriv(double * tbgrad)
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
         << "OSSHF::two_body_deriv: can't do gradient yet\n";
    abort();
  }
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:

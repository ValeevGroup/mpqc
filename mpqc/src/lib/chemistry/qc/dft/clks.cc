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

#ifdef __GNUC__
#pragma implementation
#endif

#include <math.h>

#include <util/misc/timer.h>
#include <util/misc/formio.h>
#include <util/state/stateio.h>

#include <math/optimize/scextrapmat.h>

#include <chemistry/qc/basis/petite.h>

#include <chemistry/qc/dft/clks.h>
#include <chemistry/qc/scf/lgbuild.h>
#include <chemistry/qc/scf/ltbgrad.h>

#include <chemistry/qc/dft/clkstmpl.h>

///////////////////////////////////////////////////////////////////////////
// CLKS

#define CLASSNAME CLKS
#define HAVE_STATEIN_CTOR
#define HAVE_KEYVAL_CTOR
#define PARENTS public SCF
#include <util/class/classi.h>
void *
CLKS::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = CLSCF::_castdown(cd);
  return do_castdowns(casts,cd);
}

CLKS::CLKS(StateIn& s) :
  maybe_SavableState(s)
  CLSCF(s)
{
  exc_=0;
  integrator_ = new Murray93Integrator();
  functional_ = new LSDAXFunctional();
}

CLKS::CLKS(const RefKeyVal& keyval) :
  CLSCF(keyval)
{
  exc_=0;
  integrator_ = keyval->describedclassvalue("integrator");
  if (integrator_.null()) integrator_ = new Murray93Integrator();

  functional_ = keyval->describedclassvalue("functional");
  if (functional_.null()) functional_ = new LSDAXFunctional();
}

CLKS::~CLKS()
{
}

void
CLKS::save_data_state(StateOut& s)
{
  CLSCF::save_data_state(s);
}

int
CLKS::value_implemented() const
{
  return 1;
}

int
CLKS::gradient_implemented() const
{
  return 0;
}

void
CLKS::print(ostream&o) const
{
  CLSCF::print(o);
}

RefSymmSCMatrix
CLKS::density()
{
  RefSymmSCMatrix dens(basis_dimension(), basis_matrixkit());
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

  RefSymmSCMatrix mofock = fa.clone();
  mofock.assign(0.0);

  // use eigenvectors if scf_vector_ is null
  if (scf_vector_.null())
    mofock.accumulate_transform(eigenvectors(), fa,
                                SCMatrix::TransposeTransform);
  else
    mofock.accumulate_transform(scf_vector_, fa,
                                SCMatrix::TransposeTransform);

  return mofock;
}

RefSCExtrapData
CLKS::extrap_data()
{
  RefSCExtrapData data =
    new SymmSCMatrix2SCExtrapData(cl_fock_.result_noupdate(), vxc_);
  return data;
}

//////////////////////////////////////////////////////////////////////////////

void
CLKS::ao_fock()
{
  RefPetiteList pl = integral()->petite_list(basis());
  
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
  
    LocalCLKSContribution lclc(gmat, pmat, functional_->a0());
    RefPetiteList pl = integral()->petite_list();
    LocalGBuild<LocalCLKSContribution>
      gb(lclc, tbi_, pl, basis(), scf_grp_, pmax, desired_value_accuracy()/100.0);
    gb.run();

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
  
  tim_enter("integrate");
  cl_dens_diff_ = pl->to_AO_basis(cl_dens_);
  cl_dens_diff_.scale(0.5);
  integrator_->set_wavefunction(this);
  integrator_->set_compute_potential_integrals(1);
  integrator_->integrate(functional_, cl_dens_diff_, cl_dens_diff_);
  exc_ = integrator_->value();
  RefSymmSCMatrix vxa = cl_gmat_.clone();
  vxa->assign((double*)integrator_->alpha_vmat());
  vxa = pl->to_SO_basis(vxa);
  vxc_ = vxa;
  tim_exit("integrate");

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
  cl_fock_.result_noupdate().assign(hcore_);
  cl_fock_.result_noupdate().accumulate(dd);
  accumddh_->accum(cl_fock_.result_noupdate());
  cl_fock_.computed()=1;
}

/////////////////////////////////////////////////////////////////////////////

void
CLKS::two_body_energy(double &ec, double &ex)
{
  tim_enter("clks e2");
  ec = 0.0;
  ex = 0.0;

  if (local_ || local_dens_) {
    // grab the data pointers from the G and P matrices
    double *pmat;
    tim_enter("local data");
    RefSymmSCMatrix dens = ao_density();
    dens->scale(2.0);
    dens->scale_diagonal(0.5);
    RefSymmSCMatrix ptmp = get_local_data(dens, pmat, SCF::Read);
    tim_exit("local data");

    // initialize the two electron integral classes
    tbi_ = integral()->electron_repulsion();
    tbi_->set_integral_storage(0);

    tim_enter("init pmax");
    signed char * pmax = init_pmax(pmat);
    tim_exit("init pmax");
  
    LocalCLKSEnergyContribution lclc(pmat, functional_->a0());
    RefPetiteList pl = integral()->petite_list();
    LocalGBuild<LocalCLKSEnergyContribution>
      gb(lclc, tbi_, pl, basis(), scf_grp_, pmax, desired_value_accuracy()/100.0);
    gb.run();

    delete[] pmax;

    tbi_ = 0;

    ec = lclc.ec;
    ex = lclc.ex;
  }
  else {
    cerr << node0 << indent << "Cannot yet use anything but Local matrices\n";
    abort();
  }
  tim_exit("clks e2");
}

/////////////////////////////////////////////////////////////////////////////

void
CLKS::two_body_deriv(double * tbgrad)
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

    LocalCLKSGradContribution l(pmat);
    RefTwoBodyDerivInt tbi = integral()->electron_repulsion_deriv();
    RefPetiteList pl = integral()->petite_list();
    LocalTBGrad<LocalCLKSGradContribution>
      tb(l, tbi, pl, basis(), scf_grp_, tbgrad, pmax, desired_gradient_accuracy());
    tb.run();
    scf_grp_->sum(tbgrad,3 * basis()->molecule()->natom());
  }

  // for now quit
  else {
    cerr << node0 << indent
         << "CLKS::two_body_deriv: can't do gradient yet\n";
    abort();
  }
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:

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

#ifdef __GNUC__
#pragma implementation
#endif

#include <math.h>

#include <util/misc/timer.h>
#include <util/misc/formio.h>

#include <chemistry/qc/basis/petite.h>

#include <chemistry/qc/scf/clhf.h>
#include <chemistry/qc/scf/lgbuild.h>
#include <chemistry/qc/scf/ltbgrad.h>

///////////////////////////////////////////////////////////////////////////

class LocalCLHFContribution {
  private:
    double * const gmat;
    double * const pmat;

  public:
    LocalCLHFContribution(double *g, double *p) : gmat(g), pmat(p) {}
    ~LocalCLHFContribution() {}

    inline void cont1(int ij, int kl, double val) {
      gmat[ij] += val*pmat[kl];
      gmat[kl] += val*pmat[ij];
    }
    
    inline void cont2(int ij, int kl, double val) {
      val *= -0.25;
      gmat[ij] += val*pmat[kl];
      gmat[kl] += val*pmat[ij];
    }
    
    inline void cont3(int ij, int kl, double val) {
      val *= -0.5;
      gmat[ij] += val*pmat[kl];
      gmat[kl] += val*pmat[ij];
    }
    
    inline void cont4(int ij, int kl, double val) {
      val *= 0.75;
      gmat[ij] += val*pmat[kl];
      gmat[kl] += val*pmat[ij];
    }
    
    inline void cont5(int ij, int kl, double val) {
      val *= 0.5;
      gmat[ij] += val*pmat[kl];
      gmat[kl] += val*pmat[ij];
    }
};

class LocalCLHFEnergyContribution {
  private:
    double * const pmat;

  public:
    double ec;
    double ex;
    
    LocalCLHFEnergyContribution(double *p) : pmat(p) {
      ec=ex=0;
    }
    ~LocalCLHFEnergyContribution() {}

    inline void cont1(int ij, int kl, double val) {
      ec += val*pmat[ij]*pmat[kl];
    }
    
    inline void cont2(int ij, int kl, double val) {
      ex -= 0.25*val*pmat[ij]*pmat[kl];
    }
    
    inline void cont3(int ij, int kl, double val) {
      ex -= 0.5*val*pmat[ij]*pmat[kl];
    }
    
    inline void cont4(int ij, int kl, double val) {
      ec += val*pmat[ij]*pmat[kl];
      ex -= 0.25*val*pmat[ij]*pmat[kl];
    }
    
    inline void cont5(int ij, int kl, double val) {
      ec += val*pmat[ij]*pmat[kl];
      ex -= 0.5*val*pmat[ij]*pmat[kl];
    }
};

class LocalCLHFGradContribution {
  private:
    double * const pmat;

  public:
    LocalCLHFGradContribution(double *p) : pmat(p) {}
    ~LocalCLHFGradContribution() {}

    inline double cont1(int ij, int kl) {
      return pmat[ij]*pmat[kl];
    }

    inline double cont2(int ij, int kl) {
      return pmat[ij]*pmat[kl];
    }
};

#ifdef __GNUC__
template class GBuild<LocalCLHFContribution>;
template class GBuild<LocalCLHFEnergyContribution>;

template class LocalGBuild<LocalCLHFContribution>;
template class LocalGBuild<LocalCLHFEnergyContribution>;

template class TBGrad<LocalCLHFGradContribution>;
template class LocalTBGrad<LocalCLHFGradContribution>;
#endif

///////////////////////////////////////////////////////////////////////////
// CLHF

#define CLASSNAME CLHF
#define HAVE_STATEIN_CTOR
#define HAVE_KEYVAL_CTOR
#define PARENTS public SCF
#include <util/class/classi.h>
void *
CLHF::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = CLSCF::_castdown(cd);
  return do_castdowns(casts,cd);
}

CLHF::CLHF(StateIn& s) :
  CLSCF(s)
  maybe_SavableState(s)
{
}

CLHF::CLHF(const RefKeyVal& keyval) :
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
CLHF::value_implemented()
{
  return 1;
}

int
CLHF::gradient_implemented()
{
  return 1;
}

int
CLHF::hessian_implemented()
{
  return 0;
}

void
CLHF::print(ostream&o)
{
  CLSCF::print(o);
}

//////////////////////////////////////////////////////////////////////////////

void
CLHF::ao_fock()
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
  
    LocalCLHFContribution lclc(gmat, pmat);
    LocalGBuild<LocalCLHFContribution>
      gb(lclc, tbi_, integral(), basis(), scf_grp_, pmax);
    gb.build_gmat(desired_value_accuracy()/100.0);

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
  cl_fock_.result_noupdate().assign(cl_hcore_);
  cl_fock_.result_noupdate().accumulate(dd);
  cl_fock_.computed()=1;
}

/////////////////////////////////////////////////////////////////////////////

void
CLHF::two_body_energy(double &ec, double &ex)
{
  tim_enter("clhf e2");
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
  
    LocalCLHFEnergyContribution lclc(pmat);
    LocalGBuild<LocalCLHFEnergyContribution>
      gb(lclc, tbi_, integral(), basis(), scf_grp_, pmax);
    gb.build_gmat(desired_value_accuracy()/100.0);

    delete[] pmax;

    tbi_ = 0;

    ec = lclc.ec;
    ex = lclc.ex;
  }
  else {
    cerr << node0 << indent << "Cannot yet use anything but Local matrices\n";
    abort();
  }
  tim_exit("clhf e2");
}

/////////////////////////////////////////////////////////////////////////////

void
CLHF::two_body_deriv(double * tbgrad)
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

    LocalCLHFGradContribution l(pmat);
    LocalTBGrad<LocalCLHFGradContribution>
      tb(l, integral(), basis(), scf_grp_);
    tb.build_tbgrad(tbgrad, pmax, desired_gradient_accuracy());
  }

  // for now quit
  else {
    cerr << node0 << indent
         << "CLHF::two_body_deriv: can't do gradient yet\n";
    abort();
  }
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:

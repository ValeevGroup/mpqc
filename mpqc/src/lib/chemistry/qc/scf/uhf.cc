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

#ifdef __GNUC__
#pragma implementation
#endif

#include <math.h>

#include <util/misc/timer.h>
#include <util/misc/formio.h>

#include <chemistry/qc/basis/petite.h>

#include <chemistry/qc/scf/uhf.h>
#include <chemistry/qc/scf/lgbuild.h>
#include <chemistry/qc/scf/ltbgrad.h>

///////////////////////////////////////////////////////////////////////////

class LocalUHFContribution {
  private:
    double * const gmata;
    double * const gmatb;
    double * const pmata;
    double * const pmatb;

  public:
    LocalUHFContribution(double *ga, double *pa, double *gb, double *pb) :
      gmata(ga), pmata(pa), gmatb(gb), pmatb(pb) {}
    ~LocalUHFContribution() {}

    inline void cont1(int ij, int kl, double val) {
      gmata[ij] += val*(pmata[kl]+pmatb[kl]);
      gmata[kl] += val*(pmata[ij]+pmatb[ij]);

      gmatb[ij] += val*(pmata[kl]+pmatb[kl]);
      gmatb[kl] += val*(pmata[ij]+pmatb[ij]);
    }
    
    inline void cont2(int ij, int kl, double val) {
      val *= 0.5;
      gmata[ij] -= val*pmata[kl];
      gmata[kl] -= val*pmata[ij];

      gmatb[ij] -= val*pmatb[kl];
      gmatb[kl] -= val*pmatb[ij];
    }
    
    inline void cont3(int ij, int kl, double val) {
      gmata[ij] -= val*pmata[kl];
      gmata[kl] -= val*pmata[ij];

      gmatb[ij] -= val*pmatb[kl];
      gmatb[kl] -= val*pmatb[ij];
    }
    
    inline void cont4(int ij, int kl, double val) {
      cont1(ij,kl,val);
      cont2(ij,kl,val);
    }
    
    inline void cont5(int ij, int kl, double val) {
      cont1(ij,kl,val);
      cont3(ij,kl,val);
    }
};

class LocalUHFEnergyContribution {
  private:
    double * const pmata;
    double * const pmatb;

  public:
    double ec;
    double ex;
    
    LocalUHFEnergyContribution(double *a, double *b) : pmata(a), pmatb(b) {
      ec=ex=0;
    }

    ~LocalUHFEnergyContribution() {}

    inline void cont1(int ij, int kl, double val) {
      ec += val*(pmata[ij]+pmatb[ij])*(pmata[kl]+pmatb[kl]);
    }
    
    inline void cont2(int ij, int kl, double val) {
      ex -= 0.5*val*(pmata[ij]*pmata[kl]+pmatb[ij]*pmatb[kl]);
    }
    
    inline void cont3(int ij, int kl, double val) {
      ex -= val*(pmata[ij]*pmata[kl]+pmatb[ij]*pmatb[kl]);
    }
    
    inline void cont4(int ij, int kl, double val) {
      cont1(ij,kl,val);
      cont2(ij,kl,val);
    }
    
    inline void cont5(int ij, int kl, double val) {
      cont1(ij,kl,val);
      cont3(ij,kl,val);
    }
};

class LocalUHFGradContribution {
  private:
    double * const pmata;
    double * const pmatb;

  public:
    LocalUHFGradContribution(double *a, double *b) : pmata(a), pmatb(b) {}
    ~LocalUHFGradContribution() {}

    inline double cont1(int ij, int kl) {
      return (pmata[ij]*pmata[kl])+(pmatb[ij]*pmatb[kl]) +
             (pmata[ij]*pmatb[kl])+(pmatb[ij]*pmata[kl]);
    }

    inline double cont2(int ij, int kl) {
      return 2*((pmata[ij]*pmata[kl])+(pmatb[ij]*pmatb[kl]));
    }
};

#ifdef __GNUC__
template class GBuild<LocalUHFContribution>;
template class GBuild<LocalUHFEnergyContribution>;
template class LocalGBuild<LocalUHFContribution>;
template class LocalGBuild<LocalUHFEnergyContribution>;

template class TBGrad<LocalUHFGradContribution>;
template class LocalTBGrad<LocalUHFGradContribution>;
#endif

///////////////////////////////////////////////////////////////////////////
// UHF

#define CLASSNAME UHF
#define HAVE_STATEIN_CTOR
#define HAVE_KEYVAL_CTOR
#define PARENTS public UnrestrictedSCF
#include <util/class/classi.h>
void *
UHF::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = UnrestrictedSCF::_castdown(cd);
  return do_castdowns(casts,cd);
}

UHF::UHF(StateIn& s) :
  UnrestrictedSCF(s)
  maybe_SavableState(s)
{
}

UHF::UHF(const RefKeyVal& keyval) :
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
UHF::value_implemented()
{
  return 1;
}

int
UHF::gradient_implemented()
{
  return 1;
}

int
UHF::hessian_implemented()
{
  return 0;
}

void
UHF::print(ostream&o)
{
  UnrestrictedSCF::print(o);
}

//////////////////////////////////////////////////////////////////////////////

void
UHF::two_body_energy(double &ec, double &ex)
{
  tim_enter("uhf e2");
  ec = 0.0;
  ex = 0.0;
  if (local_ || local_dens_) {
    // grab the data pointers from the G and P matrices
    double *apmat;
    double *bpmat;
    tim_enter("local data");
    RefSymmSCMatrix adens = alpha_ao_density();
    RefSymmSCMatrix bdens = beta_ao_density();
    adens->scale(2.0);
    adens->scale_diagonal(0.5);
    bdens->scale(2.0);
    bdens->scale_diagonal(0.5);
    RefSymmSCMatrix aptmp = get_local_data(adens, apmat, SCF::Read);
    RefSymmSCMatrix bptmp = get_local_data(bdens, bpmat, SCF::Read);
    tim_exit("local data");

    // initialize the two electron integral classes
    tbi_ = integral()->electron_repulsion();
    tbi_->set_integral_storage(0);

    signed char * pmax = init_pmax(apmat);
  
    LocalUHFEnergyContribution lclc(apmat, bpmat);
    LocalGBuild<LocalUHFEnergyContribution>
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
  tim_exit("uhf e2");
}

//////////////////////////////////////////////////////////////////////////////

void
UHF::ao_fock()
{
  RefPetiteList pl = integral()->petite_list(basis());
  
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

    signed char * pmax = init_pmax(pmat);
  
    LocalUHFContribution lclc(gmat, pmat, gmato, pmato);
    LocalGBuild<LocalUHFContribution>
      gb(lclc, tbi_, integral(), basis(), scf_grp_, pmax);
    gb.build_gmat(desired_value_accuracy()/100.0);

    delete[] pmax;

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
    cerr << node0 << indent << "Cannot yet use anything but Local matrices\n";
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

  da.assign(0.0);
  accumddh_->accum(dda);
  focka_.result_noupdate().accumulate(da);
  fockb_.result_noupdate().accumulate(da);

  focka_.computed()=1;
  fockb_.computed()=1;
}

/////////////////////////////////////////////////////////////////////////////

void
UHF::two_body_deriv(double * tbgrad)
{
  RefSCElementMaxAbs m = new SCElementMaxAbs;
  densa_.element_op(m);
  double pmax = m->result();
  m=0;

  // now try to figure out the matrix specialization we're dealing with.
  // if we're using Local matrices, then there's just one subblock, or
  // see if we can convert P to local matrices

  if (local_ || local_dens_) {
    // grab the data pointers from the P matrices
    double *pmata, *pmatb;
    RefSymmSCMatrix ptmpa = get_local_data(densa_, pmata, SCF::Read);
    RefSymmSCMatrix ptmpb = get_local_data(densb_, pmatb, SCF::Read);
  
    LocalUHFGradContribution l(pmata,pmatb);
    LocalTBGrad<LocalUHFGradContribution> tb(l, integral(), basis(),scf_grp_);
    tb.build_tbgrad(tbgrad, pmax, desired_gradient_accuracy());
  }

  // for now quit
  else {
    cerr << node0 << indent
         << "UHF::two_body_deriv: can't do gradient yet\n";
    abort();
  }
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:

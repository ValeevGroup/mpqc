//
// tchf.cc --- implementation of the two-configuration Hartree-Fock SCF class
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

#include <chemistry/qc/scf/tchf.h>
#include <chemistry/qc/scf/lgbuild.h>
#include <chemistry/qc/scf/ltbgrad.h>

///////////////////////////////////////////////////////////////////////////

class LocalTCContribution {
  private:
    double * const gmata;
    double * const gmatb;
    double * const kmata;
    double * const kmatb;

    double * const pmata;
    double * const pmatb;
    double * const opmata;
    double * const opmatb;

    double bound;
  public:
    LocalTCContribution(double *ga, double *pa, double *gb, double *pb,
                        double *ka, double *opa, double *kb, double *opb) :
      gmata(ga), pmata(pa), gmatb(gb), pmatb(pb),
      kmata(ka), opmata(opa), kmatb(kb), opmatb(opb) {}
    ~LocalTCContribution() {}

    void set_bound(double b) { bound = b; }

    inline void cont1(int ij, int kl, double val) {
      gmata[ij] += val*pmata[kl];
      gmata[kl] += val*pmata[ij];

      gmatb[ij] += val*pmatb[kl];
      gmatb[kl] += val*pmatb[ij];
    }
    
    inline void cont2(int ij, int kl, double val) {
      val *= 0.25;
      gmata[ij] -= val*pmata[kl];
      gmata[kl] -= val*pmata[ij];

      gmatb[ij] -= val*pmatb[kl];
      gmatb[kl] -= val*pmatb[ij];

      kmata[ij] += val*opmata[kl];
      kmata[kl] += val*opmata[ij];

      kmatb[ij] += val*opmatb[kl];
      kmatb[kl] += val*opmatb[ij];
    }
    
    inline void cont3(int ij, int kl, double val) {
      val *= 0.5;
      gmata[ij] -= val*pmata[kl];
      gmata[kl] -= val*pmata[ij];

      gmatb[ij] -= val*pmatb[kl];
      gmatb[kl] -= val*pmatb[ij];

      kmata[ij] += val*opmata[kl];
      kmata[kl] += val*opmata[ij];

      kmatb[ij] += val*opmatb[kl];
      kmatb[kl] += val*opmatb[ij];
    }
    
    inline void cont4(int ij, int kl, double val) {
      gmata[ij] += 0.75*val*pmata[kl];
      gmata[kl] += 0.75*val*pmata[ij];

      gmatb[ij] += 0.75*val*pmatb[kl];
      gmatb[kl] += 0.75*val*pmatb[ij];

      kmata[ij] += 0.25*val*opmata[kl];
      kmata[kl] += 0.25*val*opmata[ij];

      kmatb[ij] += 0.25*val*opmatb[kl];
      kmatb[kl] += 0.25*val*opmatb[ij];
    }
    
    inline void cont5(int ij, int kl, double val) {
      val *= 0.5;
      gmata[ij] += val*pmata[kl];
      gmata[kl] += val*pmata[ij];

      gmatb[ij] += val*pmatb[kl];
      gmatb[kl] += val*pmatb[ij];

      kmata[ij] += val*opmata[kl];
      kmata[kl] += val*opmata[ij];

      kmatb[ij] += val*opmatb[kl];
      kmatb[kl] += val*opmatb[ij];
    }
};

class LocalTCEnergyContribution {
  private:
    double * const pmata;
    double * const pmatb;
    double * const opmata;
    double * const opmatb;

  public:
    double eca;
    double exa;
    double ecb;
    double exb;
    double ecab;
    double exab;
    
    LocalTCEnergyContribution(double *pa, double *pb, double *opa, double *opb) :
      pmata(pa), pmatb(pb), opmata(opa), opmatb(opb) {
      exa=eca=0;
      exb=ecb=0;
      exab=ecab=0;
    }
    ~LocalTCEnergyContribution() {}

    inline void cont1(int ij, int kl, double val) {
      eca += val*pmata[ij]*pmata[kl];
      ecb += val*pmatb[ij]*pmatb[kl];
      ecab += val*(pmata[ij]*pmatb[kl]-pmatb[ij]*pmata[kl]);
    }
    
    inline void cont2(int ij, int kl, double val) {
      exa -= 0.25*val*pmata[ij]*pmata[kl];
      exb -= 0.25*val*pmatb[ij]*pmatb[kl];
      exab -= 0.25*val*(pmata[ij]*pmatb[kl]-pmatb[ij]*pmata[kl]);
    }
    
    inline void cont3(int ij, int kl, double val) {
      exa -= 0.5*val*pmata[ij]*pmata[kl];
      exb -= 0.5*val*pmatb[ij]*pmatb[kl];
      exab -= 0.5*val*(pmata[ij]*pmatb[kl]-pmatb[ij]*pmata[kl]);
    }
    
    inline void cont4(int ij, int kl, double val) {
      eca += val*pmata[ij]*pmata[kl];
      ecb += val*pmatb[ij]*pmatb[kl];
      ecab += val*(pmata[ij]*pmatb[kl]-pmatb[ij]*pmata[kl]);

      exa -= 0.25*val*pmata[ij]*pmata[kl];
      exb -= 0.25*val*pmatb[ij]*pmatb[kl];
      exab -= 0.25*val*(pmata[ij]*pmatb[kl]-pmatb[ij]*pmata[kl]);
    }
    
    inline void cont5(int ij, int kl, double val) {
      eca += val*pmata[ij]*pmata[kl];
      ecb += val*pmatb[ij]*pmatb[kl];
      ecab += val*(pmata[ij]*pmatb[kl]-pmatb[ij]*pmata[kl]);

      exa -= 0.5*val*pmata[ij]*pmata[kl];
      exb -= 0.5*val*pmatb[ij]*pmatb[kl];
      exab -= 0.5*val*(pmata[ij]*pmatb[kl]-pmatb[ij]*pmata[kl]);
    }
};

class LocalTCGradContribution {
  private:
    double * const pmat;
    double * const pmata;
    double * const pmatb;
    double c1sq;
    double c2sq;
    double c1c2;

  public:
    LocalTCGradContribution(double *p, double *pa, double *pb,
                            double c1, double c2) :
      pmat(p), pmata(pa), pmatb(pb)
    {
      c1sq = c1*c1;
      c2sq = c2*c2;
      c1c2 = c1*c2;
    }
    ~LocalTCGradContribution() {}

    inline double cont1(int ij, int kl) {
      return pmat[ij]*pmat[kl] +
        c1sq*(pmata[ij]*pmat[kl] + pmat[ij]*pmata[kl]) +
        c2sq*(pmatb[ij]*pmat[kl] + pmat[ij]*pmatb[kl]) +
        0.5*c1sq*pmata[ij]*pmata[kl] +
        0.5*c2sq*pmatb[ij]*pmatb[kl];
    }

    inline double cont2(int ij, int kl) {
      return pmat[ij]*pmat[kl] +
        c1sq*(pmata[ij]*pmat[kl] + pmat[ij]*pmata[kl]) +
        c2sq*(pmatb[ij]*pmat[kl] + pmat[ij]*pmatb[kl]) -
        c1c2*(pmata[ij]*pmatb[kl] + pmatb[ij]*pmata[kl]);
    }
};

///////////////////////////////////////////////////////////////////////////

#ifdef __GNUC__
template class GBuild<LocalTCContribution>;
template class GBuild<LocalTCEnergyContribution>;
template class LocalGBuild<LocalTCContribution>;
template class LocalGBuild<LocalTCEnergyContribution>;

template class TBGrad<LocalTCGradContribution>;
template class LocalTBGrad<LocalTCGradContribution>;
#endif

///////////////////////////////////////////////////////////////////////////
// TCHF

#define CLASSNAME TCHF
#define HAVE_STATEIN_CTOR
#define HAVE_KEYVAL_CTOR
#define PARENTS public TCSCF
#include <util/class/classi.h>
void *
TCHF::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = TCSCF::_castdown(cd);
  return do_castdowns(casts,cd);
}

TCHF::TCHF(StateIn& s) :
  TCSCF(s)
  maybe_SavableState(s)
{
}

TCHF::TCHF(const RefKeyVal& keyval) :
  TCSCF(keyval)
{
}

TCHF::~TCHF()
{
}

void
TCHF::save_data_state(StateOut& s)
{
  TCSCF::save_data_state(s);
}

int
TCHF::value_implemented()
{
  return 1;
}

int
TCHF::gradient_implemented()
{
  return 1;
}

int
TCHF::hessian_implemented()
{
  return 0;
}

void
TCHF::print(ostream&o)
{
  TCSCF::print(o);
}

//////////////////////////////////////////////////////////////////////////////

void
TCHF::ao_fock()
{
  RefPetiteList pl = integral()->petite_list(basis());
  
  // calculate G.  First transform cl_dens_diff_ to the AO basis, then
  // scale the off-diagonal elements by 2.0
  RefSymmSCMatrix da = pl->to_AO_basis(cl_dens_diff_);
  RefSymmSCMatrix db = da.copy();
  RefSymmSCMatrix oda = pl->to_AO_basis(op_densa_diff_);
  RefSymmSCMatrix odb = pl->to_AO_basis(op_densb_diff_);
  da.accumulate(oda);
  db.accumulate(odb);

  da->scale(2.0);
  da->scale_diagonal(0.5);
  
  db->scale(2.0);
  db->scale_diagonal(0.5);
  
  oda->scale(2.0);
  oda->scale_diagonal(0.5);
  
  odb->scale(2.0);
  odb->scale_diagonal(0.5);
  
  // now try to figure out the matrix specialization we're dealing with
  // if we're using Local matrices, then there's just one subblock, or
  // see if we can convert G and P to local matrices
  if (local_ || local_dens_) {

    // grab the data pointers from the G and P matrices
    double *gmata, *gmatb, *kmata, *kmatb, *pmata, *pmatb, *opmata, *opmatb;
    RefSymmSCMatrix gatmp = get_local_data(ao_gmata_, gmata, SCF::Accum);
    RefSymmSCMatrix patmp = get_local_data(da, pmata, SCF::Read);
    RefSymmSCMatrix gbtmp = get_local_data(ao_gmatb_, gmatb, SCF::Accum);
    RefSymmSCMatrix pbtmp = get_local_data(db, pmatb, SCF::Read);
    RefSymmSCMatrix katmp = get_local_data(ao_ka_, kmata, SCF::Accum);
    RefSymmSCMatrix opatmp = get_local_data(oda, opmata, SCF::Read);
    RefSymmSCMatrix kbtmp = get_local_data(ao_kb_, kmatb, SCF::Accum);
    RefSymmSCMatrix opbtmp = get_local_data(odb, opmatb, SCF::Read);
    
    signed char * pmax = init_pmax(pmata);
    signed char * pmaxb = init_pmax(pmatb);
  
    int i;
    for (i=0; i < i_offset(basis()->nshell()); i++) {
      if (pmaxb[i] > pmax[i])
        pmax[i]=pmaxb[i];
    }
    
    delete[] pmaxb;
    
    LocalTCContribution lclc(gmata, pmata, gmatb, pmatb,
                             kmata, opmata, kmatb, opmatb);
    LocalGBuild<LocalTCContribution>
      gb(lclc, tbi_, integral(), basis(), scf_grp_, pmax);
    gb.build_gmat(desired_value_accuracy()/100.0);

    delete[] pmax;

    // if we're running on multiple processors, then sum the G matrices
    if (scf_grp_->n() > 1) {
      scf_grp_->sum(gmata, i_offset(basis()->nbasis()));
      scf_grp_->sum(gmatb, i_offset(basis()->nbasis()));
      scf_grp_->sum(kmata, i_offset(basis()->nbasis()));
      scf_grp_->sum(kmatb, i_offset(basis()->nbasis()));
    }
    
    // if we're running on multiple processors, or we don't have local
    // matrices, then accumulate gtmp back into G
    if (!local_ || scf_grp_->n() > 1) {
      ao_gmata_->convert_accumulate(gatmp);
      ao_gmatb_->convert_accumulate(gbtmp);
      ao_ka_->convert_accumulate(katmp);
      ao_kb_->convert_accumulate(kbtmp);
    }
  }

  // for now quit
  else {
    cerr << node0 << indent << "Cannot yet use anything but Local matrices\n";
    abort();
  }
  
  db=0;
  oda=0;
  odb=0;

  // now symmetrize the skeleton G matrix, placing the result in dd
  RefSymmSCMatrix skel_gmat = ao_gmata_.copy();
  skel_gmat.scale(1.0/(double)pl->order());
  pl->symmetrize(skel_gmat,focka_.result_noupdate());
  
  skel_gmat = ao_gmatb_.copy();
  skel_gmat.scale(1.0/(double)pl->order());
  pl->symmetrize(skel_gmat,fockb_.result_noupdate());
  
  skel_gmat = ao_ka_.copy();
  skel_gmat.scale(1.0/(double)pl->order());
  pl->symmetrize(skel_gmat,ka_.result_noupdate());
  
  skel_gmat = ao_kb_.copy();
  skel_gmat.scale(1.0/(double)pl->order());
  pl->symmetrize(skel_gmat,kb_.result_noupdate());
  
  // Fa = H+Ga
  focka_.result_noupdate().accumulate(hcore_);

  // Fb = H+Gb
  fockb_.result_noupdate().accumulate(hcore_);

  da.assign(0.0);
  accumddh_->accum(da);
  focka_.result_noupdate().accumulate(da);
  fockb_.result_noupdate().accumulate(da);
  ka_.result_noupdate().accumulate(da);
  kb_.result_noupdate().accumulate(da);
  da=0;

  focka_.computed()=1;
  fockb_.computed()=1;
  ka_.computed()=1;
  kb_.computed()=1;
}

/////////////////////////////////////////////////////////////////////////////

void
TCHF::two_body_energy(double &ec, double &ex)
{
  cerr << node0 << indent
       << "TCHF:two_body_energy not implemented"
       << endl;
  abort();
  
  tim_enter("tchf e2");
  ec = 0.0;
  ex = 0.0;

  if (local_ || local_dens_) {
    RefPetiteList pl = integral()->petite_list(basis());

    // grab the data pointers from the G and P matrices
    double *pmata, *pmatb, *spmata, *spmatb;

    tim_enter("local data");
    RefSymmSCMatrix densa = alpha_density();
    RefSymmSCMatrix densb = beta_density();
    RefSymmSCMatrix densc = densb.clone();
    so_density(densc, 2.0);
    densc.scale(-2.0);

    RefSymmSCMatrix sdensa = densa.copy();
    sdensa.accumulate(densc);
    
    RefSymmSCMatrix sdensb = densb.copy();
    sdensb.accumulate(densc);

    densc=0;
    
    densa = pl->to_AO_basis(densa);
    densb = pl->to_AO_basis(densb);
    sdensa = pl->to_AO_basis(sdensa);
    sdensb = pl->to_AO_basis(sdensb);
    
    densa->scale(2.0);
    densa->scale_diagonal(0.5);
    densb->scale(2.0);
    densb->scale_diagonal(0.5);
    sdensa->scale(2.0);
    sdensa->scale_diagonal(0.5);
    sdensb->scale(2.0);
    sdensb->scale_diagonal(0.5);

    RefSymmSCMatrix ptmpa = get_local_data(densa, pmata, SCF::Read);
    RefSymmSCMatrix ptmpb = get_local_data(densb, pmatb, SCF::Read);
    RefSymmSCMatrix sptmpa = get_local_data(sdensa, spmata, SCF::Read);
    RefSymmSCMatrix sptmpb = get_local_data(sdensb, spmatb, SCF::Read);
    tim_exit("local data");

    // initialize the two electron integral classes
    tbi_ = integral()->electron_repulsion();
    tbi_->set_integral_storage(0);

    tim_enter("init pmax");
    signed char * pmax = init_pmax(pmata);
    tim_exit("init pmax");
  
    LocalTCEnergyContribution lclc(pmata,pmatb,spmata,spmatb);
    LocalGBuild<LocalTCEnergyContribution>
      gb(lclc, tbi_, integral(), basis(), scf_grp_, pmax);
    gb.build_gmat(desired_value_accuracy()/100.0);

    delete[] pmax;

    tbi_ = 0;

    printf("%20.10f %20.10f\n", lclc.eca, lclc.exa);
    printf("%20.10f %20.10f\n", lclc.ecb, lclc.exb);
    printf("%20.10f %20.10f\n", lclc.ecab, lclc.exab);
    
  }
  else {
    cerr << node0 << indent << "Cannot yet use anything but Local matrices\n";
    abort();
  }
  tim_exit("tchf e2");
}

/////////////////////////////////////////////////////////////////////////////

void
TCHF::two_body_deriv(double * tbgrad)
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
  
    LocalTCGradContribution l(pmat,pmata,pmatb,ci1_,ci2_);
    LocalTBGrad<LocalTCGradContribution> tb(l, integral(), basis(), scf_grp_);
    tb.build_tbgrad(tbgrad, pmax, desired_gradient_accuracy());
  }

  // for now quit
  else {
    cerr << node0 << indent
         << "TCHF::two_body_deriv: can't do gradient yet\n";
    abort();
  }
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:

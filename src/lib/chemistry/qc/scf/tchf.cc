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

#include <math.h>

#include <util/misc/regtime.h>
#include <util/misc/formio.h>
#include <util/state/stateio.h>

#include <chemistry/qc/basis/petite.h>

#include <chemistry/qc/scf/tchf.h>
#include <chemistry/qc/scf/lgbuild.h>
#include <chemistry/qc/scf/ltbgrad.h>

#include <chemistry/qc/scf/tchftmpl.h>

using namespace std;
using namespace sc;

///////////////////////////////////////////////////////////////////////////
// TCHF

static ClassDesc TCHF_cd(
  typeid(TCHF),"TCHF",1,"public TCSCF",
  0, create<TCHF>, create<TCHF>);

TCHF::TCHF(StateIn& s) :
  SavableState(s),
  TCSCF(s)
{
}

TCHF::TCHF(const Ref<KeyVal>& keyval) :
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
TCHF::value_implemented() const
{
  return 1;
}

bool
TCHF::analytic_gradient_implemented() const
{
  return true;
}

void
TCHF::print(ostream&o) const
{
  TCSCF::print(o);
}

//////////////////////////////////////////////////////////////////////////////

void
TCHF::ao_fock(double accuracy)
{
  Ref<PetiteList> pl = integral()->petite_list(basis());
  
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
    
//      LocalTCContribution lclc(gmata, pmata, gmatb, pmatb,
//                               kmata, opmata, kmatb, opmatb);
//      LocalGBuild<LocalTCContribution>
//        gb(lclc, tbi_, pl, basis(), scf_grp_, pmax,
//           desired_value_accuracy()/100.0);
//      gb.run();
    int nthread = threadgrp_->nthread();
    LocalGBuild<LocalTCContribution> **gblds =
      new LocalGBuild<LocalTCContribution>*[nthread];
    LocalTCContribution **conts = new LocalTCContribution*[nthread];
    
    double **gmatas = new double*[nthread];
    gmatas[0] = gmata;
    double **gmatbs = new double*[nthread];
    gmatbs[0] = gmatb;
    double **kmatas = new double*[nthread];
    kmatas[0] = kmata;
    double **kmatbs = new double*[nthread];
    kmatbs[0] = kmatb;
    
    Ref<GaussianBasisSet> bs = basis();
    int ntri = i_offset(bs->nbasis());

    double gmat_accuracy = accuracy;
    if (min_orthog_res() < 1.0) { gmat_accuracy *= min_orthog_res(); }

    for (i=0; i < nthread; i++) {
      if (i) {
        gmatas[i] = new double[ntri];
        memset(gmatas[i], 0, sizeof(double)*ntri);
        gmatbs[i] = new double[ntri];
        memset(gmatbs[i], 0, sizeof(double)*ntri);
        kmatas[i] = new double[ntri];
        memset(kmatas[i], 0, sizeof(double)*ntri);
        kmatbs[i] = new double[ntri];
        memset(kmatbs[i], 0, sizeof(double)*ntri);
      }
      conts[i] = new LocalTCContribution(gmatas[i], pmata, gmatbs[i], pmatb,
                                         kmatas[i], opmata, kmatbs[i], opmatb);
      gblds[i] = new LocalGBuild<LocalTCContribution>(*conts[i], tbis_[i],
        pl, bs, scf_grp_, pmax, gmat_accuracy, nthread, i
        );

      threadgrp_->add_thread(i, gblds[i]);
    }

    Timer tim("start thread");
    if (threadgrp_->start_threads() < 0) {
      ExEnv::err0() << indent
           << "TCHF: error starting threads" << endl;
      abort();
    }
    tim.exit("start thread");

    tim.enter("stop thread");
    if (threadgrp_->wait_threads() < 0) {
      ExEnv::err0() << indent
           << "TCHF: error waiting for threads" << endl;
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
          kmata[j] += kmatas[i][j];
          kmatb[j] += kmatbs[i][j];
        }
        delete[] gmatas[i];
        delete[] gmatbs[i];
        delete[] kmatas[i];
        delete[] kmatbs[i];
      }

      delete gblds[i];
      delete conts[i];
    }

    delete[] gmatas;
    delete[] gmatbs;
    delete[] kmatas;
    delete[] kmatbs;
    delete[] gblds;
    delete[] conts;

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
    ExEnv::err0() << indent << "Cannot yet use anything but Local matrices\n";
    abort();
  }
  
  da=0;
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

  RefSymmSCMatrix ddh = hcore_.clone();
  ddh.assign(0.0);
  accumddh_->accum(ddh);
  focka_.result_noupdate().accumulate(ddh);
  fockb_.result_noupdate().accumulate(ddh);
  ka_.result_noupdate().accumulate(ddh);
  kb_.result_noupdate().accumulate(ddh);
  ddh=0;

  focka_.computed()=1;
  fockb_.computed()=1;
  ka_.computed()=1;
  kb_.computed()=1;
}

/////////////////////////////////////////////////////////////////////////////

void
TCHF::two_body_energy(double &ec, double &ex)
{
  ExEnv::err0() << indent
       << "TCHF:two_body_energy not implemented"
       << endl;
  abort();
  
  Timer tim("tchf e2");
  ec = 0.0;
  ex = 0.0;

  if (local_ || local_dens_) {
    Ref<PetiteList> pl = integral()->petite_list(basis());

    // grab the data pointers from the G and P matrices
    double *pmata, *pmatb, *spmata, *spmatb;

    tim.enter("local data");
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
    tim.exit("local data");

    // initialize the two electron integral classes
    Ref<TwoBodyInt> tbi = integral()->electron_repulsion();
    tbi->set_integral_storage(0);

    tim.enter("init pmax");
    signed char * pmax = init_pmax(pmata);
    tim.exit("init pmax");
  
    LocalTCEnergyContribution lclc(pmata,pmatb,spmata,spmatb);
    LocalGBuild<LocalTCEnergyContribution>
      gb(lclc, tbi, pl, basis(), scf_grp_, pmax,
         desired_value_accuracy()/100.0);
    gb.run();

    delete[] pmax;

    printf("%20.10f %20.10f\n", lclc.eca, lclc.exa);
    printf("%20.10f %20.10f\n", lclc.ecb, lclc.exb);
    printf("%20.10f %20.10f\n", lclc.ecab, lclc.exab);
    
  }
  else {
    ExEnv::err0() << indent << "Cannot yet use anything but Local matrices\n";
    abort();
  }
  tim.exit("tchf e2");
}

/////////////////////////////////////////////////////////////////////////////

void
TCHF::two_body_deriv(double * tbgrad)
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
  
    LocalTCGradContribution l(pmat,pmata,pmatb,ci1_,ci2_);
    Ref<TwoBodyDerivInt> tbi = integral()->electron_repulsion_deriv();
    Ref<PetiteList> pl = integral()->petite_list();
    LocalTBGrad<LocalTCGradContribution> tb(l, tbi, pl, basis(), scf_grp_,
                                            tbgrad, pmax, desired_gradient_accuracy());
    tb.run();
    scf_grp_->sum(tbgrad,3 * basis()->molecule()->natom());
  }

  // for now quit
  else {
    ExEnv::err0() << indent
         << "TCHF::two_body_deriv: can't do gradient yet\n";
    abort();
  }
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:

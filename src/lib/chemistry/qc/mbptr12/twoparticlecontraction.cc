//
// twoparticlecontraction.cc
//
// Copyright (C) 2005 Edward Valeev
//
// Author: Edward Valeev <evaleev@vt.edu>
// Maintainer: EV
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

#include <util/misc/formio.h>
#include <util/class/scexception.h>
#include <util/ref/ref.h>
#include <chemistry/qc/mbptr12/twoparticlecontraction.h>
#include <math/scmat/blas.h>

using namespace sc;
using namespace sc::mbptr12;

////////
// Two ParticleContraction

static ClassDesc TwoParticleContraction_cd(
  typeid(TwoParticleContraction),"TwoParticleContraction",1,"virtual public SavableState",
  0, 0, 0);

TwoParticleContraction::TwoParticleContraction(unsigned int nrow, unsigned int ncol) :
  nrow_(nrow), ncol_(ncol)
{
}

TwoParticleContraction::TwoParticleContraction(StateIn& si) : SavableState(si)
{
  si.get(nrow_);
  si.get(ncol_);
}

void
TwoParticleContraction::save_data_state(StateOut& so)
{
  so.put(nrow_);
  so.put(ncol_);
}

unsigned int
TwoParticleContraction::nrow() const { return nrow_; }

unsigned int
TwoParticleContraction::ncol() const { return ncol_; }

double
TwoParticleContraction::dot_prod(const double* A, const double* B) const
{
  const blasint blksize = nrow_ * ncol_;
  const blasint unitstride = 1;
  return F77_DDOT(&blksize,A,&unitstride,B,&unitstride);
}


//////////////
// Direct_Contraction

static ClassDesc Direct_Contraction_cd(
  typeid(Direct_Contraction),"Direct_Contraction",1,"public TwoParticleContraction",
  0, 0, create<Direct_Contraction>);

Direct_Contraction::Direct_Contraction(unsigned int nrow, unsigned int ncol, double scale) :
  TwoParticleContraction(nrow,ncol), scale_(scale)
{
}

Direct_Contraction::Direct_Contraction(StateIn& si) : TwoParticleContraction(si)
{
  si.get(scale_);
}

void
Direct_Contraction::save_data_state(StateOut& so)
{
  so.put(scale_);
}

double
Direct_Contraction::contract(const double* A, const double* B) const
{
  return scale_ * dot_prod(A,B);
}

//////////////
// ABS_OBS_Contraction

static ClassDesc ABS_OBS_Contraction_cd(
  typeid(ABS_OBS_Contraction),"ABS_OBS_Contraction",1,"public TwoParticleContraction",
  0, 0, create<ABS_OBS_Contraction>);

ABS_OBS_Contraction::ABS_OBS_Contraction(unsigned int nobs, unsigned int nocc1, unsigned int nocc2) :
  TwoParticleContraction(nobs,nobs), nocc1_(nocc1), nocc2_(nocc2)
{
}

ABS_OBS_Contraction::ABS_OBS_Contraction(StateIn& si) : TwoParticleContraction(si)
{
  si.get(nocc1_);
  si.get(nocc2_);
}

void
ABS_OBS_Contraction::save_data_state(StateOut& so)
{
  so.put(nocc1_);
  so.put(nocc2_);
}

double
ABS_OBS_Contraction::contract(const double* A, const double* B) const
{
  const blasint unitstride = 1;
  const blasint nobs = nrow();
  const blasint nvir1 = nobs - nocc1_;
  const blasint nvir2 = nobs - nocc2_;
  const blasint nocc2 = nocc2_;
  
  const double* Aoff = A;
  const double* Boff = B;
  double result = 0.0;
  /// contracting occ-occ blocks of A and B.
  for(blasint o=0; o<nocc1_; o++,Aoff+=nobs,Boff+=nobs)
    result += F77_DDOT(&nocc2,Aoff,&unitstride,Boff,&unitstride);

  /// contracting vir-vir blocks of A and B
  Aoff = A + nocc1_*nobs + nocc2_;  Boff = B + nocc1_*nobs + nocc2_;
  for(blasint v=0; v<nvir1; v++,Aoff+=nobs,Boff+=nobs)
    result -= F77_DDOT(&nvir2,Aoff,&unitstride,Boff,&unitstride);
  
  return result;
}

//////////////
// CABS_OBS_Contraction

static ClassDesc CABS_OBS_Contraction_cd(
  typeid(CABS_OBS_Contraction),"CABS_OBS_Contraction",1,"public TwoParticleContraction",
  0, 0, create<CABS_OBS_Contraction>);

CABS_OBS_Contraction::CABS_OBS_Contraction(unsigned int nobs) :
  TwoParticleContraction(nobs,nobs)
{
}

CABS_OBS_Contraction::CABS_OBS_Contraction(StateIn& si) : TwoParticleContraction(si)
{
}

void
CABS_OBS_Contraction::save_data_state(StateOut& so)
{
}

double
CABS_OBS_Contraction::contract(const double* A, const double* B) const
{
  return -1.0 * dot_prod(A,B);
}


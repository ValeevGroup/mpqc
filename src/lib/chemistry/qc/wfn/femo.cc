//
// femo.cc
//
// Copyright (C) 2008 Edward Valeev
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

#include <chemistry/qc/wfn/femo.h>
#include <math/scmat/blocked.h>
#include <util/misc/scexception.h>

using namespace sc;

FEMO::FEMO(int nalpha, int nbeta,
           const RefDiagSCMatrix& evalsa, const RefDiagSCMatrix& evalsb) :
        E_(0.0), nalpha_(nalpha), nbeta_(nbeta)
{
  // first convert evals to something we can deal with easily
  BlockedDiagSCMatrix *bevalsa = require_dynamic_cast<BlockedDiagSCMatrix*>(evalsa, "FEMO::FEMO");
  BlockedDiagSCMatrix *bevalsb = evalsb.null() ? bevalsa : require_dynamic_cast<BlockedDiagSCMatrix*>(evalsb, "FEMO::FEMO");

  const int nirrep_ = bevalsa->nblocks();
  if (nirrep_ != bevalsb->nblocks())
    throw ProgrammingError("FEMO::FEMO -- alpha and beta eigenvalues have different number of symmetry blocks",__FILE__,__LINE__);

  const RefSCDimension modim = bevalsa->dim();
  if (modim.n() != bevalsb->dim().n())
    throw ProgrammingError("FEMO::FEMO -- different number of alpha and beta eigenvalues given",__FILE__,__LINE__);
  
  double **evalsa_ = new double*[nirrep_];
  double **evalsb_ = new double*[nirrep_];
  for (int i=0; i < nirrep_; i++) {
    const int nf = modim->blocks()->size(i);
    if (nf) {
      evalsa_[i] = new double[nf];
      evalsb_[i] = new double[nf];
      bevalsa->block(i)->convert(evalsa_[i]);
      bevalsb->block(i)->convert(evalsb_[i]);
    } else {
      evalsa_[i] = 0;
      evalsb_[i] = 0;
    }
  }

  // Populate alpha orbitals
  alphapi_.resize(nirrep_);
  for(int g=0; g<nirrep_; ++g) alphapi_[g] = 0;

  for (int i=0; i < nalpha_; i++) {
    // find lowest eigenvalue
    int lir=0,ln=0;
    double lowest=999999999;

    for (int ir=0; ir < nirrep_; ir++) {
      int nf=modim->blocks()->size(ir);
      if (!nf)
        continue;
      for (int j=0; j < nf; j++) {
        if (evalsa_[ir][j] < lowest) {
          lowest=evalsa_[ir][j];
          lir=ir;
          ln=j;
        }
      }
    }
    alphapi_[lir]++;
    E_ += evalsa_[lir][ln];
    evalsa_[lir][ln]=999999999;
  }

  // populate the beta orbitals
  betapi_.resize(nirrep_);
  for(int g=0; g<nirrep_; ++g) betapi_[g] = 0;

  for (int i=0; i < nbeta_; i++) {
    // find lowest eigenvalue
    int lir=0,ln=0;
    double lowest=999999999;

    for (int ir=0; ir < nirrep_; ir++) {
      int nf=modim->blocks()->size(ir);
      if (!nf)
        continue;
      for (int j=0; j < nf; j++) {
        if (evalsb_[ir][j] < lowest) {
          lowest=evalsb_[ir][j];
          lir=ir;
          ln=j;
        }
      }
    }
    betapi_[lir]++;
    E_ += evalsb_[lir][ln];
    evalsb_[lir][ln]=999999999;
  }

  // get rid of evals_
  for (int i=0; i < nirrep_; i++) {
    if (evalsa_[i])
      delete[] evalsa_[i];
    if (evalsb_[i])
      delete[] evalsb_[i];
  }
  delete[] evalsa_;
  delete[] evalsb_;
}

double
FEMO::E() const { return E_; }

int
FEMO::nalpha() const { return nalpha_; }

int
FEMO::nbeta() const { return nbeta_; }

int
FEMO::nalpha(int irrep) const { return alphapi_.at(irrep); }

int
FEMO::nbeta(int irrep) const { return betapi_.at(irrep); }

////

HundsFEMOSeeker::HundsFEMOSeeker(int nelectron, double Etol, bool allow_closed_shell,
                                 const RefDiagSCMatrix& evalsa, const RefDiagSCMatrix& evalsb)
{
  if (Etol < 0.)
    Etol = 0.0;
  
  // start with the lowest multiplicity until the energy changes by more than Etol
  int nalpha = (nelectron+1)/2;
  int nbeta = nelectron/2;
  // is this configuration allowed?
  if (!allow_closed_shell && nalpha == nbeta) {
    ++nalpha;
    --nbeta;
  }
  Ref<FEMO> femo_prev;
  Ref<FEMO> femo_current = new FEMO(nalpha,nbeta,evalsa,evalsb);
  double deltaE = Etol - 1.0;
  while (deltaE < Etol && nbeta >= 1) {
    femo_prev = femo_current;
    ++nalpha;
    --nbeta;
    femo_current = new FEMO(nalpha,nbeta,evalsa,evalsb);
    deltaE = femo_current->E() - femo_prev->E();
  }
  // If we ran out of electrons, grab the last configuration
  if (nbeta == 0)
    result_ = femo_current;
  // else grab the one before
  else
    result_ = femo_prev;
}

double
HundsFEMOSeeker::tolerance(0.001);

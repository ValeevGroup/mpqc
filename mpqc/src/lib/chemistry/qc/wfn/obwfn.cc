//
// obwfn.cc
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@limitpt.com>
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

#include <math/scmat/local.h>
#include <math/scmat/blocked.h>

#include <chemistry/qc/basis/integral.h>
#include <chemistry/qc/basis/obint.h>
#include <chemistry/qc/basis/petite.h>
#include <chemistry/qc/wfn/obwfn.h>


#ifndef DBL_EPSILON
#define DBL_EPSILON 1.0e-15
#endif

SavableState_REF_def(OneBodyWavefunction);

#define CLASSNAME OneBodyWavefunction
#define PARENTS public Wavefunction
#include <util/state/statei.h>
#include <util/class/classia.h>
void *
OneBodyWavefunction::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = Wavefunction::_castdown(cd);
  return do_castdowns(casts,cd);
}

OneBodyWavefunction::OneBodyWavefunction(const RefKeyVal&keyval):
  Wavefunction(keyval),
  eigenvectors_(this),
  eigenvalues_(this),
  density_(this),
  nirrep_(0),
  nvecperirrep_(0),
  occupations_(0),
  alpha_occupations_(0),
  beta_occupations_(0)
{
  double acc = keyval->doublevalue("eigenvector_accuracy");
  if (keyval->error() != KeyVal::OK)
    acc = 1.0e-7;
  
  eigenvectors_.set_desired_accuracy(acc);
  eigenvalues_.set_desired_accuracy(acc);

  if (eigenvectors_.desired_accuracy() < DBL_EPSILON) {
    eigenvectors_.set_desired_accuracy(DBL_EPSILON);
    eigenvalues_.set_desired_accuracy(DBL_EPSILON);
  }
}

OneBodyWavefunction::OneBodyWavefunction(StateIn&s):
  Wavefunction(s),
  eigenvectors_(this),
  eigenvalues_(this),
  density_(this),
  nirrep_(0),
  nvecperirrep_(0),
  occupations_(0),
  alpha_occupations_(0),
  beta_occupations_(0)
  maybe_SavableState(s)
{
  eigenvectors_.result_noupdate() =
    basis_matrixkit()->matrix(basis_dimension(), basis_dimension());
  eigenvectors_.restore_state(s);
  eigenvectors_.result_noupdate().restore(s);

  eigenvalues_.result_noupdate() =
    basis_matrixkit()->diagmatrix(basis_dimension());
  eigenvalues_.restore_state(s);
  eigenvalues_.result_noupdate().restore(s);

  density_.result_noupdate() =
    basis_matrixkit()->symmmatrix(basis_dimension());
  density_.restore_state(s);
  density_.result_noupdate().restore(s);
}

OneBodyWavefunction::~OneBodyWavefunction()
{
  if (nvecperirrep_) {
    delete[] nvecperirrep_;
    delete[] occupations_;
    delete[] alpha_occupations_;
    delete[] beta_occupations_;
  }
  nirrep_=0;
  nvecperirrep_=0;
  occupations_=0;
  alpha_occupations_=0;
  beta_occupations_=0;
}

void
OneBodyWavefunction::save_data_state(StateOut&s)
{
  Wavefunction::save_data_state(s);

  eigenvectors_.save_data_state(s);
  eigenvectors_.result_noupdate().save(s);

  eigenvalues_.save_data_state(s);
  eigenvalues_.result_noupdate().save(s);

  density_.save_data_state(s);
  density_.result_noupdate().save(s);
}

// at some point this will have to check for zero eigenvalues and not
// invert them
static void
form_m_half(RefSymmSCMatrix& M)
{
  // Diagonalize M to get m and U
  RefSCMatrix U(M.dim(), M.dim(), M.kit());
  RefDiagSCMatrix m(M.dim(), M.kit());
  M.diagonalize(m,U);

  // take square root of all elements of m
  RefSCElementOp op = new SCElementSquareRoot;
  m.element_op(op);

  // invert m
  op = new SCElementInvert;
  m.element_op(op);

  // back transform m^-1/2 to get M^-1/2 ( U*m^-1/2*U~ = M^-1/2)
  M.assign(0.0);
  M.accumulate_transform(U,m);
}

RefSCMatrix
OneBodyWavefunction::projected_eigenvectors(const RefOneBodyWavefunction& owfn,
                                            int alp)
{
  // first let's calculate the overlap between the new and old basis set
  integral()->set_basis(basis(), owfn->basis());
  RefSCMatrix s2(basis()->basisdim(), owfn->basis()->basisdim(),
                 basis()->matrixkit());

  RefSCElementOp op = new OneBodyIntOp(integral()->overlap());
  s2.assign(0.0);
  s2.element_op(op);
  op = 0;

  integral()->set_basis(basis());

  // now transform the AO s2 to the SO basis. this probably doesn't really
  // work.
  RefPetiteList pl = integral()->petite_list();
  RefPetiteList plo = integral()->petite_list(owfn->basis());
  RefSCMatrix s2b(pl->AO_basisdim(), plo->AO_basisdim(),
                  basis()->so_matrixkit());
  s2b->convert(s2);
  s2 = pl->aotoso().t() * s2b * plo->aotoso();

  //s2.print("overlap between new basis and old");
  s2b=0;

  BlockedSCMatrix *s2p = BlockedSCMatrix::require_castdown(s2.pointer(),
    "OneBodyWavefunction::projected_eigenvectors: s2p"
    );
  
  // now get the old eigenvectors and the new core hamiltonian guess
  RefSCMatrix ovec;
  if (owfn->spin_unrestricted()) {
    if (alp)
      ovec = owfn->alpha_eigenvectors();
    else
      ovec = owfn->beta_eigenvectors();
  }
  else
    ovec = owfn->eigenvectors();
    
  BlockedSCMatrix *ovecp = BlockedSCMatrix::require_castdown(ovec.pointer(),
    "OneBodyWavefunction::projected_eigenvectors: ovec"
    );
  //ovec.print("old wavefunction");
    
  RefSCMatrix vec = hcore_guess();
  BlockedSCMatrix *vecp = BlockedSCMatrix::require_castdown(vec.pointer(),
    "OneBodyWavefunction::projected_eigenvectors: vec"
    );
  //vec.print("hcore guess wavefunction");

  // we'll need the inverse of S eventually
  RefSymmSCMatrix s = overlap();
  BlockedSymmSCMatrix *sp = BlockedSymmSCMatrix::require_castdown(
    s.pointer(),
    "OneBodyWavefunction::projected_eigenvectors: s"
    );

  RefSymmSCMatrix sinv = s.i();
  BlockedSymmSCMatrix *sinvp = BlockedSymmSCMatrix::require_castdown(
    sinv.pointer(),
    "OneBodyWavefunction::projected_eigenvectors: sinv"
    );
  //sinv.print("Sinv");

  for (int irrep=0; irrep < vecp->nblocks(); irrep++) {
    // find out how many occupied orbitals there should be
    int nocc = 0;
    if (owfn->spin_unrestricted()) {
      if (alp)
        while (owfn->alpha_occupation(irrep,nocc) &&
               nocc < plo->nfunction(irrep)) nocc++;
      else
        while (owfn->beta_occupation(irrep,nocc) &&
               nocc < plo->nfunction(irrep)) nocc++;
    } else
      while (owfn->occupation(irrep,nocc) &&
             nocc < plo->nfunction(irrep)) nocc++;
  
    if (!nocc)
      continue;
    
    // get subblock of old vector
    RefSCMatrix ovecsb = ovecp->block(irrep).get_subblock(
      0, plo->nfunction(irrep)-1, 0, nocc-1);
    //ovecsb.print("old vec subblock");
  
    // form C' = S2 * Cold
    RefSCMatrix cprime = s2p->block(irrep) * ovecsb;
    //cprime.print("C' matrix");
  
    // we're done with ovecsb, free up some memory
    ovecsb=0;
  
    // now we need a matrix D = S^-1 * C' 
    RefSCMatrix D = sinvp->block(irrep) * cprime;
    //D.print("D matrix");
  
    // we also need X = C'~ * S^-1 * C'
    RefSymmSCMatrix X(cprime.coldim(),basis()->matrixkit());
    X.assign(0.0);
    X.accumulate_transform(cprime,sinvp->block(irrep),
                           SCMatrix::TransposeTransform);
    //X.print("X matrix");
  
    // we're done with cprime, free up some memory
    cprime=0;

    // now form X^-1/2
    form_m_half(X);
    //X.print("X^-1/2 matrix");

    // and form C'' = D * X^-1/2
    RefSCMatrix Cpp = D * X;
    //Cpp.print("new vector (occupied bits)");
  
    // we're done with X, free up some memory
    X=0;
    D=0;
  
    // stuff the projected bit into guess vector
    vecp->block(irrep).assign_subblock(Cpp,0,Cpp.nrow()-1,0,nocc-1);
    vecp->block(irrep)->schmidt_orthog(sp->block(irrep).pointer(),
                                       pl->nfunction(irrep));

    Cpp=0;
  }

  return vec;
}

// this is a hack for big basis sets where the core hamiltonian eigenvalues
// are total garbage.  Use the old wavefunction's occupied eigenvalues, and
// set all others to 99.

RefDiagSCMatrix
OneBodyWavefunction::projected_eigenvalues(const RefOneBodyWavefunction& owfn,
                                           int alp)
{
  // get the old eigenvalues and the new core hamiltonian evals
  RefDiagSCMatrix oval;
  if (owfn->spin_unrestricted()) {
    if (alp)
      oval = owfn->alpha_eigenvalues();
    else
      oval = owfn->beta_eigenvalues();
  } else
    oval = owfn->eigenvalues();

  BlockedDiagSCMatrix *ovalp =
    BlockedDiagSCMatrix::require_castdown(oval.pointer(),
      "OneBodyWavefunction::projected_eigenvalues: oval"
      );
    
  RefDiagSCMatrix val = core_hamiltonian().eigvals();
  BlockedDiagSCMatrix *valp =
    BlockedDiagSCMatrix::require_castdown(val.pointer(),
      "OneBodyWavefunction::projected_eigenvalues: val"
      );

  RefPetiteList pl = integral()->petite_list(basis());
  RefPetiteList plo = integral()->petite_list(owfn->basis());

  for (int irrep=0; irrep < valp->nblocks(); irrep++) {
    // find out how many occupied orbitals there should be
    int nocc = 0;
    if (owfn->spin_unrestricted()) {
      if (alp)
        while (owfn->alpha_occupation(irrep,nocc) &&
               nocc < plo->nfunction(irrep)) nocc++;
      else
        while (owfn->beta_occupation(irrep,nocc) &&
               nocc < plo->nfunction(irrep)) nocc++;
    } else
      while (owfn->occupation(irrep,nocc) &&
             nocc < plo->nfunction(irrep)) nocc++;
  
    int nf = pl->nfunction(irrep);
    int nfo = plo->nfunction(irrep);

    if (!nf)
      continue;
    
    double *vals = new double[pl->nfunction(irrep)];
    valp->block(irrep)->convert(vals);

    double *ovals;
    if (nfo) {
      ovals = new double[plo->nfunction(irrep)];
      ovalp->block(irrep)->convert(ovals);
    }
    
    int i;
    for (i=0; i < nocc; i++)
      vals[i] = ovals[i];
    for (i=nocc; i < pl->nfunction(irrep); i++)
      vals[i] = 99.0;

    valp->block(irrep)->assign(vals);

    delete[] vals;
    if (nfo)
      delete[] ovals;
  }

  return val;
}

RefSCMatrix
OneBodyWavefunction::hcore_guess()
{
  RefSymmSCMatrix hcore = core_hamiltonian();
  RefSCMatrix vec(basis_dimension(), basis_dimension(), basis_matrixkit());
  RefDiagSCMatrix val(basis_dimension(), basis_matrixkit());
  hcore.diagonalize(val,vec);

  vec = ao_to_orthog_ao() * vec;

  return vec;
}

// Function for returning an orbital value at a point
double
OneBodyWavefunction::orbital(const SCVector3& r, int iorb)
{
  return Wavefunction::orbital(r,iorb,eigenvectors());
}

// Function for returning an orbital value at a point
double
OneBodyWavefunction::orbital_density(const SCVector3& r,
                                            int iorb,
                                            double* orbval)
{
  return Wavefunction::orbital_density(r,iorb,eigenvectors(),orbval);
}

void
OneBodyWavefunction::print(ostream&o)
{
  Wavefunction::print(o);
}

void
OneBodyWavefunction::init_sym_info()
{
  RefPetiteList pl = integral()->petite_list();

  nirrep_ = pl->nirrep();
  nvecperirrep_ = new int[nirrep_];
  occupations_ = new double[basis()->nbasis()];
  alpha_occupations_ = new double[basis()->nbasis()];
  beta_occupations_ = new double[basis()->nbasis()];

  int ij=0;
  for (int i=0; i < nirrep_; i++) {
    nvecperirrep_[i] = pl->nfunction(i);

    // this is wrong...we should reorder the occupations so that this
    // vector represents the AO occupations...one day
    for (int j=0; j < nvecperirrep_[i]; j++, ij++) {
      occupations_[ij] = occupation(i,j);
      alpha_occupations_[ij] = alpha_occupation(i,j);
      beta_occupations_[ij] = beta_occupation(i,j);
    }
  }
}

double
OneBodyWavefunction::occupation(int vectornum)
{
  if (!nirrep_) init_sym_info();
  return occupations_[vectornum];
}

double
OneBodyWavefunction::alpha_occupation(int vectornum)
{
  if (!nirrep_) init_sym_info();
  return alpha_occupations_[vectornum];
}

double
OneBodyWavefunction::beta_occupation(int vectornum)
{
  if (!nirrep_) init_sym_info();
  return beta_occupations_[vectornum];
}

double
OneBodyWavefunction::alpha_occupation(int irrep, int vectornum)
{
  if (!spin_polarized())
    return 0.5*occupation(irrep, vectornum);
  
  cerr << class_name() << "::alpha_occupation not implemented" << endl;
  abort();
  return 0;
}

double
OneBodyWavefunction::beta_occupation(int irrep, int vectornum)
{
  if (!spin_polarized())
    return 0.5*occupation(irrep, vectornum);
  
  cerr << class_name() << "::beta_occupation not implemented" << endl;
  abort();
  return 0;
}

RefSCMatrix
OneBodyWavefunction::alpha_eigenvectors()
{
  if (!spin_unrestricted())
    return eigenvectors().copy();

  cerr << class_name() << "::alpha_eigenvectors not implemented" << endl;
  abort();
  return 0;
}

RefSCMatrix
OneBodyWavefunction::beta_eigenvectors()
{
  if (!spin_unrestricted())
    return eigenvectors().copy();

  cerr << class_name() << "::beta_eigenvectors not implemented" << endl;
  abort();
  return 0;
}

RefDiagSCMatrix
OneBodyWavefunction::alpha_eigenvalues()
{
  if (!spin_unrestricted())
    return eigenvalues().copy();

  cerr << class_name() << "::alpha_eigenvalues not implemented" << endl;
  abort();
  return 0;
}

RefDiagSCMatrix
OneBodyWavefunction::beta_eigenvalues()
{
  if (!spin_unrestricted())
    return eigenvalues().copy();

  cerr << class_name() << "::beta_eigenvalues not implemented" << endl;
  abort();
  return 0;
}

int
OneBodyWavefunction::nelectron()
{
  int nbasis = basis()->nbasis();
  double tocc = 0.0;
  if (!spin_polarized()) {
    for (int i=0; i<nbasis; i++) {
      tocc += occupation(i);
    }
  }
  else {
    for (int i=0; i<nbasis; i++) {
      tocc += alpha_occupation(i) + beta_occupation(i);
    }
  }
  return int(tocc+0.5);
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:

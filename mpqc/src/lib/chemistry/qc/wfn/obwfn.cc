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

#include <util/state/stateio.h>
#include <util/misc/formio.h>

#include <math/symmetry/corrtab.h>
#include <math/scmat/local.h>
#include <math/scmat/blocked.h>

#include <chemistry/qc/basis/integral.h>
#include <chemistry/qc/basis/obint.h>
#include <chemistry/qc/basis/petite.h>
#include <chemistry/qc/wfn/obwfn.h>

#define DEBUG 0

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
  density_(this),
  eigenvectors_(this),
  eigenvalues_(this),
  nirrep_(0),
  nvecperirrep_(0),
  occupations_(0),
  alpha_occupations_(0),
  beta_occupations_(0)
{
  double acc = keyval->doublevalue("eigenvector_accuracy");
  if (keyval->error() != KeyVal::OK)
    acc = value_.desired_accuracy();
  
  eigenvectors_.set_desired_accuracy(acc);
  eigenvalues_.set_desired_accuracy(acc);

  if (eigenvectors_.desired_accuracy() < DBL_EPSILON) {
    eigenvectors_.set_desired_accuracy(DBL_EPSILON);
    eigenvalues_.set_desired_accuracy(DBL_EPSILON);
  }
}

OneBodyWavefunction::OneBodyWavefunction(StateIn&s):
  maybe_SavableState(s)
  Wavefunction(s),
  density_(this),
  eigenvectors_(this),
  eigenvalues_(this),
  nirrep_(0),
  nvecperirrep_(0),
  occupations_(0),
  alpha_occupations_(0),
  beta_occupations_(0)
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

#if 1
RefSCMatrix
OneBodyWavefunction::projected_eigenvectors(const RefOneBodyWavefunction& owfn,
                                            int alp)
{
  //............................................................
  // first obtain the guess density matrix

  // The old density in the old SO basis
  RefSymmSCMatrix oldP_so;
  if (owfn->spin_unrestricted()) {
    if (alp) oldP_so = owfn->alpha_density();
    else oldP_so = owfn->beta_density();
  }
  else oldP_so = owfn->density();

  cout << node0 << endl << indent
       << "Projecting the guess density.\n"
       << endl;
  cout << incindent;

  // The old overlap
  RefSymmSCMatrix oldS = owfn->overlap();
  cout << node0 << indent
       << "The number of electrons in the guess density = "
       << (oldP_so*oldS).trace() << endl;

  // Transform the old SO overlap into the orthogonal SO basis, oSO
  RefSymmSCMatrix old_so_to_oso = owfn->so_to_orthog_so();
  RefSymmSCMatrix oldP_oso = oldP_so->clone();
  oldP_oso->assign(0.0);
  oldP_oso->accumulate_transform(old_so_to_oso, oldP_so);

  //............................................................
  // transform the guess density into the current basis

  // the transformation matrix is the new basis/old basis overlap
  integral()->set_basis(owfn->basis(), basis());
  RefSCMatrix old_to_new_ao(owfn->basis()->basisdim(), basis()->basisdim(),
                            basis()->matrixkit());
  RefSCElementOp op = new OneBodyIntOp(integral()->overlap());
  old_to_new_ao.assign(0.0);
  old_to_new_ao.element_op(op);
  op = 0;
  integral()->set_basis(basis());

  // now must transform the transform into the SO basis
  RefPetiteList pl = integral()->petite_list();
  RefPetiteList oldpl = owfn->integral()->petite_list();
  RefSCMatrix blocked_old_to_new_ao(oldpl->AO_basisdim(), pl->AO_basisdim(),
                                    basis()->so_matrixkit());
  blocked_old_to_new_ao->convert(old_to_new_ao);
  RefSCMatrix old_to_new_so
    = oldpl->sotoao() * blocked_old_to_new_ao * pl->aotoso();

  // now must transform the transform into the orthogonal SO basis
  RefSymmSCMatrix so_to_oso = so_to_orthog_so();
  RefSCMatrix old_to_new_oso = old_so_to_oso.gi()
                             * old_to_new_so
                             * so_to_oso;
  old_so_to_oso = 0;
  old_to_new_so = 0;
  
  // The old density transformed to the new orthogonal SO basis
  RefSymmSCMatrix newP_oso(basis_dimension(), basis_matrixkit());
  newP_oso->assign(0.0);
  newP_oso->accumulate_transform(old_to_new_oso.t(), oldP_oso);
  old_to_new_oso = 0;
  oldP_oso = 0;
  //newP_oso.print("projected orthoSO density");

  cout << node0 << indent
       << "The number of electrons in the projected density = "
       << newP_oso.trace() << endl;

  //............................................................

  // reverse the sign of the density so the eigenvectors will
  // be ordered in the right way
  newP_oso.scale(-1.0);

  // use the guess density in the new basis to find the orbitals
  // (the density should be diagonal in the MO basis--this proceedure
  // will not give canonical orbitals, but they should at least give
  // a decent density)
  RefDiagSCMatrix newP_oso_vals(newP_oso.dim(), basis_matrixkit());
  RefSCMatrix newP_oso_vecs(newP_oso.dim(), newP_oso.dim(), basis_matrixkit());
  newP_oso.diagonalize(newP_oso_vals, newP_oso_vecs);
  //newP_oso_vals.print("eigenvalues of projected density");

  // Reordering of the vectors isn't needed because of the way
  // the density was scaled above.
  RefSCMatrix newvec_oso = newP_oso_vecs;

  // transform the orbitals from the orthogonal to nonorthogonal SO basis
  RefSCMatrix newvec_so = so_to_oso * newvec_oso;

  if (debug_ >= 2) {
      newvec_oso.print("projected ortho SO vector");
      so_to_oso.print("SO to ortho SO transformation");
      newvec_so.print("projected SO vector");
    }

  cout << decindent;
  return newvec_so;
}
#else
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
#if DEBUG
  ovec.print("old wavefunction");
#endif

  RefSCMatrix vec = hcore_guess();
  BlockedSCMatrix *vecp = BlockedSCMatrix::require_castdown(vec.pointer(),
    "OneBodyWavefunction::projected_eigenvectors: vec"
    );
#if DEBUG
  vec.print("hcore guess wavefunction");
#endif

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
#if DEBUG
  sinv.print("Sinv");
#endif

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
#if DEBUG
    ovecsb.print("old vec subblock");
#endif

    // form C' = S2 * Cold
    RefSCMatrix cprime = s2p->block(irrep) * ovecsb;
#if DEBUG
    cprime.print("C' matrix");
#endif

    // we're done with ovecsb, free up some memory
    ovecsb=0;
  
    // now we need a matrix D = S^-1 * C' 
    RefSCMatrix D = sinvp->block(irrep) * cprime;
#if DEBUG
    D.print("D matrix");
#endif

    // we also need X = C'~ * S^-1 * C'
    RefSymmSCMatrix X(cprime.coldim(),basis()->matrixkit());
    X.assign(0.0);
    X.accumulate_transform(cprime,sinvp->block(irrep),
                           SCMatrix::TransposeTransform);
#if DEBUG
    X.print("X matrix");
#endif

    // we're done with cprime, free up some memory
    cprime=0;

    // now form X^-1/2
    form_m_half(X);
#if DEBUG
    X.print("X^-1/2 matrix");
#endif

    // and form C'' = D * X^-1/2
    RefSCMatrix Cpp = D * X;
#if DEBUG
    Cpp.print("new vector (occupied bits)");
#endif
  
    // we're done with X, free up some memory
    X=0;
    D=0;
  
    // stuff the projected bit into guess vector
    vecp->block(irrep).assign_subblock(Cpp,0,Cpp.nrow()-1,0,nocc-1);
    vecp->block(irrep)->schmidt_orthog(sp->block(irrep).pointer(),
                                       pl->nfunction(irrep));

    Cpp=0;
  }

#if DEBUG
  vec.print("projected vectors");
#endif

  return vec;
}
#endif

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

    int i;
    if (nfo) {
      double *ovals = new double[plo->nfunction(irrep)];
      ovalp->block(irrep)->convert(ovals);
      for (i=0; i < nocc; i++) vals[i] = ovals[i];
      delete[] ovals;
    }
    
    for (i=nocc; i < pl->nfunction(irrep); i++)
      vals[i] = 99.0;

    valp->block(irrep)->assign(vals);

    delete[] vals;
  }

#if DEBUG
  val.print("projected values");
#endif

  return val;
}

RefSCMatrix
OneBodyWavefunction::hcore_guess()
{
  RefSymmSCMatrix hcore = core_hamiltonian();
  RefSCMatrix vec(basis_dimension(), basis_dimension(), basis_matrixkit());
  RefDiagSCMatrix val(basis_dimension(), basis_matrixkit());
  hcore.diagonalize(val,vec);

  vec = so_to_orthog_so() * vec;

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
OneBodyWavefunction::print(ostream&o) const
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

void
OneBodyWavefunction::symmetry_changed()
{
  Wavefunction::symmetry_changed();

  // junk the old occupation information
  delete[] nvecperirrep_;
  delete[] occupations_;
  delete[] alpha_occupations_;
  delete[] beta_occupations_;
  nirrep_ = 0;
  nvecperirrep_=0;
  occupations_=0;
  alpha_occupations_=0;
  beta_occupations_=0;
  
  // for now, delete old eigenvectors...later we'll transform to new
  // pointgroup
  eigenvectors_.result_noupdate() = 0;
}

int
OneBodyWavefunction::form_occupations(int *&newocc, const int *oldocc)
{
  delete[] newocc;
  newocc = 0;

  CorrelationTable corrtab;
  if (corrtab.initialize_table(initial_pg_, molecule()->point_group()))
    return 0;

  newocc = new int[corrtab.subn()];
  memset(newocc,0,sizeof(int)*corrtab.subn());

  for (int i=0; i<corrtab.n(); i++) {
      for (int j=0; j<corrtab.ngamma(i); j++) {
          int gam = corrtab.gamma(i,j);
          newocc[gam] += (corrtab.subdegen(gam)*oldocc[i])/corrtab.degen(i);
        }
    }

  return 1;
}
  
/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:

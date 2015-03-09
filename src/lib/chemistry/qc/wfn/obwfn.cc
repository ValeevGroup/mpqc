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

#include <cassert>
#include <util/state/stateio.h>
#include <util/misc/formio.h>

#include <math/symmetry/corrtab.h>
#include <math/scmat/local.h>
#include <math/scmat/blocked.h>

#include <chemistry/qc/basis/integral.h>
#include <chemistry/qc/basis/obint.h>
#include <chemistry/qc/basis/petite.h>
#include <chemistry/qc/wfn/obwfn.h>
#include <chemistry/qc/wfn/orbitalspace.h>
#include <util/misc/scexception.h>


using namespace std;
using namespace sc;

#define DEBUG 0

#ifndef DBL_EPSILON
#define DBL_EPSILON 1.0e-15
#endif

static ClassDesc OneBodyWavefunction_cd(
  typeid(OneBodyWavefunction),"OneBodyWavefunction",1,"public Wavefunction",
  0, 0, 0);

OneBodyWavefunction::OneBodyWavefunction(const Ref<KeyVal>&keyval):
  Wavefunction(keyval),
  density_(this),
  oso_eigenvectors_(this),
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

  oso_eigenvectors_.set_desired_accuracy(acc);
  eigenvalues_.set_desired_accuracy(acc);

  if (oso_eigenvectors_.desired_accuracy() < DBL_EPSILON) {
    oso_eigenvectors_.set_desired_accuracy(DBL_EPSILON);
    eigenvalues_.set_desired_accuracy(DBL_EPSILON);
  }

}

OneBodyWavefunction::OneBodyWavefunction(StateIn&s):
  SavableState(s),
  Wavefunction(s),
  density_(this),
  oso_eigenvectors_(this),
  eigenvalues_(this),
  nirrep_(0),
  nvecperirrep_(0),
  occupations_(0),
  alpha_occupations_(0),
  beta_occupations_(0)
{
  oso_eigenvectors_.result_noupdate() =
    basis_matrixkit()->matrix(oso_dimension(), oso_dimension());
  oso_eigenvectors_.restore_state(s);
  oso_eigenvectors_.result_noupdate().restore(s);

  eigenvalues_.result_noupdate() =
    basis_matrixkit()->diagmatrix(oso_dimension());
  eigenvalues_.restore_state(s);
  eigenvalues_.result_noupdate().restore(s);

  density_.result_noupdate() =
    basis_matrixkit()->symmmatrix(so_dimension());
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

  oso_eigenvectors_.save_data_state(s);
  oso_eigenvectors_.result_noupdate().save(s);

  eigenvalues_.save_data_state(s);
  eigenvalues_.result_noupdate().save(s);

  density_.save_data_state(s);
  density_.result_noupdate().save(s);
}

RefSCMatrix
OneBodyWavefunction::projected_eigenvectors(const Ref<OneBodyWavefunction>& owfn,
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

  ExEnv::out0() << endl << indent
       << "Projecting the guess density.\n"
       << endl;
  ExEnv::out0() << incindent;

  // The old overlap
  RefSymmSCMatrix oldS = owfn->overlap();
  ExEnv::out0() << indent
       << "The number of electrons in the guess density = "
       << (oldP_so*oldS).trace() << endl;

  // Transform the old SO overlap into the orthogonal SO basis, oSO
  RefSCMatrix old_so_to_oso = owfn->so_to_orthog_so();
  RefSymmSCMatrix oldP_oso(owfn->oso_dimension(), owfn->basis_matrixkit());
  oldP_oso->assign(0.0);
  oldP_oso->accumulate_transform(old_so_to_oso, oldP_so);

  //............................................................
  // transform the guess density into the current basis

  // the transformation matrix is the new basis/old basis overlap
  integral()->set_basis(owfn->basis(), basis());
  RefSCMatrix old_to_new_ao(owfn->basis()->basisdim(), basis()->basisdim(),
                            basis()->matrixkit());
  Ref<SCElementOp> op = new OneBodyIntOp(integral()->overlap());
  old_to_new_ao.assign(0.0);
  old_to_new_ao.element_op(op);
  op = 0;

  // now must transform the transform into the SO basis
  owfn->integral()->set_basis(owfn->basis());
  Ref<PetiteList> oldpl = owfn->integral()->petite_list();
  integral()->set_basis(basis());
  Ref<PetiteList> pl = integral()->petite_list();
  RefSCMatrix blocked_old_to_new_ao(oldpl->AO_basisdim(), pl->AO_basisdim(),
                                    basis()->so_matrixkit());
  blocked_old_to_new_ao->convert(old_to_new_ao);
  RefSCMatrix old_to_new_so
    = oldpl->sotoao() * blocked_old_to_new_ao * pl->aotoso();

  // now must transform the transform into the orthogonal SO basis
  RefSCMatrix so_to_oso = so_to_orthog_so();
  RefSCMatrix old_to_new_oso = owfn->so_to_orthog_so_inverse().t()
                             * old_to_new_so
                             * so_to_oso.t();
  old_so_to_oso = 0;
  old_to_new_so = 0;

  // The old density transformed to the new orthogonal SO basis
  RefSymmSCMatrix newP_oso(oso_dimension(), basis_matrixkit());
  newP_oso->assign(0.0);
  newP_oso->accumulate_transform(old_to_new_oso.t(), oldP_oso);
  old_to_new_oso = 0;
  oldP_oso = 0;
  //newP_oso.print("projected orthoSO density");

  ExEnv::out0() << indent
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

  if (debug_ >= 2) {
      newvec_oso.print("projected ortho SO vector");
      so_to_oso.print("SO to ortho SO transformation");
    }

  ExEnv::out0() << decindent;
  return newvec_oso;
}

// this is a hack for big basis sets where the core hamiltonian eigenvalues
// are total garbage.  Use the old wavefunction's occupied eigenvalues, and
// set all others to 99.

RefDiagSCMatrix
OneBodyWavefunction::projected_eigenvalues(const Ref<OneBodyWavefunction>& owfn,
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
    require_dynamic_cast<BlockedDiagSCMatrix*>(oval.pointer(),
      "OneBodyWavefunction::projected_eigenvalues: oval"
      );

  // get the core hamiltonian eigenvalues
  RefDiagSCMatrix val;
  hcore_guess(val);
  BlockedDiagSCMatrix *valp =
    require_dynamic_cast<BlockedDiagSCMatrix*>(val.pointer(),
      "OneBodyWavefunction::projected_eigenvalues: val"
      );

  RefSCDimension oso = oso_dimension();
  RefSCDimension ooso = owfn->oso_dimension();

  for (int irrep=0; irrep < valp->nblocks(); irrep++) {
    // find out how many occupied orbitals there should be

    int nf = oso->blocks()->size(irrep);
    int nfo = ooso->blocks()->size(irrep);

    int nocc = 0;
    if (owfn->spin_unrestricted()) {
      if (alp)
        while (owfn->alpha_occupation(irrep,nocc) &&
               nocc < nfo) nocc++;
      else
        while (owfn->beta_occupation(irrep,nocc) &&
               nocc < nfo) nocc++;
    } else
      while (owfn->occupation(irrep,nocc) &&
             nocc < nfo) nocc++;

    if (!nf)
      continue;

    double *vals = new double[nf];
    valp->block(irrep)->convert(vals);

    int i;
    if (nfo) {
      double *ovals = new double[nfo];
      ovalp->block(irrep)->convert(ovals);
      for (i=0; i < nocc; i++) vals[i] = ovals[i];
      delete[] ovals;
    }

    for (i=nocc; i < nf; i++)
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
OneBodyWavefunction::so_to_mo()
{
  // works for transforming H, S, etc (covariant)
  return orthog_so_to_mo() * so_to_orthog_so();
  // works for transforming the Density (contravariant)
  //return orthog_so_to_mo() * so_to_orthog_so_inverse().t();
}

RefSCMatrix
OneBodyWavefunction::orthog_so_to_mo()
{
  return oso_eigenvectors().t();
}

RefSCMatrix
OneBodyWavefunction::mo_to_so()
{
  // works for transforming H, S, etc (covariant)
  return so_to_orthog_so_inverse() * mo_to_orthog_so();
  // works for transforming the Density (contravariant)
  //return so_to_orthog_so().t() * mo_to_orthog_so();
}

RefSCMatrix
OneBodyWavefunction::mo_to_orthog_so()
{
  return oso_eigenvectors();
}

RefSCMatrix
OneBodyWavefunction::eigenvectors()
{
  return so_to_orthog_so().t() * oso_eigenvectors();
}

RefSCMatrix
OneBodyWavefunction::hcore_guess()
{
  RefDiagSCMatrix val;
  return hcore_guess(val);
}

RefSCMatrix
OneBodyWavefunction::hcore_guess(RefDiagSCMatrix &val)
{
  RefSCMatrix vec(oso_dimension(), oso_dimension(), basis_matrixkit());
  val = basis_matrixkit()->diagmatrix(oso_dimension());

  // I'm about to do something strange, but it will only work
  // if the SO and orthogonal SO dimensions are equivalent.  This
  // is not the case for canonical orthogonalization when there
  // are linear dependencies.
  if (so_dimension()->equiv(oso_dimension())) {
    // Yes, this is diagonalizing Hcore in a nonorthogonal basis
    // and does not really make any sense--except it seems to
    // always give a better initial guess.  I don't understand
    // why it works better.
    core_hamiltonian().diagonalize(val,vec);
    if (debug_ > 1) {
	  val.print("hcore eigenvalues in SO basis");
      vec.print("hcore eigenvectors in SO basis");
    }
  }
  else {
    RefSymmSCMatrix hcore_oso(oso_dimension(), basis_matrixkit());
    hcore_oso->assign(0.0);
    hcore_oso->accumulate_transform(so_to_orthog_so(), core_hamiltonian());

    if (debug_ > 1) {
      hcore_oso.print("hcore in ortho SO basis");
    }

    hcore_oso.diagonalize(val,vec);

    if (debug_ > 1) {
      val.print("hcore eigenvalues in ortho SO basis");
      vec.print("hcore eigenvectors in ortho SO basis");
    }
  }

  return vec;
}

// Function for returning an orbital value at a point
double
OneBodyWavefunction::orbital(const SCVector3& r, int iorb)
{
  Ref<PetiteList> pl = integral()->petite_list();
  RefSCMatrix ao_orbital_coeff = pl->evecs_to_AO_basis(this->so_to_mo().t());
  return Wavefunction::orbital(r,iorb, ao_orbital_coeff);
}

void OneBodyWavefunction::orbitals(const std::vector<SCVector3> & Points,
                                   std::vector<double>& Vals,
                                   unsigned int first, unsigned int last,
                                   bool energy_ordered)
{
  Ref<PetiteList> pl = integral()->petite_list();
  RefSCMatrix aocoefs_full = pl->evecs_to_AO_basis(this->so_to_mo().t());
  RefSCMatrix aocoefs;
  if (energy_ordered) {
    Ref<OrbitalSpace> mospace = new OrbitalSpace("p", "energy-ordered MOs to evaluate",
                                                 aocoefs_full, this->basis(), this->integral(),
                                                 this->eigenvalues(),
                                                 first, aocoefs_full.ncol() - last - 1);
    aocoefs = mospace->coefs_nb();
  }
  else {
    aocoefs = this->matrixkit()->matrix(aocoefs_full.rowdim(),
                                        new SCDimension(last - first + 1));
    for(int r=0; r<aocoefs.nrow(); ++r)
      for(int c=0; c<aocoefs.ncol(); ++c)
        aocoefs(r,c) = aocoefs_full(r,c+first);  // select mos \in [first,last]
  }

  // Wavefunction::orbitals wants transposed MOs (nmo by nao)
  aocoefs = aocoefs.t();
  const int numpoints = Points.size();
  const int nmo = aocoefs.nrow();
  Vals.resize(numpoints * nmo);
  int count = 0;
  RefSCVector values = aocoefs->kit()->vector(aocoefs.rowdim());
  for(int i1 = 0; i1 < numpoints; ++i1) {
    Wavefunction::orbitals(Points[i1], aocoefs, values);
    values.convert(&(Vals[count]));
    count += nmo;
  }
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
  RefSCDimension d = oso_dimension();
  nirrep_ = d->blocks()->nblock();
  nvecperirrep_ = new int[nirrep_];
  occupations_ = new double[d->n()];
  alpha_occupations_ = new double[d->n()];
  beta_occupations_ = new double[d->n()];

  int ij=0;
  for (int i=0; i < nirrep_; i++) {
    nvecperirrep_[i] = d->blocks()->size(i);

    for (int j=0; j < nvecperirrep_[i]; j++, ij++) {
      if (!spin_unrestricted()) occupations_[ij] = occupation(i,j);
      else occupations_[ij] = 0.0;
      alpha_occupations_[ij] = alpha_occupation(i,j);
      beta_occupations_[ij] = beta_occupation(i,j);
    }
  }
}

double
OneBodyWavefunction::occupation(int vectornum)
{
  if (spin_unrestricted()) {
    ExEnv::errn() << "OneBodyWavefunction::occupation: called for USCF case"
                 << endl;
    abort();
  }
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

  ExEnv::errn() << class_name() << "::alpha_occupation not implemented" << endl;
  abort();
  return 0;
}

double
OneBodyWavefunction::beta_occupation(int irrep, int vectornum)
{
  if (!spin_polarized())
    return 0.5*occupation(irrep, vectornum);

  ExEnv::errn() << class_name() << "::beta_occupation not implemented" << endl;
  abort();
  return 0;
}

RefSCMatrix
OneBodyWavefunction::oso_alpha_eigenvectors()
{
  if (!spin_unrestricted())
    return oso_eigenvectors().copy();

  ExEnv::errn() << class_name() << "::oso_alpha_eigenvectors not implemented" << endl;
  abort();
  return 0;
}

RefSCMatrix
OneBodyWavefunction::oso_beta_eigenvectors()
{
  if (!spin_unrestricted())
    return oso_eigenvectors().copy();

  ExEnv::errn() << class_name() << "::oso_beta_eigenvectors not implemented" << endl;
  abort();
  return 0;
}

RefSCMatrix
OneBodyWavefunction::alpha_eigenvectors()
{
  if (!spin_unrestricted())
    return eigenvectors().copy();

  ExEnv::errn() << class_name() << "::alpha_eigenvectors not implemented" << endl;
  abort();
  return 0;
}

RefSCMatrix
OneBodyWavefunction::beta_eigenvectors()
{
  if (!spin_unrestricted())
    return eigenvectors().copy();

  ExEnv::errn() << class_name() << "::beta_eigenvectors not implemented" << endl;
  abort();
  return 0;
}

RefDiagSCMatrix
OneBodyWavefunction::alpha_eigenvalues()
{
  if (!spin_unrestricted())
    return eigenvalues().copy();

  ExEnv::errn() << class_name() << "::alpha_eigenvalues not implemented" << endl;
  abort();
  return 0;
}

RefDiagSCMatrix
OneBodyWavefunction::beta_eigenvalues()
{
  if (!spin_unrestricted())
    return eigenvalues().copy();

  ExEnv::errn() << class_name() << "::beta_eigenvalues not implemented" << endl;
  abort();
  return 0;
}

int
OneBodyWavefunction::nelectron()
{
  int noso = oso_dimension()->n();
  double tocc = 0.0;
  if (not spin_unrestricted()) {
    for (int i=0; i<noso; i++) {
      tocc += occupation(i);
    }
  }
  else {
    for (int i=0; i<noso; i++) {
      tocc += alpha_occupation(i) + beta_occupation(i);
    }
  }
  MPQC_ASSERT((tocc - static_cast<int>(tocc)) < 1e-8);
  return static_cast<int>(tocc);
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
  oso_eigenvectors_.result_noupdate() = 0;
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

void
OneBodyWavefunction::set_desired_value_accuracy(double eps)
{
  Function::set_desired_value_accuracy(eps);
  oso_eigenvectors_.set_desired_accuracy(eps);
  eigenvalues_.set_desired_accuracy(eps);
}



/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:

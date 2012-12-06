//
// taylor.cc
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

#include <util/misc/scexception.h>
#include <util/misc/formio.h>
#include <util/state/stateio.h>
#include <math/scmat/local.h>
#include <chemistry/molecule/taylor.h>

using namespace std;
using namespace sc;

static ClassDesc TaylorMolecularEnergy_cd(
  typeid(TaylorMolecularEnergy),"TaylorMolecularEnergy",1,"public MolecularEnergy",
  0, create<TaylorMolecularEnergy>, create<TaylorMolecularEnergy>);

// Note:  this gets the values of the coordinates from the current molecule
// rather than the coordinates.
TaylorMolecularEnergy::TaylorMolecularEnergy(const Ref<KeyVal>&keyval):
  MolecularEnergy(keyval)
{
  coordinates_ << keyval->describedclassvalue("coordinates");
  // if coordinates is nonnull use cartesian coordinates
  if (coordinates_.nonnull()) {
      dim_ = new SCDimension(coordinates_->n());
    }
  else {
      dim_ = moldim();
    }
  if (coordinates_.nonnull()) {
      expansion_point_ = matrixkit()->vector(dim_);
      coordinates_->update_values(molecule());
      coordinates_->values_to_vector(expansion_point_);
    }
  else {
      expansion_point_ = get_cartesian_x();
    }

  e0_ = keyval->doublevalue("e0");

  // count the number of force constants
  int i;
  int n_fc1 = keyval->count("force_constants_1");
  int n_fc = keyval->count("force_constants");

  int use_guess_hessian = 0;
  if (coordinates_.null() && n_fc == 0) {
      use_guess_hessian = 1;
      n_fc = (moldim().n()*(moldim().n()+1))/2;
      maxorder_ = 2;
    }

  force_constant_index_.resize(n_fc1+n_fc);
  force_constant_value_.resize(n_fc1+n_fc);
  maxorder_ = 0;
  if (n_fc1 > 0) maxorder_ = 1;

  // first read in the short hand notation for first derivatives
  for (i=0; i<n_fc1; i++) {
      force_constant_value_[i] = keyval->doublevalue("force_constants_1", i);
      force_constant_index_[i].resize(1);
      force_constant_index_[i][0] = i;
    }

  if (use_guess_hessian) {
      RefSymmSCMatrix hess(moldim(), matrixkit());
      guess_hessian(hess);
      int ifc,j;
      for (ifc=i=0; i<moldim().n(); i++) {
          for (j=0; j<=i; j++, ifc++) {
              force_constant_index_[n_fc1+ifc].resize(2);
              force_constant_index_[n_fc1+ifc][0] = i;
              force_constant_index_[n_fc1+ifc][1] = j;
              force_constant_value_[n_fc1+ifc] = hess->get_element(i,j);
            }
        }
    }
  else {
      // read in the general force constants
      for (i=0; i<n_fc; i++) {
          int order = keyval->count("force_constants", i) - 1;
          force_constant_value_[n_fc1+i]
              = keyval->doublevalue("force_constants", i,order);
          force_constant_index_[n_fc1+i].resize(order);
          if (maxorder_ < order) maxorder_ = order;
          for (int j=0; j<order; j++) {
              force_constant_index_[n_fc1+i][j]
                  = keyval->intvalue("force_constants",i,j) - 1;
            }
        }
    }
}

TaylorMolecularEnergy::~TaylorMolecularEnergy()
{
}

TaylorMolecularEnergy::TaylorMolecularEnergy(StateIn&s):
  SavableState(s),
  MolecularEnergy(s)
{
  throw ProgrammingError("cannot save state for this class",
                         __FILE__, __LINE__, class_desc());
}

void
TaylorMolecularEnergy::save_data_state(StateOut&s)
{
  MolecularEnergy::save_data_state(s);
  throw ProgrammingError("cannot save state for this class",
                         __FILE__, __LINE__, class_desc());
}

void
TaylorMolecularEnergy::print(ostream&o) const
{
  MolecularEnergy::print(o);
  if (coordinates_.nonnull()) coordinates_->print_details(molecule(), o);
  int nfc = force_constant_index_.size();
  o << indent << "Force Constants:" << endl;
  o << incindent;
  for (int i=0; i<nfc; i++) {
      int order = force_constant_index_[i].size();
      for (int j=0; j<order; j++) {
          o << indent << scprintf("%5d",force_constant_index_[i][j]+1);
        }
      o << indent
        << scprintf(" %*.*f",14,10,force_constant_value_[i])
        << endl;
    }
  o << decindent;
}

// this is used by the factor function
static int
factorial(int i)
{
  if (i<2) return 1;
  return i*factorial(i-1);
}

// Compute the factors such as 1/4!, etc. assuming that only unique
// force constants we given.
static double
factor(std::vector<int>&indices)
{
  std::map<int,int> n_occur;
  int i;
  for (i=0; i<indices.size(); i++) {
      n_occur[indices[i]] = 0;
    }
  for (i=0; i<indices.size(); i++) {
      n_occur[indices[i]]++;
    }
  int n_indices = indices.size();
  int int_factor = 1;
  std::map<int,int>::iterator I;
  for (I=n_occur.begin(); I!=n_occur.end(); I++) {
      int n = I->second;
      int_factor *= factorial(n_indices)
                   /(factorial(n)*factorial(n_indices-n));
      n_indices -= n;
    }
  double term = ((double)int_factor) / factorial(indices.size());
  return term;
}

void
TaylorMolecularEnergy::compute()
{
  RefSCVector geometry;

  if (coordinates_.nonnull()) {
      coordinates_->update_values(molecule());
      geometry = expansion_point_.clone();
      coordinates_->values_to_vector(geometry);
    }
  else {
      geometry = get_cartesian_x();
    }

  RefSCVector displacement = geometry - expansion_point_;

  if (value_needed()) {
      double e = e0_;
      for (int i=0; i<force_constant_index_.size(); i++) {
          double term =  force_constant_value_[i]
                       * factor(force_constant_index_[i]);
          for (int j=0; j<force_constant_index_[i].size(); j++) {
              term *= displacement(force_constant_index_[i][j]);
            }
          e += term;
        }
      set_energy(e);
      set_actual_value_accuracy(desired_value_accuracy());
    }
  if (gradient_needed()) {
      RefSCVector gradient = expansion_point_.clone();
      gradient.assign(0.0);

      for (int i=0; i<force_constant_index_.size(); i++) {
          double f =  force_constant_value_[i]
                    * factor(force_constant_index_[i]);
          for (int k=0; k<force_constant_index_[i].size(); k++) {
              double t = 1.0;
              for (int l=0; l<force_constant_index_[i].size(); l++) {
                  if (l == k) continue;
                  t *= displacement(force_constant_index_[i][l]);
                }
              gradient.accumulate_element(force_constant_index_[i][k],f*t);
            }
        }

      // this will only work for cartesian coordinates
      set_gradient(gradient);
      set_actual_gradient_accuracy(desired_gradient_accuracy());
    }
}

int
TaylorMolecularEnergy::value_implemented() const
{
  return 1;
}

bool
TaylorMolecularEnergy::analytic_gradient_implemented() const
{
  return coordinates_.null() && maxorder_ >= 1;
}

bool
TaylorMolecularEnergy::analytic_hessian_implemented() const
{
  return false;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:

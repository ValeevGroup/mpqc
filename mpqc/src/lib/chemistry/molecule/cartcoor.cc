//
// cartcoor.cc
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

#include <math.h>

#include <math/scmat/matrix.h>
#include <chemistry/molecule/molecule.h>
#include <chemistry/molecule/coor.h>
#include <chemistry/molecule/simple.h>

///////////////////////////////////////////////////////////////////////////
// members of CartMolecularCoor

#define CLASSNAME CartMolecularCoor
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#define PARENTS public MolecularCoor
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
CartMolecularCoor::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = MolecularCoor::_castdown(cd);
  return do_castdowns(casts,cd);
}

CartMolecularCoor::CartMolecularCoor(RefMolecule&mol):
  MolecularCoor(mol)
{
  init();
}

CartMolecularCoor::CartMolecularCoor(const RefKeyVal& keyval):
  MolecularCoor(keyval)
{
  init();
}

CartMolecularCoor::CartMolecularCoor(StateIn& s):
  MolecularCoor(s)
{
  dim_.restore_state(s);
}

void
CartMolecularCoor::init()
{
  Molecule& m = *molecule_.pointer();

  // compute needed dimensions
  dim_ = dnatom3_;
}

CartMolecularCoor::~CartMolecularCoor()
{
}

void
CartMolecularCoor::save_data_state(StateOut&s)
{
  MolecularCoor::save_data_state(s);

  dim_.save_state(s);
}

RefSCDimension
CartMolecularCoor::dim()
{
  return dim_;
}


// presumably this will actually be passed the new cartesian coords in
// new_internal, so do almost nothing
int
CartMolecularCoor::to_cartesian(const RefSCVector&new_internal)
{
  // get a reference to Molecule for convenience
  Molecule& molecule = *(molecule_.pointer());

  // update the geometry
  for(int i=0; i < dim_.n(); i++) {
    molecule.r(i/3,i%3) = new_internal(i);
  }

  return 0;
}

// again, the coordinates we want to use are cartesian, so just copy
// the cartesian coords into internal
int
CartMolecularCoor::to_internal(RefSCVector&internal)
{
  // get a reference to Molecule for convenience
  Molecule& molecule = *(molecule_.pointer());
  
  int n = dim_.n();
  for (int i=0; i < n; i++) {
    internal(i) = molecule.r(i/3,i%3);
  }

  return 0;
}

int
CartMolecularCoor::to_cartesian(RefSCVector&gradient,RefSCVector&internal)
{
  gradient->assign(internal.pointer());
  return 0;
}

// converts the gradient in cartesian coordinates to internal coordinates
int
CartMolecularCoor::to_internal(RefSCVector&internal,RefSCVector&gradient)
{
  internal->assign(gradient.pointer());
  return 0;
}

int
CartMolecularCoor::to_cartesian(RefSymmSCMatrix&cart,RefSymmSCMatrix&internal)
{
  cart->assign(internal.pointer());
  return 0;
}

int
CartMolecularCoor::to_internal(RefSymmSCMatrix&internal,RefSymmSCMatrix&cart)
{
  internal->assign(cart.pointer());
  return 0;
}

void
CartMolecularCoor::print(ostream& os)
{
  molecule_->print(os);
}

void
CartMolecularCoor::print_simples(ostream& os)
{
}

void
CartMolecularCoor::guess_hessian(RefSymmSCMatrix&hessian)
{
  SymmMolecularCoor imcoor(molecule_);
  RefSymmSCMatrix ihessian(imcoor.dim(),matrixkit_);
  imcoor.guess_hessian(ihessian);
  imcoor.to_cartesian(hessian,ihessian);

  RefSCMatrix evecs(hessian.dim(),hessian.dim(),matrixkit_);
  RefDiagSCMatrix evals(hessian.dim(),matrixkit_);

  hessian.diagonalize(evals,evecs);
  hessian.assign(0.0);

  // get rid of the 3 translations and 3 rotations
  for (int i=0; i < evals.n(); i++) {
    if (fabs(evals.get_element(i)) < 1.0e-6) {
      for (int j=0; j < evals.n(); j++)
        evecs.set_element(j,i,0.0);
      evals.set_element(i,0.0);
    }
  }

  hessian.accumulate_transform(evecs,evals);
}

RefSymmSCMatrix
CartMolecularCoor::inverse_hessian(RefSymmSCMatrix& hessian)
{
  return hessian.gi();
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:

//
// redund.cc
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

#include <util/state/stateio.h>
#include <math/scmat/matrix.h>
#include <math/scmat/local.h>
#include <chemistry/molecule/molecule.h>
#include <chemistry/molecule/coor.h>
#include <chemistry/molecule/simple.h>

#include <util/container/bitarray.h>

using namespace sc;

///////////////////////////////////////////////////////////////////////////
// members of RedundMolecularCoor

static ClassDesc RedundMolecularCoor_cd(
  typeid(RedundMolecularCoor),"RedundMolecularCoor",1,"public IntMolecularCoor",
  0, create<RedundMolecularCoor>, create<RedundMolecularCoor>);

RedundMolecularCoor::RedundMolecularCoor(Ref<Molecule>&mol):
  IntMolecularCoor(mol)
{
  init();
}

RedundMolecularCoor::RedundMolecularCoor(const Ref<KeyVal>& keyval):
  IntMolecularCoor(keyval)
{
  init();
}

RedundMolecularCoor::RedundMolecularCoor(StateIn& s):
  IntMolecularCoor(s)
{
}

RedundMolecularCoor::~RedundMolecularCoor()
{
}

void
RedundMolecularCoor::save_data_state(StateOut&s)
{
  IntMolecularCoor::save_data_state(s);
}

void
RedundMolecularCoor::form_coordinates(int keep_variable)
{
  if (!keep_variable) variable_ = all_;

  if (form_print_simples_) print_simples(ExEnv::out0());
  if (form_print_variable_) print_variable(ExEnv::out0());
  if (form_print_constant_) print_constant(ExEnv::out0());
}

void
RedundMolecularCoor::guess_hessian(RefSymmSCMatrix&hessian)
{
  variable_->guess_hessian(molecule_,hessian);
}

RefSymmSCMatrix
RedundMolecularCoor::inverse_hessian(RefSymmSCMatrix& hessian)
{
  RefSCDimension dredun = hessian.dim();
  
  // form bmat for variable coordinates (ie all the simples)
  RefSCMatrix bmat(dredun,dnatom3_,matrixkit_);
  variable_->bmat(molecule_,bmat);
  
  // and form G = (B*B+)
  RefSymmSCMatrix bmbt(dredun,matrixkit_);
  bmbt.assign(0.0);
  bmbt.accumulate_symmetric_product(bmat);

  // free bmat, and allocate storage for the projection matrix p
  bmat = 0;

  RefSCMatrix p(dredun,dredun,matrixkit_);
  p.assign(0.0);

  // form p = G- * G
  for (int i=0; i < dredun->n(); i++)
    p.set_element(i,i,1.0);
  p = bmbt * p;
  p = bmbt.gi()*p;
  
  // accumulate (p*hessian*p).gi() into bmbt
  bmbt.assign(0.0);
  bmbt.accumulate_transform(p,hessian);
  bmbt = bmbt.gi();
  
  // finally return hinv = p*(p*h*p)-*p
  RefSymmSCMatrix thess = hessian.clone();
  thess.assign(0.0);
  thess.accumulate_transform(p,bmbt);
  return thess;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:

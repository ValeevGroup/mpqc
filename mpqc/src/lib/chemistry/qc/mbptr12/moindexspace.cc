//
// moindexspace.cc
//
// Copyright (C) 2003 Edward Valeev
//
// Author: Edward Valeev <edward.valeev@chemistry.gatech.edu>
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

#ifdef __GNUG__
#pragma implementation
#endif

#include <stdexcept>
#include <stdlib.h>
#include <util/misc/formio.h>
#include <chemistry/qc/basis/petite.h>
#include <chemistry/qc/mbptr12/moindexspace.h>

using namespace std;
using namespace sc;

inline int max(int a,int b) { return (a > b) ? a : b;}

/*---------------
  MOIndexSpace
 ---------------*/
static ClassDesc MOIndexSpace_cd(
  typeid(MOIndexSpace),"MOIndexSpace",1,"virtual public SavableState",
  0, 0, create<MOIndexSpace>);

MOIndexSpace::MOIndexSpace(const RefSCMatrix& coefs, const Ref<GaussianBasisSet> basis, const vector<int>& offsets) :
  coefs_(coefs), basis_(basis), offsets_(offsets)
{
  // Compute MO syms array
  
  init();
}

MOIndexSpace::MOIndexSpace(const RefSCMatrix& coefs, const Ref<GaussianBasisSet> basis, const vector<int>& offsets,
                           const vector<int>& mosym, IndexOrder moorder) :
  coefs_(coefs), basis_(basis), offsets_(offsets), mosym_(mosym)
{
  check_mosym();
  moorder_ = moorder;
  
  init();
}

MOIndexSpace::MOIndexSpace(StateIn& si) : SavableState(si)
{
  coefs_.restore(si);
  basis_ << SavableState::restore_state(si);

  si.get(offsets_);
  si.get(mosym_);

  int moorder; si.get(moorder); moorder = (int) moorder_;
  
  init();
}

MOIndexSpace::~MOIndexSpace()
{
}

void
MOIndexSpace::save_data_state(StateOut& so)
{
  coefs_.save(so);
  SavableState::save_state(basis_.pointer(),so);

  so.put(offsets_);
  so.put(mosym_);

  so.put((int)moorder_);
}


Ref<GaussianBasisSet>
MOIndexSpace::basis() const { return basis_; }

RefSCMatrix
MOIndexSpace::coefs() const { return coefs_; }

vector<int>
MOIndexSpace::mosym() const { return mosym_; }

MOIndexSpace::IndexOrder
MOIndexSpace::moorder() const { return moorder_; }

int
MOIndexSpace::rank() const { return rank_; }

int
MOIndexSpace::full_rank() const { return full_rank_; }

int
MOIndexSpace::nblocks() const { return nblocks_; }

vector<int>
MOIndexSpace::first_mo() const { return first_mo_; }

vector<int>
MOIndexSpace::nmo() const { return nmo_; }

vector<int>
MOIndexSpace::offsets() const { return offsets_; }

void
MOIndexSpace::check_mosym() const
{
  int ng = basis_->molecule()->point_group()->char_table().order();
  
  for(vector<int>::const_iterator p=mosym_.begin(); p != mosym_.end(); ++p) {
    if (*p < 0 || *p >= ng)
      throw std::runtime_error("MOIndexSpace::check_mosym() -- invalid value in the list of orbital irreps");
  }
}

void
MOIndexSpace::init()
{

}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:

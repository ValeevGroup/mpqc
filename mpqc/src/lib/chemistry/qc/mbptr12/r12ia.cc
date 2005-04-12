//
// r12ia.cc
//
// Copyright (C) 2002 Edward Valeev
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

#ifdef __GNUC__
#pragma implementation
#endif

#include <stdexcept>
#include <stdlib.h>
#include <util/misc/formio.h>
#include <util/misc/exenv.h>
#include <chemistry/qc/mbptr12/r12ia.h>

using namespace std;
using namespace sc;

/*--------------------------------
  R12IntsAcc
 --------------------------------*/
static ClassDesc R12IntsAcc_cd(
  typeid(R12IntsAcc),"R12IntsAcc",1,"virtual public SavableState",
  0, 0, 0);

R12IntsAcc::R12IntsAcc(int num_te_types, int ni, int nj, int nx, int ny) :
  num_te_types_(num_te_types), ni_(ni), nj_(nj), nx_(nx), ny_(ny)
{
  nxy_ = nx_*ny_;
  blksize_ = nxy_*sizeof(double);
  blocksize_ = blksize_*num_te_types_;
  committed_ = false;
  active_ = false;
  next_orbital_ = 0;
}

R12IntsAcc::R12IntsAcc(StateIn& si) : SavableState(si)
{
  si.get(num_te_types_);
  si.get(ni_);
  si.get(nj_);
  si.get(nx_);
  si.get(ny_);
  int committed; si.get(committed); committed_ = (bool) committed;
  int active; si.get(active); active_ = (bool) active;
  si.get(next_orbital_);

  nxy_ = nx_ * ny_;
  blksize_ = nxy_*sizeof(double);
  blocksize_ = blksize_*num_te_types_;
}

R12IntsAcc::~R12IntsAcc()
{
  deactivate();
}

void R12IntsAcc::save_data_state(StateOut& so)
{
  so.put(num_te_types_);
  so.put(ni_);
  so.put(nj_);
  so.put(nx_);
  so.put(ny_);
  so.put((int)committed_);
  so.put((int)active_);
  so.put(next_orbital_);
}

void
R12IntsAcc::commit()
{
  if (!committed_)
    committed_ = true;
  else
    throw std::runtime_error("R12IntsAcc::commit() -- accumulator has already been committed");
  activate();
}

void
R12IntsAcc::activate()
{
  if (active_)
    throw std::runtime_error("R12IntsAcc::activate() -- accumulator is already active");
  if (is_committed())
    active_ = true;
  else
    throw std::runtime_error("R12IntsAcc::activate() -- accumulator hasn't been committed yet");
}

void
R12IntsAcc::deactivate()
{
  active_ = false;
}

int R12IntsAcc::next_orbital() const
{
  return next_orbital_;
}

void R12IntsAcc::inc_next_orbital(int ni)
{
  next_orbital_ += ni;
}

int
R12IntsAcc::tasks_with_access(vector<int>& twa_map) const
{
  int nproc = ntasks();
  
  // Compute the number of tasks that have full access to the integrals
  // and split the work among them
  int nproc_with_ints = 0;
  for(int proc=0;proc<nproc;proc++)
    if (has_access(proc)) nproc_with_ints++;
 
  twa_map.resize(nproc);
  int count = 0;
  for(int proc=0;proc<nproc;proc++)
    if (has_access(proc)) {
      twa_map[proc] = count;
      count++;
    }
    else
      twa_map[proc] = -1;

  return nproc_with_ints;
}

///////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:

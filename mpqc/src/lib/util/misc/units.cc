//
// units.cc
//
// Copyright (C) 1997 Limit Point Systems, Inc.
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

#ifdef __GNUG__
#pragma implementation
#endif

#include <math.h>
#include <string.h>
#include <util/misc/units.h>

//////////////////////////////////////////////////////////////////////
// Units class definition

#define CLASSNAME Units
#define PARENTS public SavableState
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
Units::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] =  SavableState::_castdown(cd);
  return do_castdowns(casts,cd);
}

Units::Units(const char *strrep)
{
  if (strrep) strrep_ = strcpy(new char[strlen(strrep)+1], strrep);
  else strrep_ = 0;
  parse_unit();
}

Units::Units(char *strrep, Units::Storage action)
{
  if (action == Copy) {
      if (strrep) strrep_ = strcpy(new char[strlen(strrep)+1], strrep);
      else strrep_ = 0;
    }
  else {
      strrep_ = strrep;
    }
  parse_unit();
}

Units::Units(StateIn&s):
  SavableState(s)
{
  s.getstring(strrep_);
  parse_unit();
}

Units::~Units()
{
  delete[] strrep_;
}

void
Units::save_data_state(StateOut&s)
{
  s.putstring(strrep_);
}

double
Units::to_atomic_units() const
{
  return to_atomic_units_;
}

double
Units::from_atomic_units() const
{
  return 1.0/to_atomic_units_;
}

const char *
Units::string_rep() const
{
  return strrep_;
}

void
Units::parse_unit()
{
  to_atomic_units_ = 1.0;

  if (!strrep_) {
    }
  else if (!strcmp(strrep_, "bohr")
      ||!strcmp(strrep_, "bohrs")) {
    }
  else if (!strcmp(strrep_, "radian")
      ||!strcmp(strrep_, "radians")) {
    }
  else if (!strcmp(strrep_, "angstrom")
      ||!strcmp(strrep_, "angstroms")
      ||!strcmp(strrep_, "aangstrom")
      ||!strcmp(strrep_, "aangstroms")) {
      to_atomic_units_ *= 1.0/0.52917706;
    }
  else if (!strcmp(strrep_, "degree")
      ||!strcmp(strrep_, "degrees")) {
      to_atomic_units_ *= 3.14159265358979323846/180.0;
    }
  else {
      cerr << "Units: Cannot handle \"" << strrep_ << "\"" << endl;
      abort();
    }
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:

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
Units::to(const RefUnits &units) const
{
  return to_atomic_units_/units->to_atomic_units_;
}

double
Units::from(const RefUnits &units) const
{
  return 1.0/to(units);
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

  int invert = 0;
  const char *rest = strrep_;

  while (rest) {
      const char *end = strpbrk(rest, " */");
      int nchar;
      if (end) {
          nchar = end - rest;
        }
      else {
          nchar = strlen(rest);
        }

      // physical constants used for atomic unit conversion factors
      const double a0 = 5.2917706e-11; // m
      const double hbar = 1.0545887e-34; // J s
      const double e = 1.6021892e-19; // C
      const double me = 9.109534e-31; // kg
      const double e0 = 8.854187818e-12; // F/m

      // derived au conversion factors
      const double Ea = e*e/((4.0*M_PI*e0)*a0); // J
      const double time = hbar/Ea; // s

      // other conversions
      const double amu = 1.6605655e-27; // kg

      double factor = 1.0;
      if (!strncmp(rest, "bohr", nchar)
          ||!strncmp(rest, "bohrs", nchar)) {
        }
      else if (!strncmp(rest, "radian", nchar)
               ||!strncmp(rest, "radians", nchar)) {
        }
      else if (!strncmp(rest, "N", nchar)
               ||!strncmp(rest, "newton", nchar)) {
          factor = a0/Ea;
        }
      else if (!strncmp(rest, "dyne", nchar)) {
          factor = 1.0e-5*a0/Ea;
        }
      else if (!strncmp(rest, "m", nchar)
               ||!strncmp(rest, "meter", nchar)) {
          factor = 1.0/a0;
        }
      else if (!strncmp(rest, "cm", nchar)
               ||!strncmp(rest, "centimeter", nchar)) {
          factor = 1.0e-2/a0;
        }
      else if (!strncmp(rest, "angstrom", nchar)
               ||!strncmp(rest, "angstroms", nchar)
               ||!strncmp(rest, "aangstrom", nchar)
               ||!strncmp(rest, "aangstroms", nchar)) {
          factor = 1.0e-10/a0;
        }
      else if (!strncmp(rest, "amu", nchar)) {
          factor = amu/me;
        }
      else if (!strncmp(rest, "degree", nchar)
               ||!strncmp(rest, "degrees", nchar)) {
          factor = M_PI/180.0;
        }
      else {
          cerr << "Units: Cannot handle \"" << rest << "\"" << endl;
          abort();
        }
      if (invert) factor = 1.0/factor;
      to_atomic_units_ *= factor;
      rest = strpbrk(rest, " */");
      while (rest && (*rest == ' ' || *rest == '*' || *rest == '/')) {
          if (*rest == '/') invert = !invert;
          rest++;
        }
    }
  //cout << "FOR " << strrep_ << " got " << to_atomic_units_ << endl;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:

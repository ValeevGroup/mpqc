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

#include <util/misc/math.h>
#include <string.h>
#include <util/misc/units.h>
#include <util/state/stateio.h>

#include <util/state/linkage.h>

using namespace std;
using namespace sc;

namespace sc {

//////////////////////////////////////////////////////////////////////
// Utility functions

static inline int eq(const char* a, const char* b)
{
  return !strcmp(a,b);
}

//////////////////////////////////////////////////////////////////////
// Units class definition

static ClassDesc Units_cd(typeid(Units),"Units",1,"public SavableState",
                          0, 0, create<Units>);

Units::Units(const char *strrep)
{
  if (strrep) {
      if (strrep[0] == '\0') strrep_ = 0;
      else strrep_ = strcpy(new char[strlen(strrep)+1], strrep);
    }
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
Units::to(const Ref<Units> &units) const
{
  return to_atomic_units_/units->to_atomic_units_;
}

double
Units::from(const Ref<Units> &units) const
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
      const char *end = ::strpbrk(rest, " */");
      int nchar;
      if (end) {
          nchar = end - rest;
        }
      else {
          nchar = strlen(rest);
        }
      char *unitstring = new char[nchar+1];
      memcpy(unitstring,rest,nchar);
      unitstring[nchar] = '\0';

      // physical constants used for atomic unit conversion factors
      // from CRC Handbook 77th Ed. Tables 1&2 (CODATA 1986)
      const double a0 = 5.29177249e-11; // m
      //const double hbar = 1.05457266e-34; // J s
      const double e = 1.60217733e-19; // C
      const double me = 9.1093897e-31; // kg
      const double e0 = 8.854187817e-12; // F/m
      const double NA = 6.0221367e23; // mol-1 (from CRC Handbook)

      // derived au conversion factors
      const double Ea = e*e/((4.0*M_PI*e0)*a0); // J
      //const double time = hbar/Ea; // s

      // other conversions
      const double amu = 1.6605655e-27; // kg
      const double cal = 4.184; // J

      double factor = 1.0;
      if (eq(unitstring, "bohr")
          ||eq(unitstring, "bohrs")) {
        } 
      else if (eq(unitstring, "hartree")
               ||eq(unitstring, "hartrees")
               ||eq(unitstring, "Hartree")
               ||eq(unitstring, "Hartrees")) {
        }
      else if (eq(unitstring, "ev")
               ||eq(unitstring, "eV")) {
          factor = 1.0/27.2113834; // physics.nist.gov/constants
        }
      else if (eq(unitstring, "debye")) {
          factor = 1.0/2.541765; // several WWW sources
        }
      else if (eq(unitstring, "radian")
               ||eq(unitstring, "radians")) {
        }
      else if (eq(unitstring, "mol")
               ||eq(unitstring, "mole")) {
          factor = NA;
        }
      else if (eq(unitstring, "kcal")) {
          factor = 1000.0*cal/Ea;     
        }
      else if (eq(unitstring, "kcal_per_mol")) {
          factor = 1000.0*cal/(Ea*NA);     
        }
      else if (eq(unitstring, "N")
               ||eq(unitstring, "newton")) {
          factor = a0/Ea;
        }
      else if (eq(unitstring, "dyne")) {
          factor = 1.0e-5*a0/Ea;
        }
      else if (eq(unitstring, "m")
               ||eq(unitstring, "meter")) {
          factor = 1.0/a0;
        }
      else if (eq(unitstring, "cm")
               ||eq(unitstring, "centimeter")) {
          factor = 1.0e-2/a0;
        }
      else if (eq(unitstring, "angstrom")
               ||eq(unitstring, "angstroms")
               ||eq(unitstring, "aangstrom")
               ||eq(unitstring, "aangstroms")) {
          factor = 1.0e-10/a0;
        }
      else if (eq(unitstring, "amu")) {
          factor = amu/me;
        }
      else if (eq(unitstring, "degree")
               ||eq(unitstring, "degrees")) {
          factor = M_PI/180.0;
        }
      else {
          ExEnv::errn() << "Units: Cannot handle \"" << unitstring << "\"" << endl;
          abort();
        }
      delete[] unitstring;
      if (invert) factor = 1.0/factor;
      to_atomic_units_ *= factor;
      rest = ::strpbrk(rest, " */");
      while (rest && (*rest == ' ' || *rest == '*' || *rest == '/')) {
          if (*rest == '/') invert = !invert;
          rest++;
        }
    }
}

/////////////////////////////////////////////////////////////////////////////

}

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:

//
// molinfo.cc
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include <util/misc/units.h>
#include <util/misc/formio.h>
#include <util/group/message.h>
#include <chemistry/molecule/atominfo.h>

////////////////////////////////////////////////////////////////////////
// AtomInfo

struct AtomInfo::atomname
AtomInfo::names_[MaxZ] = 
  {{"",           ""}, // 0
   {"hydrogen",   "H"}, // 1
   {"helium",     "He"}, // 2
   {"lithium",    "Li"}, // 3
   {"beryllium",  "Be"}, // 4
   {"boron",      "B"}, // 5
   {"carbon",     "C"}, // 6
   {"nitrogen",   "N"}, // 7
   {"oxygen",     "O"}, // 8
   {"fluorine",   "F"}, // 9
   {"neon",       "Ne"}, // 10
   {"sodium",     "Na"}, // 11
   {"magnesium",  "Mg"}, // 12
   {"aluminum",   "Al"}, // 13
   {"silicon",    "Si"}, // 14
   {"phosphorus" ,"P"}, // 15
   {"sulfur",     "S"}, // 16
   {"chlorine",   "Cl"}, // 17
   {"argon",      "Ar"}, // 18
   {"potassium",  "K"}, // 19
   {"calcium",    "Ca"}, // 20
   {"scandium",   "Sc"}, // 21
   {"titanium",   "Ti"}, // 22
   {"vanadium",   "V"}, // 23
   {"chromium",   "Cr"}, // 24
   {"manganese",  "Mn"}, // 25
   {"iron",       "Fe"}, // 26
   {"cobalt",     "Co"}, // 27
   {"nickel",     "Ni"}, // 28
   {"copper",     "Cu"}, // 29
   {"zinc",       "Zn"}, // 30
   {"gallium",    "Ga"}, // 31
   {"germanium",  "Ge"}, // 32
   {"arsenic",    "As"}, // 33
   {"selenium",   "Se"}, // 34
   {"bromine",    "Br"}, // 35
   {"krypton",    "Kr"}, // 36
   {"rubidium",     "Rb"}, // 37
   {"strontium",    "Sr"}, // 38
   {"yttrium",      "Y"}, // 39
   {"zirconium",    "Zr"}, // 40
   {"niobium",      "Nb"}, // 41
   {"molybdenum",   "Mo"}, // 42
   {"technetium",   "Tc"}, // 43
   {"ruthenium",    "Ru"}, // 44
   {"rhodium",      "Rh"}, // 45
   {"palladium",    "Pd"}, // 46
   {"silver",       "Ag"}, // 47
   {"cadminium",    "Cd"}, // 48
   {"indium",       "In"}, // 49
   {"tin",          "Sn"}, // 50
   {"antimony",     "Sb"}, // 51
   {"tellurium",    "Te"}, // 52
   {"iodine",       "I"}, // 53
   {"xenon",        "Xe"}, // 54
   {"cesium",       "Cs"}, // 55
   {"barium",       "Ba"}, // 56
   {"lanthanium",   "La"}, // 57
   {"cerium",       "Ce"}, // 58
   {"praseodymium", "Pr"}, // 59
   {"neodymium",    "Nd"}, // 60
   {"promethium",   "Pm"}, // 61
   {"samarium",     "Sm"}, // 62
   {"europium",     "Eu"}, // 63
   {"gadolinium",   "Gd"}, // 64
   {"terbium",      "Tb"}, // 65
   {"dysprosium",   "Dy"}, // 66
   {"holmium",      "Ho"}, // 67
   {"erbium",       "Er"}, // 68
   {"thulium",      "Tm"}, // 69
   {"ytterbium",    "Yb"}, // 70
   {"lutetium",     "Lu"}, // 71
   {"hafnium",      "Hf"}, // 72
   {"tantalum",     "Ta"}, // 73
   {"tungsten",     "W"}, // 74
   {"rhenium",      "Re"}, // 75
   {"osmium",       "Os"}, // 76
   {"iridium",      "Ir"}, // 77
   {"platinum",     "Pt"}, // 78
   {"gold",         "Au"}, // 79
   {"mercury",      "Hg"}, // 80
   {"thallium",     "Tl"}, // 81
   {"lead",         "Pb"}, // 82
   {"bismuth",      "Bi"}, // 83
   {"polonium",     "Po"}, // 84
   {"astatine",     "At"}, // 85
   {"radon",        "Rn"}, // 86
   {"francium",     "Fr"}, // 87
   {"radium",       "Ra"}, // 88
   {"actinium",     "Ac"}, // 89
   {"thorium",      "Th"}, // 90
   {"protactinium", "Pa"}, // 91
   {"uranium",      "U"}, // 92
   {"neptunium",    "Np"}, // 93
   {"plutonium",    "Pu"}, // 94
   {"americium",    "Am"}, // 95
   {"curium",       "Cm"}, // 96
   {"berkelium",    "Bk"}, // 97
   {"californium",  "Cf"}, // 98
   {"einsteinum",   "Es"}, // 99
   {"fermium",      "Fm"}, // 100
   {"mendelevium",  "Md"}, // 101
   {"nobelium",     "No"}, // 102
   {"lawrencium",   "Lr"}, // 103
   {"rutherfordium","Rf"}, // 104
   {"hahnium",      "Ha"}, // 105
   {"Unnamed",      "Un"}, // 106
   {"Unnamed",      "Un"} // 107
  };

#define CLASSNAME AtomInfo
#define PARENTS public SavableState
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
AtomInfo::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SavableState::_castdown(cd);
  return do_castdowns(casts,cd);
}

AtomInfo::AtomInfo()
{
  overridden_values_ = 0;
  load_library_values();
}

AtomInfo::AtomInfo(const RefKeyVal& keyval)
{
  overridden_values_ = 0;
  load_library_values();
  override_library_values(keyval);
}

AtomInfo::AtomInfo(StateIn& s):
  SavableState(s)
{
  if (s.node_to_node()) {
      s.get_array_double(mass_,MaxZ);
      s.get_array_double(atomic_radius_,MaxZ);
      s.get_array_double(vdw_radius_,MaxZ);
      s.get_array_double(bragg_radius_,MaxZ);
      for (int i=0; i<MaxZ; i++) s.get_array_double(rgb_[i],3);
      s.getstring(overridden_values_);
    }
  else {
      overridden_values_ = 0;
      load_library_values();
      char *overrides;
      s.getstring(overrides);
      if (overrides) {
          RefParsedKeyVal keyval = new ParsedKeyVal;
          keyval->parse_string(overrides);
          override_library_values(keyval.pointer());
          delete[] overrides;
        }
    }
}

AtomInfo::~AtomInfo()
{
  delete[] overridden_values_;
}

void
AtomInfo::save_data_state(StateOut& s)
{
  if (s.node_to_node()) {
      s.put_array_double(mass_,MaxZ);
      s.put_array_double(atomic_radius_,MaxZ);
      s.put_array_double(vdw_radius_,MaxZ);
      s.put_array_double(bragg_radius_,MaxZ);
      for (int i=0; i<MaxZ; i++) s.put_array_double(rgb_[i],3);
      s.putstring(overridden_values_);
    }
  else {
      s.putstring(overridden_values_);
    }
}

void
AtomInfo::load_library_values()
{
  RefMessageGrp grp = MessageGrp::get_default_messagegrp();
  if (grp->me() == 0) {
      const char* libdir;
      RefKeyVal keyval;
      if (libdir = getenv("SCLIBDIR")) {
          const char* atominfo = "/atominfo.kv";
          const char *eq = strchr(libdir,'=');
          if (eq) libdir = eq + 1;
          char *filename = new char[strlen(libdir) + strlen(atominfo) + 1];
          strcpy(filename, libdir);
          strcat(filename, atominfo);
          keyval = new ParsedKeyVal(filename);
          delete[] filename;
        }
      else {
          keyval = new ParsedKeyVal(SRCLIBDIR "atominfo.kv");
        }
      RefKeyVal pkeyval = new PrefixKeyVal(keyval, "atominfo");
      load_values(pkeyval,0);
    }
  grp->bcast(mass_,MaxZ);
  grp->bcast(atomic_radius_,MaxZ);
  grp->bcast(vdw_radius_,MaxZ);
  grp->bcast(bragg_radius_,MaxZ);
  for (int i=0; i<MaxZ; i++) grp->bcast(rgb_[i],3);
}

void
AtomInfo::override_library_values(const RefKeyVal &keyval)
{
  load_values(keyval, 1);
}

void
AtomInfo::load_values(const RefKeyVal& keyval, int override)
{
  int i;

  RefUnits amu = new Units("amu");
  RefUnits bohr = new Units("bohr");

  load_values(mass_, "mass", keyval, override, amu);
  load_values(atomic_radius_, "atomic_radius", keyval, override, bohr);
  load_values(vdw_radius_, "vdw_radius", keyval, override, bohr);
  load_values(bragg_radius_, "bragg_radius", keyval, override, bohr);
  load_values(rgb_, "rgb", keyval, override);
}

void
AtomInfo::load_values(double *array, const char *keyword,
                      const RefKeyVal &keyval, int override,
                      const RefUnits &units)
{
  RefKeyVal pkeyval = new PrefixKeyVal(keyval,keyword);
  RefUnits fileunits = new Units(pkeyval->pcharvalue("unit"), Units::Steal);
  double f = 1.0;
  if (fileunits.nonnull() && units.nonnull()) {
      f = fileunits->to(units);
    }
  double def;
  if (!override) {
      def = pkeyval->doublevalue("default");
      array[0] = def;
    }
  int have_overridden = 0;
  int i;
  for (i=1; i<MaxZ; i++) {
      double val = f * pkeyval->doublevalue(names_[i].symbol);
      if (pkeyval->error() != KeyVal::OK) {
          if (!override) array[i] = def;
        }
      else {
          array[i] = val;
          if (override) {
              const char *prefix = " ";
              if (!have_overridden) {
                  add_overridden_value(keyword);
                  add_overridden_value(":(");
                  if (fileunits.nonnull() && fileunits->string_rep()) {
                      char ustring[256];
                      sprintf(ustring,"unit=\"%s\"",fileunits->string_rep());
                      add_overridden_value(ustring);
                    }
                  else {
                      prefix = "";
                    }
                  have_overridden = 1;
                }
              char *strval = pkeyval->pcharvalue(names_[i].symbol);
              char assignment[256];
              sprintf(assignment,"%s%s=%s", prefix, names_[i].symbol, strval);
              delete[] strval;
              add_overridden_value(assignment);
            }
        }
    }
  if (have_overridden) {
      add_overridden_value(")");
    }
}

void
AtomInfo::load_values(double array[][3], const char *keyword,
                      const RefKeyVal &keyval, int override)
{
  int i,j;
  RefKeyVal pkeyval = new PrefixKeyVal(keyval,keyword);
  double def[3];
  if (!override) {
      for (i=0; i<3; i++) {
          def[i] = pkeyval->doublevalue("default",i);
          array[0][i] = def[i];
        }
    }
  int have_overridden = 0;
  for (i=1; i<MaxZ; i++) {
      double val;
      for (j=0; j<3; j++) {
          val = pkeyval->doublevalue(names_[i].symbol,j);
          if (pkeyval->error() != KeyVal::OK) {
              if (!override) array[i][j] = def[j];
            }
          else {
              array[i][j] = val;
              if (override) {
                  const char *prefix = " ";
                  if (!have_overridden) {
                      add_overridden_value(keyword);
                      add_overridden_value(":(");
                      prefix = "";
                      have_overridden = 1;
                    }
                  char *strval = pkeyval->pcharvalue(names_[i].symbol,j);
                  char assignment[256];
                  sprintf(assignment,"%s%s:%d=%s",
                          prefix, names_[i].symbol, j, strval);
                  delete[] strval;
                  add_overridden_value(assignment);
                }
            }
        }
    }
  if (have_overridden) {
      add_overridden_value(")");
    }
}

void
AtomInfo::add_overridden_value(const char *assignment)
{
  int length = strlen(assignment)+1;
  if (overridden_values_) length += strlen(overridden_values_);
  char *new_overridden_values = new char[length];
  new_overridden_values[0] = '\0';
  if (overridden_values_) strcat(new_overridden_values, overridden_values_);
  strcat(new_overridden_values, assignment);
  delete[] overridden_values_;
  overridden_values_ = new_overridden_values;
}

int
AtomInfo::string_to_Z(const char *name)
{
  int i,j;
  int Z;

  // see if the name is a atomic number
  Z = atoi(name);

  // name is not an atomic number--must be a symbol or atom name
  if (!Z) {

      // convert the name to lower case
      char* tmpname = strdup(name);
      for (j=0; j<strlen(tmpname); j++) {
	  if (isupper(tmpname[j])) tmpname[j] = tolower(tmpname[j]);
	}

      // loop thru the elements, looking for a match
      for (i=0; i<MaxZ; i++) {

          // see if an atom name matches
          // if we find a match we are done looping
	  if (!strcmp(names_[i].name,tmpname)) {
	      Z = i;
	      break;
	    }

          // see if an atomic symbol (converted to lower case) matches
          char* tmpsymbol = strdup(names_[i].symbol);
          for (j=0; j<strlen(tmpsymbol); j++) {
	      if (isupper(tmpsymbol[j])) tmpsymbol[j] = tolower(tmpsymbol[j]);
            }
          // if we find a match we are done looping
	  if (!strcmp(tmpsymbol,tmpname)) {
	      Z = i;
              free(tmpsymbol);
	      break;
	    }
          free(tmpsymbol);
	}

      // free the lowercase version of the name
      free(tmpname);
    }  

  // check to see if z value is OK, if not then the name must have been
  // invalid
  if (Z < 1 || Z > MaxZ) {
      cerr << node0 << sprintf("AtomInfo: invalid name: %s\n",name);
      exit(1);
    }

  return Z;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:

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
#include <util/misc/formio.h>
#include <chemistry/molecule/molinfo.h>

////////////////////////////////////////////////////////////////////////
// MolInfo

#define CLASSNAME MolInfo
#define PARENTS public DescribedClass
#define HAVE_KEYVAL_CTOR
#include <util/class/classi.h>
void *
MolInfo::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = DescribedClass::_castdown(cd);
  return do_castdowns(casts,cd);
}

MolInfo::MolInfo(const RefKeyVal& pkeyval)
{
  const char* libdir;
  if (pkeyval->exists("molinfofiles")) {
      RefKeyVal libkeyval = new ParsedKeyVal("molinfo",pkeyval);
      keyval = new AggregateKeyVal(pkeyval,libkeyval);
    }
  else if (libdir = getenv("SCLIBDIR")) {
      const char* molinfo = "/molinfo.ipv2";
      const char *eq = strchr(libdir,'=');
      if (eq) libdir = eq + 1;
      char *filename = new char[strlen(libdir) + strlen(molinfo) + 1];
      strcpy(filename, libdir);
      strcat(filename, molinfo);
      keyval = new ParsedKeyVal(filename);
      delete[] filename;
    }
  else {
      keyval = new ParsedKeyVal(SRCLIBDIR "molinfo.ipv2");
    }
}

MolInfo::~MolInfo()
{
}

KeyVal::KeyValError
MolInfo::error()
{
  return keyval->error();
}

double
MolInfo::doublevalue(const char* sym,
                     const char* property)
{
  char keyword[KeyVal::MaxKeywordLength];
  sprintf(keyword, "%s:%s", sym, property);
  double val = keyval->doublevalue(keyword);
  if (keyval->error() != KeyVal::OK) {
      sprintf(keyword, "default:%s", property);
      val = keyval->doublevalue(keyword);
    }
  return val;
}

double
MolInfo::doublevalue(const char* sym,
                     const char* property, int i)
{
  char keyword[KeyVal::MaxKeywordLength];
  sprintf(keyword, "%s:%s:%d", sym, property, i);
  double val = keyval->doublevalue(keyword);
  if (keyval->error() != KeyVal::OK) {
      sprintf(keyword, "default:%s:%d", property, i);
      val = keyval->doublevalue(keyword);
    }
  return val;
}

////////////////////////////////////////////////////////////////////////
// AtomInfo

#define CLASSNAME AtomInfo
#define PARENTS public DescribedClass
#define HAVE_KEYVAL_CTOR
#include <util/class/classi.h>
void *
AtomInfo::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = MolInfo::_castdown(cd);
  return do_castdowns(casts,cd);
}

AtomInfo::AtomInfo(const RefKeyVal& pkeyval):
  MolInfo(pkeyval)
{
  keyval = new PrefixKeyVal(":molinfo:atoms", keyval);
  for (int i=0; i<ATOMINFO_MAXZ; i++) {
      have_radius[i] = have_rgb[i] = 0;
      rgb_vals[i][0] = rgb_vals[i][1] = rgb_vals[i][2] = 0.0;
      radius_vals[i] = 0.0;
    }
  preload_values();

  // see if pkeyval has a radius scale factor that should override
  // the one in the database
  double tmp = pkeyval->doublevalue("radius_scale_factor");
  if (pkeyval->error() == KeyVal::OK) radius_scale_factor_ = tmp;
}

void
AtomInfo::preload_values()
{
  // radius unit
  char* radius_unit = keyval->pcharvalue("radius_unit");
  if (radius_unit) {
      if (!strcmp(radius_unit, "angstrom")) {
          radius_to_bohr_ = 1.0/0.52917706;
        }
      else if (!strcmp(radius_unit, "bohr")) {
          radius_to_bohr_ = 1.0;
        }
      else {
          cerr << node0 << indent
               << scprintf("AtomInfo: unknown radius_unit: \"%s\"\n",
                           radius_unit);
          abort();
        }
      delete[] radius_unit;
    }
  else {
      radius_to_bohr_ = 1.0;
    }

  // radius scale factor
  radius_scale_factor_ = keyval->doublevalue("radius_scale_factor");
  if (keyval->error() != KeyVal::OK) radius_scale_factor_ = 1.0;
}

AtomInfo::~AtomInfo()
{
}

double
AtomInfo::doublevalue(const ChemicalElement& atom,
                      const char* property)
{
  return MolInfo::doublevalue(atom.symbol(), property);
}

double
AtomInfo::doublevalue(const ChemicalElement& atom,
                      const char* property, int i)
{
  return MolInfo::doublevalue(atom.symbol(), property, i);
}

int
AtomInfo::get_zindex(const ChemicalElement& atom)
{
  int Zindex = atom.number() - 1;
  if (Zindex >= ATOMINFO_MAXZ) {
      cerr << node0 << indent << "AtomInfo: ATOMINFO_MAXZ exceeded\n";
      abort();
    }
  return Zindex;
}

double
AtomInfo::radius(const ChemicalElement& atom)
{
  int zindex = get_zindex(atom);
  if (!have_radius[zindex]) {
      radius_vals[zindex] = doublevalue(atom, "radius")
                          * radius_to_bohr_ * radius_scale_factor_;
      if (error() != KeyVal::OK) {
          cerr << node0 << indent << "AtomInfo couldn't find a radius\n";
          keyval->dump();
          abort();
        }
      if (radius_vals[zindex] <= 1.0e-6) {
          cerr << node0 << indent << "AtomInfo: got a tiny radius\n";
          abort();
        }
      have_radius[zindex] = 1;
    }
  return radius_vals[zindex];
}

double
AtomInfo::rgb(const ChemicalElement& atom, int color)
{
  int zindex = get_zindex(atom);
  if (!have_rgb[zindex]) {
      rgb_vals[zindex][0] = doublevalue(atom, "rgb", 0);
      rgb_vals[zindex][1] = doublevalue(atom, "rgb", 1);
      rgb_vals[zindex][2] = doublevalue(atom, "rgb", 2);
      have_rgb[zindex] = 1;
    }
  return rgb_vals[zindex][color];
}

double
AtomInfo::red(const ChemicalElement& atom)
{
  return rgb(atom, 0);
}

double
AtomInfo::green(const ChemicalElement& atom)
{
  return rgb(atom, 1);
}

double
AtomInfo::blue(const ChemicalElement& atom)
{
  return rgb(atom, 2);
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:

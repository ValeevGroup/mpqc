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

#include <stdlib.h>
#include <util/misc/string.h>
#include <ctype.h>
#include <sys/stat.h>

#include <sstream>
#include <stdexcept>

#include <util/misc/units.h>
#include <util/misc/autovec.h>
#include <util/misc/formio.h>
#include <util/state/stateio.h>
#include <util/group/message.h>
#include <chemistry/molecule/atominfo.h>

using namespace std;
using namespace sc;

////////////////////////////////////////////////////////////////////////
// AtomInfo

struct AtomInfo::atom
AtomInfo::elements_[Nelement] =
  {{1, "hydrogen",   "H"},
   {2, "helium",     "He"},
   {3, "lithium",    "Li"},
   {4, "beryllium",  "Be"},
   {5, "boron",      "B"},
   {6, "carbon",     "C"},
   {7, "nitrogen",   "N"},
   {8, "oxygen",     "O"},
   {9, "fluorine",   "F"},
   {10, "neon",       "Ne"},
   {11, "sodium",     "Na"},
   {12, "magnesium",  "Mg"},
   {13, "aluminum",   "Al"},
   {14, "silicon",    "Si"},
   {15, "phosphorus" ,"P"},
   {16, "sulfur",     "S"},
   {17, "chlorine",   "Cl"},
   {18, "argon",      "Ar"},
   {19, "potassium",  "K"},
   {20, "calcium",    "Ca"},
   {21, "scandium",   "Sc"},
   {22, "titanium",   "Ti"},
   {23, "vanadium",   "V"},
   {24, "chromium",   "Cr"},
   {25, "manganese",  "Mn"},
   {26, "iron",       "Fe"},
   {27, "cobalt",     "Co"},
   {28, "nickel",     "Ni"},
   {29, "copper",     "Cu"},
   {30, "zinc",       "Zn"},
   {31, "gallium",    "Ga"},
   {32, "germanium",  "Ge"},
   {33, "arsenic",    "As"},
   {34, "selenium",   "Se"},
   {35, "bromine",    "Br"},
   {36, "krypton",    "Kr"},
   {37, "rubidium",     "Rb"},
   {38, "strontium",    "Sr"},
   {39, "yttrium",      "Y"},
   {40, "zirconium",    "Zr"},
   {41, "niobium",      "Nb"},
   {42, "molybdenum",   "Mo"},
   {43, "technetium",   "Tc"},
   {44, "ruthenium",    "Ru"},
   {45, "rhodium",      "Rh"},
   {46, "palladium",    "Pd"},
   {47, "silver",       "Ag"},
   {48, "cadminium",    "Cd"},
   {49, "indium",       "In"},
   {50, "tin",          "Sn"},
   {51, "antimony",     "Sb"},
   {52, "tellurium",    "Te"},
   {53, "iodine",       "I"},
   {54, "xenon",        "Xe"},
   {55, "cesium",       "Cs"},
   {56, "barium",       "Ba"},
   {57, "lanthanium",   "La"},
   {58, "cerium",       "Ce"},
   {59, "praseodymium", "Pr"},
   {60, "neodymium",    "Nd"},
   {61, "promethium",   "Pm"},
   {62, "samarium",     "Sm"},
   {63, "europium",     "Eu"},
   {64, "gadolinium",   "Gd"},
   {65, "terbium",      "Tb"},
   {66, "dysprosium",   "Dy"},
   {67, "holmium",      "Ho"},
   {68, "erbium",       "Er"},
   {69, "thulium",      "Tm"},
   {70, "ytterbium",    "Yb"},
   {71, "lutetium",     "Lu"},
   {72, "hafnium",      "Hf"},
   {73, "tantalum",     "Ta"},
   {74, "tungsten",     "W"},
   {75, "rhenium",      "Re"},
   {76, "osmium",       "Os"},
   {77, "iridium",      "Ir"},
   {78, "platinum",     "Pt"},
   {79, "gold",         "Au"},
   {80, "mercury",      "Hg"},
   {81, "thallium",     "Tl"},
   {82, "lead",         "Pb"},
   {83, "bismuth",      "Bi"},
   {84, "polonium",     "Po"},
   {85, "astatine",     "At"},
   {86, "radon",        "Rn"},
   {87, "francium",     "Fr"},
   {88, "radium",       "Ra"},
   {89, "actinium",     "Ac"},
   {90, "thorium",      "Th"},
   {91, "protactinium", "Pa"},
   {92, "uranium",      "U"},
   {93, "neptunium",    "Np"},
   {94, "plutonium",    "Pu"},
   {95, "americium",    "Am"},
   {96, "curium",       "Cm"},
   {97, "berkelium",    "Bk"},
   {98, "californium",  "Cf"},
   {99, "einsteinum",   "Es"},
   {100, "fermium",      "Fm"},
   {101, "mendelevium",  "Md"},
   {102, "nobelium",     "No"},
   {103, "lawrencium",   "Lr"},
   {104, "rutherfordium","Rf"},
   {105, "hahnium",      "Ha"},
   {106, "seaborgium",   "Sg"},
   {107, "bohrium",      "Bh"},
   {108, "hassium",      "Hs"},
   {109, "meitnerium",   "Mt"},
   {110, "darmstadtium", "Ds"},
   {111, "roentgenium",  "Rg"},
   {112, "ununbium",     "Uub"},
   {113, "ununtrium",    "Uut"},
   {114, "ununquadium",  "Uuq"},
   {115, "ununpentium",  "Uup"},
   {116, "ununhexium",   "Uuh"},
   {117, "ununseptium",  "Uus"},
   {118, "ununoctium",   "Uuo"}
  };

static ClassDesc AtomInfo_cd(
  typeid(AtomInfo),"AtomInfo",3,"public SavableState",
  0, create<AtomInfo>, create<AtomInfo>);

AtomInfo::AtomInfo()
{
  initialize_names();
  overridden_values_ = 0;
  load_library_values();
}

AtomInfo::AtomInfo(const Ref<KeyVal>& keyval)
{
  initialize_names();
  overridden_values_ = 0;
  load_library_values();
  override_library_values(keyval);
}

AtomInfo::AtomInfo(StateIn& s):
  SavableState(s)
{
  initialize_names();
  overridden_values_ = 0;
  load_library_values();
  char *overrides;
  s.getstring(overrides);
  if (overrides) {
      Ref<ParsedKeyVal> keyval = new ParsedKeyVal;
      keyval->parse_string(overrides);
      override_library_values(keyval.pointer());
      delete[] overrides;
    }
  if (s.version(::class_desc<AtomInfo>()) < 2) {
      atomic_radius_scale_ = 1.0;
      vdw_radius_scale_ = 1.0;
      bragg_radius_scale_ = 1.0;
      maxprob_radius_scale_ = 1.0;
    }
  else {
      s.get(atomic_radius_scale_);
      s.get(vdw_radius_scale_);
      s.get(bragg_radius_scale_);
      s.get(maxprob_radius_scale_);
    }
}

AtomInfo::~AtomInfo()
{
  delete[] overridden_values_;
}

void
AtomInfo::save_data_state(StateOut& s)
{
  s.putstring(overridden_values_);
  s.put(atomic_radius_scale_);
  s.put(vdw_radius_scale_);
  s.put(bragg_radius_scale_);
  s.put(maxprob_radius_scale_);
}

void
AtomInfo::initialize_names()
{
  for (int i=0; i<Nelement; i++) {
      symbol_to_Z_[elements_[i].symbol] = elements_[i].Z;
      name_to_Z_[elements_[i].name] = elements_[i].Z;
    }

  // Z == DefaultZ is reserved for default values
  symbol_to_Z_["Def"] = DefaultZ;
  name_to_Z_["default"] = DefaultZ;

  // Z == 1000 is reserved for point charges
  name_to_Z_["charge"] = 1000;
  symbol_to_Z_["Q"] = 1000;

  for (std::map<std::string,int>::iterator i = symbol_to_Z_.begin();
       i != symbol_to_Z_.end(); i++) {
      Z_to_symbols_[i->second] = i->first;
    }

  for (std::map<std::string,int>::iterator i = name_to_Z_.begin();
       i != name_to_Z_.end(); i++) {
      Z_to_names_[i->second] = i->first;
    }
}

bool
AtomInfo::has_announced_library_source_ = false;

void
AtomInfo::load_library_values()
{
  Ref<MessageGrp> grp = MessageGrp::get_default_messagegrp();
  sc::auto_vec<char> in_char_array;
  if (grp->me() == 0) {
      const char* libdir;
      std::string filename;
      if ((libdir = getenv("MPQC_DATA_PATH")) != 0) {
          const char* atominfo = "/atominfo.kv";
          const char *eq = ::strchr(libdir,'=');
          if (eq) libdir = eq + 1;
          filename = std::string(libdir) + atominfo;
        }
      else {
          struct stat sb;
          const char *ainfo = MPQCDATAPATH "/atominfo.kv";
#ifdef SRC_MPQC_DATA_PATH
          if (stat(ainfo, &sb) != 0) {
              ainfo = SRC_MPQC_DATA_PATH "/atominfo.kv";
            }
#endif
          filename = ainfo;
        }
      if (!has_announced_library_source_) {
        ExEnv::out0() << indent << "Reading file " << filename << "." << endl;
        has_announced_library_source_ = true;
      }
      ifstream is(filename.c_str());
      ostringstream ostrs;
      is >> ostrs.rdbuf();
      int n = 1 + strlen(ostrs.str().c_str());
      in_char_array.reset(strcpy(new char[n],ostrs.str().c_str()));
      grp->bcast(n);
      grp->bcast(in_char_array.get(), n);
    }
  else {
      int n;
      grp->bcast(n);
      in_char_array.reset(new char[n]);
      grp->bcast(in_char_array.get(), n);
    }

  Ref<ParsedKeyVal> keyval = new ParsedKeyVal();
  keyval->parse_string(in_char_array.get());
  Ref<KeyVal> pkeyval = new PrefixKeyVal(keyval.pointer(), "atominfo");
  load_values(pkeyval,0);
}

void
AtomInfo::override_library_values(const Ref<KeyVal> &keyval)
{
  load_values(keyval, 1);
}

void
AtomInfo::load_values(const Ref<KeyVal>& keyval, int override)
{
  Ref<Units> amu = new Units("amu");
  Ref<Units> bohr = new Units("bohr");
  Ref<Units> hartree = new Units("hartree");

  load_values(Z_to_mass_, 0, "mass", keyval, override, amu);
  load_values(Z_to_atomic_radius_, &atomic_radius_scale_, "atomic_radius",
              keyval, override, bohr);
  load_values(Z_to_vdw_radius_, &vdw_radius_scale_, "vdw_radius",
              keyval, override, bohr);
  load_values(Z_to_bragg_radius_, &bragg_radius_scale_,
              "bragg_radius", keyval, override, bohr);
  load_values(Z_to_maxprob_radius_, &maxprob_radius_scale_,
              "maxprob_radius", keyval, override, bohr);
  load_values(Z_to_rgb_, "rgb", keyval, override);
  load_values(Z_to_ip_, 0, "ip", keyval, override, hartree);
}

void
AtomInfo::load_values(std::map<int,double>&values,
                      double *scale, const char *keyword,
                      const Ref<KeyVal> &keyval, int override,
                      const Ref<Units> &units)
{
  Ref<KeyVal> pkeyval = new PrefixKeyVal(keyval,keyword);
  std::string unit = pkeyval->stringvalue("unit");
  Ref<Units> fileunits = new Units(unit.c_str());
  double f = 1.0;
  if (fileunits && units) {
      f = fileunits->to(units);
    }
  double def = 0.0;
  if (!override) {
      def = pkeyval->doublevalue("default");
      values[DefaultZ] = def;
    }
  int have_overridden = 0;
  for (int elem=0; elem<Nelement; elem++) {
      int Z = elements_[elem].Z;
      double val = f * pkeyval->doublevalue(elements_[elem].symbol);
      if (pkeyval->error() != KeyVal::OK) {
          if (!override) values[Z] = def;
        }
      else {
          values[Z] = val;
          if (override) {
              const char *prefix = " ";
              if (!have_overridden) {
                  add_overridden_value(keyword);
                  add_overridden_value(":(");
                  if (fileunits && fileunits->string_rep()) {
                      char ustring[256];
                      sprintf(ustring,"unit=\"%s\"",fileunits->string_rep());
                      add_overridden_value(ustring);
                    }
                  else {
                      prefix = "";
                    }
                  have_overridden = 1;
                }
              std::string assignment;
              assignment += prefix;
              assignment += elements_[elem].symbol;
              assignment += '=';
              assignment += pkeyval->stringvalue(elements_[elem].symbol);
              add_overridden_value(assignment.c_str());
            }
        }
    }
  if (scale) {
      KeyValValuedouble kvvscale(1.0);
      *scale = pkeyval->doublevalue("scale_factor", kvvscale);
      if (pkeyval->error() == KeyVal::OK) {
          if (override) {
              const char *prefix = " ";
              if (!have_overridden) {
                  add_overridden_value(keyword);
                  add_overridden_value(":(");
                  have_overridden = 1;
                  prefix = "";
                }
              std::string assignment;
              assignment += prefix;
              assignment += "scale_factor=";
              assignment += pkeyval->stringvalue("scale_factor");
              add_overridden_value(assignment.c_str());
            }
        }
    }
  if (have_overridden) {
      add_overridden_value(")");
    }
}

void
AtomInfo::load_values(std::map<int,std::vector<double> >&values,
                      const char *keyword,
                      const Ref<KeyVal> &keyval, int override)
{
  Ref<KeyVal> pkeyval = new PrefixKeyVal(keyval,keyword);
  double def[3];
  if (!override) {
      values[DefaultZ].resize(3);
      for (int i=0; i<3; i++) {
          def[i] = pkeyval->doublevalue("default",i);
          values[DefaultZ][i] = def[i];
        }
    }
  int have_overridden = 0;
  for (int elem=0; elem<Nelement; elem++) {
      double val;
      int Z = elements_[elem].Z;
      values[Z].resize(3);
      for (int j=0; j<3; j++) {
          val = pkeyval->doublevalue(elements_[elem].symbol,j);
          if (pkeyval->error() != KeyVal::OK) {
              if (!override) values[Z][j] = def[j];
            }
          else {
              values[Z][j] = val;
              if (override) {
                  const char *prefix = " ";
                  if (!have_overridden) {
                      add_overridden_value(keyword);
                      add_overridden_value(":(");
                      prefix = "";
                      have_overridden = 1;
                    }
                  std::string strval
                      = pkeyval->stringvalue(elements_[elem].symbol,j);
                  char assignment[256];
                  sprintf(assignment,"%s%s:%d=%s",
                          prefix, elements_[elem].symbol, j, strval.c_str());
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
AtomInfo::string_to_Z(const std::string &name, int allow_exceptions)
{
  int Z;

  // see if the name is a atomic number
  Z = atoi(name.c_str());
  if (Z) return Z;

  // convert the name to lower case
  std::string tmpname(name);
  for (int j=0; j<tmpname.size(); j++) {
      if (isupper(tmpname[j])) tmpname[j] = tolower(tmpname[j]);
    }

  std::map<std::string,int>::const_iterator iname
      = name_to_Z_.find(tmpname);
  if (iname != name_to_Z_.end()) return iname->second;

  if (tmpname.size() > 0) {
      if (islower(tmpname[0])) tmpname[0] = toupper(tmpname[0]);
    }

  iname = symbol_to_Z_.find(tmpname);
  if (iname != symbol_to_Z_.end()) return iname->second;

  if (allow_exceptions) {
      ExEnv::err0() << scprintf("AtomInfo: invalid name: %s\n",name.c_str());
      throw std::runtime_error("invalid atom name");
    }

  return 0;
}

double
AtomInfo::lookup_value(const std::map<int,double>& values, int Z) const
{
  std::map<int,double>::const_iterator found = values.find(Z);
  if (found == values.end()) {
      found = values.find(DefaultZ);
    }
  return found->second;
}

double
AtomInfo::lookup_array_value(const std::map<int,std::vector<double> >& values,
                             int Z, int i) const
{
  std::map<int,std::vector<double> >::const_iterator found = values.find(Z);
  if (found == values.end()) {
      found = values.find(DefaultZ);
    }
  return found->second[i];
}

double
AtomInfo::vdw_radius(int Z) const
{
  return lookup_value(Z_to_vdw_radius_,Z)*vdw_radius_scale_;
}

double
AtomInfo::bragg_radius(int Z) const
{
  return lookup_value(Z_to_bragg_radius_,Z)*bragg_radius_scale_;
}

double
AtomInfo::atomic_radius(int Z) const
{
  return lookup_value(Z_to_atomic_radius_,Z)*atomic_radius_scale_;
}

double
AtomInfo::maxprob_radius(int Z) const
{
  return lookup_value(Z_to_maxprob_radius_,Z)*maxprob_radius_scale_;
}

double
AtomInfo::ip(int Z) const
{
  return lookup_value(Z_to_ip_,Z);
}

double
AtomInfo::rgb(int Z, int color) const
{
  return lookup_array_value(Z_to_rgb_,Z,color);
}

double
AtomInfo::red(int Z) const
{
  return lookup_array_value(Z_to_rgb_,Z,0);
}

double
AtomInfo::green(int Z) const
{
  return lookup_array_value(Z_to_rgb_,Z,1);
}

double
AtomInfo::blue(int Z) const
{
  return lookup_array_value(Z_to_rgb_,Z,2);
}

double
AtomInfo::mass(int Z) const
{
  return lookup_value(Z_to_mass_,Z);
}

std::string
AtomInfo::name(int Z) const
{
  std::map<int,std::string>::const_iterator found = Z_to_names_.find(Z);
  if (found == Z_to_names_.end()) {
      ExEnv::err0() << scprintf("AtomInfo: invalid Z: %d\n",Z);
      throw std::runtime_error("invalid Z");
    }
  return found->second;
}

std::string
AtomInfo::symbol(int Z) const
{
  std::map<int,std::string>::const_iterator found = Z_to_symbols_.find(Z);
  if (found == Z_to_symbols_.end()) {
      ExEnv::err0() << scprintf("AtomInfo: invalid Z: %d\n",Z);
      throw std::runtime_error("invalid Z");
    }
  return found->second;
}

void
AtomInfo::print(std::ostream& os) const {
  os << indent << "AtomInfo:" << endl << incindent;
  os << indent << "vdw_radius_scale = " << vdw_radius_scale_ << endl;
  os << indent << "bragg_radius_scale = " << bragg_radius_scale_ << endl;
  os << indent << "atomic_radius_scale = " << atomic_radius_scale_ << endl;
  os << indent << "maxprob_radius_scale = " << maxprob_radius_scale_ << endl;
  for(int Z=1; Z<=Nelement; ++Z) {
    os << indent << "atomic number " << Z << ":" << endl << incindent;

    os << indent << "symbol               = " << symbol(Z) << endl;
    os << indent << "name                 = " << name(Z) << endl;
    os << indent << "mass                 = " << mass(Z) << endl;
    os << indent << "radius (VdW)         = " << vdw_radius(Z) << endl;
    os << indent << "radius (Bragg)       = " << bragg_radius(Z) << endl;
    os << indent << "radius (atomic)      = " << atomic_radius(Z) << endl;
    os << indent << "radius (maxprob)     = " << maxprob_radius(Z) << endl;
    os << indent << "ionization potential = " << ip(Z) << endl;

    os << decindent;
  }
  os << decindent;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:

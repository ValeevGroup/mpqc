//
// gaussbas.cc
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

#include <cstdio>
#include <stdexcept>
#include <sstream>
#include <algorithm>
#include <numeric>
#include <cassert>

#include <mpqc_config.h>
#include <util/keyval/keyval.h>
#include <util/misc/formio.h>
#include <util/misc/newstring.h>
#include <util/state/stateio.h>
#include <util/misc/scexception.h>
#include <math/scmat/blocked.h>
#include <chemistry/molecule/molecule.h>
#include <chemistry/qc/basis/gaussshell.h>
#include <chemistry/qc/basis/gaussbas.h>
#include <chemistry/qc/basis/files.h>
#include <chemistry/qc/basis/cartiter.h>
#include <chemistry/qc/basis/transform.h>
#include <chemistry/qc/basis/integral.h>
#include <chemistry/qc/basis/union.h>

using namespace std;
using namespace sc;

//////////////////////////////////////////////////

ClassDesc GaussianBasisSet::Shell::class_desc_(
  typeid(sc::GaussianBasisSet::Shell),"GaussianBasisSet::Shell",1,"public GaussianShell",
  0, 0, 0);

GaussianBasisSet::Shell::Shell(const GaussianBasisSet* basis,
                               unsigned int center,
                               const GaussianShell& shell) :
                                 GaussianShell(shell),
                                 basis_(basis),
                                 center_(center)
{
  MPQC_ASSERT(center < basis->ncenter());
}

bool sc::GaussianBasisSet::Shell::equiv(const Shell& s) const {
  const double* r_this = this->basis()->molecule()->r(this->center());
  const double* r_s = s.basis()->molecule()->r(s.center());
  for(unsigned int xyz=0; xyz<3; ++xyz) {
    if (fabs(r_this[xyz] - r_s[xyz]) > DBL_EPSILON)
      return false;
  }

  return static_cast<const GaussianShell&>(*this).equiv(static_cast<const GaussianShell&>(s));
}

//////////////////////////////////////////////////

static ClassDesc GaussianBasisSet_cd(
  typeid(GaussianBasisSet),"GaussianBasisSet",4,
  "virtual public SavableState",
  0, create<GaussianBasisSet>, create<GaussianBasisSet>);

static bool
skip_atom(bool skip_ghosts, bool include_q,
          const Ref<Molecule> &mol, int iatom)
{
  if (skip_ghosts && mol->charge(iatom) == 0.0) return true;
  // charges do not have basis functions
  if (!include_q && mol->atom_symbol(iatom) == "Q") return true;
  return false;
}

GaussianBasisSet::GaussianBasisSet()
{
}

GaussianBasisSet::GaussianBasisSet(const Ref<KeyVal>&topkeyval)
{
  molecule_ << topkeyval->describedclassvalue("molecule");
  if (molecule_.null()) {
      throw InputError("GaussianBasisSet: no molecule specified",
                       __FILE__, __LINE__, "molecule", 0);
    }

  // see if the user requests pure am or cartesian functions
  int pure;
  pure = topkeyval->booleanvalue("puream");
  if (topkeyval->error() != KeyVal::OK) pure = -1;

  // construct a keyval that contains the basis library
  Ref<KeyVal> keyval;

  if (topkeyval->exists("basisfiles")) {
      Ref<MessageGrp> grp = MessageGrp::get_default_messagegrp();
      Ref<ParsedKeyVal> parsedkv = new ParsedKeyVal();
      char *in_char_array;
      if (grp->me() == 0) {
          ostringstream ostrs;
          // Look at the basisdir and basisfiles variables to find out what
          // basis set files are to be read in.  The files are read on node
          // 0 only.
          ParsedKeyVal::cat_files("basis",topkeyval,ostrs);
          int n = 1 + strlen(ostrs.str().c_str());
          in_char_array = strcpy(new char[n],ostrs.str().c_str());
          grp->bcast(n);
          grp->bcast(in_char_array, n);
        }
      else {
          int n;
          grp->bcast(n);
          in_char_array = new char[n];
          grp->bcast(in_char_array, n);
        }
      parsedkv->parse_string(in_char_array);
      delete[] in_char_array;
      Ref<KeyVal> libkeyval = parsedkv.pointer();
      keyval = new AggregateKeyVal(topkeyval,libkeyval);
    }
  else {
      keyval = topkeyval;
    }

  // Bases keeps track of what basis set data bases have already
  // been read in.  It also handles the conversion of basis
  // names to file names.
  BasisFileSet bases(keyval);
  init(molecule_,keyval,bases,1,pure);
}

Ref<GaussianBasisSet>
GaussianBasisSet::unit()
{
  static Ref<GaussianBasisSet> unit_;

  if (unit_.null()) {
    // make a Molecule composed of 1 classical zero charge at the origin
    Ref<Molecule> molecule = new Molecule;
    molecule->add_atom(molecule->atominfo()->string_to_Z("Q"),
                       0.0, 0.0, 0.0, // xyz
                       "dummy",       // label
                       0.0,           // no mass
                       1, 0.0         // no charge
    );
    std::string name("Unit");
    std::string label = name;
    std::vector<GaussianShell> shells;
    shells.push_back(GaussianShell::unit());
    std::vector<unsigned int> shell_to_center(1,0u);

    unit_ = new GaussianBasisSet(molecule, shells, shell_to_center, name, label);
  }

  return unit_;
}

GaussianBasisSet::GaussianBasisSet(const GaussianBasisSet& gbs) :
  name_(gbs.name_),
  label_(gbs.label_),
  molecule_(gbs.molecule_),
  shells_(gbs.shells_)
{
  init2();
}

GaussianBasisSet::GaussianBasisSet(const Ref<Molecule>& molecule,
                                   const std::vector<GaussianShell>& shell,
                                   const std::vector<unsigned int>& shell_to_center,
                                   std::string name,
                                   std::string label) :
  molecule_(molecule),
  name_(name),
  label_(label)
{
  for(size_t s=0; s<shell.size(); ++s)
    shells_.push_back(Shell(this, shell_to_center[s], shell[s]));

  init2();
}

Ref<GaussianBasisSet>
GaussianBasisSet::operator+(const Ref<GaussianBasisSet>& B)
{
  return new UnionBasisSet(this,B);
}

GaussianBasisSet::GaussianBasisSet(StateIn&s):
  SavableState(s)
{
  if (s.version(::class_desc<GaussianBasisSet>()) < 4) {
      throw FileOperationFailed("cannot restore from old GaussianBasisSet archives",
                                __FILE__, __LINE__, 0,
                                FileOperationFailed::Corrupt,
                                class_desc());
  }

  molecule_ << SavableState::restore_state(s);

  {
    std::vector<GaussianShell> gshells;
    int counter = 0;
    FromStateIn(gshells, s, counter);
    std::vector<unsigned int> shell_to_center;
    FromStateIn(shell_to_center, s, counter);
    MPQC_ASSERT(shell_to_center.size() == gshells.size());

    for(size_t s=0; s<gshells.size(); ++s)
      shells_.push_back(Shell(this, shell_to_center[s], gshells[s]));
  }
  s.get(name_);
  s.get(label_);

  init2();
}

void
GaussianBasisSet::save_data_state(StateOut&s)
{
  SavableState::save_state(molecule_.pointer(),s);
  { // save GaussianShell and center info separately
    std::vector<GaussianShell> gshells;
    std::vector<unsigned int> shell_to_center;
    for(size_t s=0; s<shells_.size(); ++s) {
      gshells.push_back(static_cast<GaussianShell>(shells_[s]));
      shell_to_center.push_back(shells_[s].center());
    }
    int counter = 0;
    ToStateOut(gshells, s, counter);
    ToStateOut(shell_to_center, s, counter);
  }
  s.put(name_);
  s.put(label_);
}

GaussianBasisSet&
GaussianBasisSet::operator=(const GaussianBasisSet& gbs) {
  molecule_ = gbs.molecule_;
  shells_ = gbs.shells_;
  name_ = gbs.name_;
  label_ = gbs.label_;

  init2();
  return *this;
}

void
GaussianBasisSet::init(Ref<Molecule>&molecule,
                       Ref<KeyVal>&keyval,
                       BasisFileSet& bases,
                       int have_userkeyval,
                       int pur)
{
  int pure, havepure;
  pure = pur;
  if (pur == -1) {
      havepure = 0;
    }
  else {
      havepure = 1;
    }

  int skip_ghosts = keyval->booleanvalue("skip_ghosts");
  bool missing_ok = keyval->booleanvalue("missing_ok");
  bool include_q = keyval->booleanvalue("include_q");

  name_ = keyval->stringvalue("name");
  int have_custom, nelement;

  if (keyval->exists("basis")) {
      have_custom = 1;
      nelement = keyval->count("element");
    }
  else {
      have_custom = 0;
      nelement = 0;
      if (name_.empty()) {
          throw InputError("GaussianBasisSet: No name given for basis set\n",
                           __FILE__, __LINE__);
      }
    }

  // Construct label_
  if (!name_.empty())
    label_ = name_;
  else {
    if (have_custom) {
      ostringstream oss;
      Ref<AtomInfo> atominfo = molecule->atominfo();
      // If element is given then construct label_ using element symbol and basis name
      // combinations, e.g. "{ [Fe S1] [Ni S2] [C aug-cc-pVDZ] }"
      if (nelement) {
        oss << "{ ";
        for(int e=0; e<nelement; e++) {
          std::string elementname = keyval->stringvalue("element", e);
          int Z = atominfo->string_to_Z(elementname);
          std::string elemsymbol = atominfo->symbol(Z);
          std::string basisname = keyval->stringvalue("basis", e);
          oss << "[" << elemsymbol << " " << basisname << "] ";
        }
        oss << "}";
      }
      // If element is not given then construct label_ using basis names for each atom
      // e.g. "[ aug-cc-pVDZ cc-pVDZ cc-pVDZ ]"
      else {
        int natom = molecule->natom();
        oss << "[ ";
        for(int a=0; a<natom; a++) {
          std::string basisname = keyval->stringvalue("basis", a);
          oss << basisname << " ";
        }
        oss << "]";
      }
      label_ = oss.str();
    }
  }


  // construct prefixes for each atom: :basis:<atom>:<basisname>:<shell#>
  // and read in the shell
  nbasis_ = 0;
  int ishell = 0;
  const int ncenter = molecule->natom();
  int iatom;
  for (iatom=0; iatom < ncenter; iatom++) {
      if (skip_atom(skip_ghosts,include_q,molecule,iatom)) continue;
      int Z = molecule->Z(iatom);
      // see if there is a specific basisname for this atom
      std::string sbasisname;
      if (have_custom && !nelement) {
          sbasisname = keyval->stringvalue("basis",iatom);
        }
      else if (nelement) {
          int i;
          for (i=0; i<nelement; i++) {
              std::string elementname = keyval->stringvalue("element", i);
              int tmpZ = molecule->atominfo()->string_to_Z(elementname);
              if (tmpZ == Z) {
                  sbasisname = keyval->stringvalue("basis", i);
                  break;
                }
            }
        }
      if (sbasisname.size() == 0) {
          if (name_.empty()) {
            std::ostringstream oss;
            oss << indent << "GaussianBasisSet: no basis name for atom "
                  << iatom
                  << " (Z=" <<molecule->atominfo()->name(Z) << ")"
                  << std::endl;
            throw ProgrammingError(oss.str().c_str(),
                                   __FILE__, __LINE__);
            }
          sbasisname = name_;
        }
      std::string name(molecule->atominfo()->name(Z));
      ishell += count_shells_(keyval, name.c_str(),
                              sbasisname.c_str(), bases, havepure, pure, missing_ok);
    }

  ishell = 0;
  center_to_nshell_.resize(ncenter);
  for (iatom=0; iatom<ncenter; iatom++) {
      if (skip_atom(skip_ghosts,include_q,molecule,iatom)) {
          center_to_nshell_[iatom] = 0;
          continue;
        }
      int Z = molecule->Z(iatom);
      // see if there is a specific basisname for this atom
      std::string sbasisname;
      if (have_custom && !nelement) {
          sbasisname = keyval->stringvalue("basis",iatom);
        }
      else if (nelement) {
          int i;
          for (i=0; i<nelement; i++) {
              std::string elementname = keyval->stringvalue("element", i);
              int tmpZ = molecule->atominfo()->string_to_Z(elementname);
              if (tmpZ == Z) {
                  sbasisname = keyval->stringvalue("basis", i);
                  break;
                }
            }
        }
      if (sbasisname.size() == 0) {
          if (name_.empty()) {
              std::ostringstream oss;
              oss
                  << indent << "GaussianBasisSet: no basis name for atom "
                  << iatom
                  << " (Z=" <<molecule->atominfo()->name(Z) << ")"
                  << std::endl;
              throw ProgrammingError(oss.str().c_str(),
                                     __FILE__, __LINE__);
            }
          sbasisname = name_;
        }

      int ishell_old = ishell;
      std::string name(molecule->atominfo()->name(Z));
      get_shells_(iatom, keyval, name.c_str(),
                  sbasisname.c_str(), bases, havepure, pure, missing_ok);

      center_to_nshell_[iatom] = ishell - ishell_old;
     }

  // delete the name_ if the basis set is customized
  if (have_custom) {
      name_.resize(0);
    }

  // finish with the initialization steps that don't require any
  // external information
  init2(skip_ghosts,include_q);

  // verify that symmetry-equivalent atoms have the same basis set
  // this may help to squash hard-to-detect bugs due to the user providing an erroneous basis keyword
  if (molecule_->point_group()->order() > 1) {
    validate_point_group();
  }
}

void
GaussianBasisSet::validate_point_group() const {
  // loop over all unique centers
  // for each center that has a non-unit orbit, loop over the orbit
  // if this center is not same as the reference center:
  // 1) assert that the number of shells is the same
  // 2) assert that shells on the two centers match in-order
  const int nuniqatoms = molecule_->nunique();
  for(int ua=0; ua<nuniqatoms; ++ua) {
    const int a0 = molecule_->unique(ua);
    const int nshell0 = this->nshell_on_center(a0);
    const int nequivatoms = molecule_->nequivalent(ua);
    if (nequivatoms > 1) {
      for(int ea=0; ea<nequivatoms; ++ea) {
        const int a1 = molecule_->equivalent(ua,ea);
        if (a0 != a1) {
          const int nshell1 = this->nshell_on_center(a1);
          if (nshell0 == nshell1) {
            for(int s=0; s<nshell0; ++s) {
              const int s0 = this->shell_on_center(a0, s);
              const int s1 = this->shell_on_center(a1, s);
              if (static_cast<const GaussianShell&>(this->shell(s0)).equiv(static_cast<const GaussianShell&>(this->shell(s1))) == false) {
                std::ostringstream oss;
                oss << "atoms " << a0 << " and " << a1 << " are symmetry-equivalent but host different basis sets\n"
                    << "the most likely cause is invalid basis vector";
                throw InputError(oss.str().c_str(), __FILE__, __LINE__);
              }
            }
          }
          else {
            std::ostringstream oss;
            oss << "atoms " << a0 << " and " << a1 << " are symmetry-equivalent but host basis sets with different number of shells\n"
                << "the most likely cause is invalid basis vector";
            throw InputError(oss.str().c_str(), __FILE__, __LINE__);
          }
        }
      }
    }
  }
}

double
GaussianBasisSet::r(int icenter, int xyz) const
{
  return molecule_->r(icenter,xyz);
}

void
GaussianBasisSet::init2(int skip_ghosts,bool include_q)
{
  // validate shells_
  // shells on the same atom are consecutive
  if (not shells_.empty()){
    std::vector<bool> center_flags(ncenter(), false); // if true, has found shells on this center

    int current_center = shells_[0].center();
    center_flags[current_center] = true;
    for(size_t s=1; s<shells_.size(); ++s) {
      if (shells_[s].center() != current_center) {
        current_center = shells_[s].center();
        if (center_flags[current_center] == true){ // oops, this shell sits on a center that already has another block of shells
          throw ProgrammingError("GaussianBasisSet: shells must be blocked by center",
                                 __FILE__, __LINE__);
        }
        center_flags[current_center] = true;
      }
    }
  }

  // compute center -> shell map and related
  center_to_shell_.resize(ncenter(), -1);
  center_to_nshell_.resize(ncenter(), 0);
  center_to_nfunction_.resize(ncenter(), 0);
  for(size_t s=0; s<shells_.size(); ++s) {
    const unsigned int center = shells_[s].center();
    if (center_to_shell_[center] == -1)
      center_to_shell_[center] = s;

    center_to_nshell_[center]    += 1;
    center_to_nfunction_[center] += shells_[s].nfunction();
  }

  // compute shell -> function map and related
  shell_to_function_.resize(nshell());
  nbasis_ = 0;
  nprim_ = 0;
  for (size_t s=0; s<nshell(); s++) {
    shell_to_function_[s] = nbasis_;

    nbasis_ += shells_[s].nfunction();
    nprim_ += shells_[s].nprimitive();
  }

  // would like to do this in function_to_shell(), but it is const
  int n = nbasis();
  int nsh = nshell();
  function_to_shell_.resize(n);
  int ifunc = 0;
  for (int i=0; i<nsh; i++) {
      int nfun = operator()(i).nfunction();
      for (int j=0; j<nfun; j++) {
          function_to_shell_[ifunc] = i;
          ifunc++;
        }
    }

  // figure out if any shells are spherical harmonics
  has_pure_ = false;
  for(int i=0; i<nshell(); i++)
    has_pure_ = has_pure_ || shells_[i].has_pure();

  if (matrixkit_.null())
    matrixkit_ = SCMatrixKit::default_matrixkit();
  so_matrixkit_ = new BlockedSCMatrixKit(matrixkit_);
  if (basisdim_.null()) {
    int nb = nshell();
    int *bs = new int[nb];
    for (int s=0; s < nb; s++)
      bs[s] = shell(s).nfunction();
    basisdim_ = new SCDimension(nbasis(), nb, bs, "basis set dimension");
    delete[] bs;
  }
}

void
GaussianBasisSet::set_matrixkit(const Ref<SCMatrixKit>& mk)
{
  matrixkit_ = mk;
  so_matrixkit_ = new BlockedSCMatrixKit(matrixkit_);
}


int
GaussianBasisSet::count_shells_(Ref<KeyVal>& keyval, const char* element, const char* basisname, BasisFileSet& bases,
				int havepure, int pure, bool missing_ok)
{
  int nshell = 0;
  char keyword[KeyVal::MaxKeywordLength];

  sprintf(keyword,":basis:%s:%s",element,basisname);
  bool exists = keyval->exists(keyword);
  if (!exists) {
    keyval = bases.keyval(keyval, basisname);
    exists = keyval->exists(keyword);
    if (!exists) {
      if (missing_ok) return 0;
      std::ostringstream oss;
      oss << scprintf("GaussianBasisSet::count_shells_ couldn't find \"%s\":\n", keyword);
      ExEnv::err0() << indent << oss.str();
      keyval->errortrace(ExEnv::err0());
      throw InputError(oss.str().c_str(), __FILE__, __LINE__);
    }
  }

  // Check if the basis set is an array of shells
  keyval->count(keyword);
  if (keyval->error() != KeyVal::OK) {
    nshell = count_even_temp_shells_(keyval, element, basisname, havepure, pure);
  }
  else {
    recursively_get_shell(nshell, keyval, element, basisname,
                          bases, havepure, pure, 0, false);
  }

  return nshell;
}

void
GaussianBasisSet::get_shells_(unsigned int atom, Ref<KeyVal>& keyval, const char* element, const char* basisname, BasisFileSet& bases,
			      int havepure, int pure, bool missing_ok)
{
  char keyword[KeyVal::MaxKeywordLength];

  sprintf(keyword,":basis:%s:%s:am",
          element,basisname);
  if (keyval->exists(keyword)) {
    get_even_temp_shells_(atom, keyval, element, basisname, havepure, pure);
  }
  else {
    recursively_get_shell(atom, keyval, element, basisname,
                          bases, havepure, pure, 1, missing_ok);
  }
}


int
GaussianBasisSet::count_even_temp_shells_(Ref<KeyVal>& keyval, const char* element, const char* basisname,
                                          int havepure, int pure)
{
  int nshell = 0;
  char keyword[KeyVal::MaxKeywordLength];

  sprintf(keyword,":basis:%s:%s:am",element,basisname);
  if (!keyval->exists(keyword)) {
    sprintf(keyword,":basis:%s:%s",element,basisname);
    ExEnv::err0() << indent
                  << scprintf("GaussianBasisSet::count_even_temp_shells_ -- couldn't read \"%s\":\n", keyword);
    throw std::runtime_error("GaussianBasisSet::count_even_temp_shells_ -- basis set specification is invalid");
  }

  // count the number of even-tempered primitive blocks
  int nblocks = keyval->count(keyword) - 1;
  if (keyval->error() != KeyVal::OK) {
    ExEnv::err0() << indent
                  << scprintf("GaussianBasisSet::count_even_temp_shells_ -- couldn't read \"%s\":\n", keyword);
    throw std::runtime_error("GaussianBasisSet::count_even_temp_shells_ -- failed to read am");
  }
  if (nblocks == -1)
    return 0;

  sprintf(keyword,":basis:%s:%s:nprim", element, basisname);
  int j = keyval->count(keyword) - 1;
  if (nblocks != j) {
    ExEnv::err0() << indent
                  << scprintf("GaussianBasisSet::count_even_temp_shells_ -- problem reading \"%s\":\n", keyword);
    throw std::runtime_error("GaussianBasisSet::count_even_temp_shells_ -- am and nprim have different dimensions");
  }

  for(int b=0; b<=nblocks; b++) {
    sprintf(keyword,":basis:%s:%s:nprim:%d", element, basisname, b);
    int nprim = keyval->intvalue(keyword);
    if (nprim <= 0) {
      ExEnv::err0() << indent
                    << scprintf("GaussianBasisSet::count_even_temp_shells_ -- problem with \"%s\":\n", keyword);
      throw std::runtime_error("GaussianBasisSet::count_shells_ -- the number of primitives has to be positive");
    }
    nshell += nprim;
  }

  return nshell;
}


void
GaussianBasisSet::get_even_temp_shells_(unsigned int atom, Ref<KeyVal>& keyval, const char* element, const char* basisname,
                                        int havepure, int pure)
{
  char keyword[KeyVal::MaxKeywordLength];

  // count the number of even-tempered primitive blocks
  sprintf(keyword,":basis:%s:%s:am",
          element,basisname);
  int nblocks = keyval->count(keyword) - 1;
  if (keyval->error() != KeyVal::OK) {
    ExEnv::err0() << indent
                  << scprintf("GaussianBasisSet::get_even_temp_shells_ -- couldn't read \"%s\":\n", keyword);
    throw std::runtime_error("GaussianBasisSet::get_even_temp_shells_ -- failed to read am");
  }
  if (nblocks == -1)
    return;

  sprintf(keyword,":basis:%s:%s:nprim", element, basisname);
  int j = keyval->count(keyword) - 1;
  if (nblocks != j) {
    ExEnv::err0() << indent
                  << scprintf("GaussianBasisSet::get_even_temp_shells_ -- problem reading \"%s\":\n", keyword);
    throw std::runtime_error("GaussianBasisSet::get_even_temp_shells_ -- am and nprim have different dimensions");
  }

  sprintf(keyword,":basis:%s:%s:last_exp", element, basisname);
  bool have_last_exp = keyval->exists(keyword);

  sprintf(keyword,":basis:%s:%s:first_exp", element, basisname);
  bool have_first_exp = keyval->exists(keyword);

  sprintf(keyword,":basis:%s:%s:exp_ratio", element, basisname);
  bool have_exp_ratio = keyval->exists(keyword);

  if ( !have_first_exp && !have_last_exp )
    throw std::runtime_error("GaussianBasisSet::get_even_temp_shells_ -- neither last_exp nor first_exp has been specified");

  if ( have_first_exp && have_last_exp && have_exp_ratio)
    throw std::runtime_error("GaussianBasisSet::get_even_temp_shells_ -- only two of (last_exp,first_exp,exp_ratio) can be specified");

  if ( !have_first_exp && !have_last_exp && have_exp_ratio)
    throw std::runtime_error("GaussianBasisSet::get_even_temp_shells_ -- any two of (last_exp,first_exp,exp_ratio) must be specified");

  for(int b=0; b<=nblocks; b++) {

    sprintf(keyword,":basis:%s:%s:nprim:%d", element, basisname, b);
    int nprim = keyval->intvalue(keyword);
    if (nprim <= 0) {
      ExEnv::err0() << indent
                    << scprintf("GaussianBasisSet::get_even_temp_shells_ -- problem with \"%s\":\n", keyword);
      throw std::runtime_error("GaussianBasisSet::get_even_temp_shells_ -- the number of primitives has to be positive");
    }

    sprintf(keyword,":basis:%s:%s:am:%d", element, basisname, b);
    int l = keyval->intvalue(keyword);
    if (l < 0) {
      ExEnv::err0() << indent
                    << scprintf("GaussianBasisSet::get_even_temp_shells_ -- problem with \"%s\":\n", keyword);
      throw std::runtime_error("GaussianBasisSet::get_even_temp_shells_ -- angular momentum has to be non-negative");
    }

    double alpha0, alphaN, beta;
    if (have_first_exp) {
      sprintf(keyword,":basis:%s:%s:first_exp:%d", element, basisname, b);
      alpha0 = keyval->doublevalue(keyword);
      if (alpha0 <= 0.0) {
        ExEnv::err0() << indent
                      << scprintf("GaussianBasisSet::get_even_temp_shells_ -- problem with \"%s\":\n", keyword);
        throw std::runtime_error("GaussianBasisSet::get_even_temp_shells_ -- orbital exponents have to be positive");
      }
    }

    if (have_last_exp) {
      sprintf(keyword,":basis:%s:%s:last_exp:%d", element, basisname, b);
      alphaN = keyval->doublevalue(keyword);
      if (alphaN <= 0.0) {
        ExEnv::err0() << indent
                      << scprintf("GaussianBasisSet::get_even_temp_shells_ -- problem with \"%s\":\n", keyword);
        throw std::runtime_error("GaussianBasisSet::get_even_temp_shells_ -- orbital exponents have to be positive");
      }
    }

    if (have_last_exp && have_first_exp) {
      if (alphaN > alpha0) {
        ExEnv::err0() << indent
                      << scprintf("GaussianBasisSet::get_even_temp_shells_ -- problem with \"%s\":\n", keyword);
        throw std::runtime_error("GaussianBasisSet::get_even_temp_shells_ -- last_exps[i] must be smaller than first_exp[i]");
      }
      if (nprim > 1)
        beta = pow(alpha0/alphaN,1.0/(nprim-1));
      else
        beta = 1.0;
    }
    else {
      sprintf(keyword,":basis:%s:%s:exp_ratio:%d", element, basisname, b);
      beta = keyval->doublevalue(keyword);
      if (beta <= 1.0) {
        ExEnv::err0() << indent
                      << scprintf("GaussianBasisSet::get_even_temp_shells_ -- problem with \"%s\":\n", keyword);
        throw std::runtime_error("GaussianBasisSet::get_even_temp_shells_ -- exponent ratio has to be greater than 1.0");
      }
      if (have_last_exp)
        alpha0 = alphaN * pow(beta,nprim-1);
    }

    double alpha = alpha0;
    for(int p=0; p<nprim; p++, alpha /= beta ) {
      std::vector<unsigned int> am(1, l);
      std::vector<double> exps(1, alpha);
      std::vector<double> coeffs(1, 1.0);

      std::vector<bool> puream(1, (l<1) ? false : (havepure ? true : false) );
      shells_.push_back( Shell( this, atom, GaussianShell(am, puream, exps, coeffs)) );

      // Recompute beta at each step to maintain accuracy
      if (have_last_exp && have_first_exp) {
        int nprim_left = nprim - p;
        if (nprim_left > 1)
          beta = pow(alpha/alphaN,1.0/(nprim_left-1));
        else
          beta = 1.0;
      }
    }
  }
}

void
GaussianBasisSet::
  recursively_get_shell(unsigned int atom,Ref<KeyVal>&keyval,
                        const char*element,
                        const char*basisname,
                        BasisFileSet &bases,
                        int havepure,int pure,
                        int get, bool missing_ok)
{
  char keyword[KeyVal::MaxKeywordLength],prefix[KeyVal::MaxKeywordLength];

  sprintf(keyword,":basis:%s:%s",
	  element,basisname);
  int count = keyval->count(keyword);
  if (keyval->error() != KeyVal::OK) {
      keyval = bases.keyval(keyval, basisname);
    }
  count = keyval->count(keyword);
  if (keyval->error() != KeyVal::OK) {
      if (missing_ok) return;
      std::ostringstream oss;
      oss << scprintf("GaussianBasisSet::count_shells_ couldn't find \"%s\":\n", keyword);
      ExEnv::err0() << indent << oss.str();
      keyval->errortrace(ExEnv::err0());
      throw InputError(oss.str().c_str(), __FILE__, __LINE__);
    }
  if (!count) return;
  for (int j=0; j<count; j++) {
      sprintf(prefix,":basis:%s:%s",
	      element,basisname);
      Ref<KeyVal> prefixkeyval = new PrefixKeyVal(keyval,prefix,j);
      if (prefixkeyval->exists("get")) {
          std::string newbasis = prefixkeyval->stringvalue("get");
          if (newbasis.empty()) {
	      ExEnv::err0() << indent << "GaussianBasisSet: "
                   << scprintf("error processing get for \"%s\"\n", prefix);
              keyval->errortrace(ExEnv::err0());
	      exit(1);
	    }
	  recursively_get_shell(atom,keyval,element,newbasis.c_str(),bases,
                                havepure,pure,get,missing_ok);
	}
      else {
          if (get) {
	      shells_.push_back(Shell(this, atom, GaussianShell(prefixkeyval, havepure ? pure : -1)) );
	    }
	}
    }
}

GaussianBasisSet::~GaussianBasisSet()
{
}

unsigned int
GaussianBasisSet::max_nfunction_in_shell() const
{
  unsigned int max = 0;
  for (size_t i=0; i<nshell(); i++) {
      if (max < shells_[i].nfunction()) max = shells_[i].nfunction();
    }
  return max;
}

unsigned int
GaussianBasisSet::max_ncontraction() const
{
  unsigned int max = 0;
  for (size_t i=0; i<nshell(); i++) {
      if (max < shells_[i].ncontraction()) max = shells_[i].ncontraction();
    }
  return max;
}

unsigned int
GaussianBasisSet::max_angular_momentum() const
{
  unsigned int max = 0;
  for (size_t i=0; i<nshell(); i++) {
      unsigned int maxshi = shells_[i].max_angular_momentum();
      if (max < maxshi) max = maxshi;
    }
  return max;
}

unsigned int
GaussianBasisSet::max_cartesian() const
{
  unsigned int max = 0;
  for (size_t i=0; i<nshell(); i++) {
      unsigned int maxshi = shells_[i].max_cartesian();
      if (max < maxshi) max = maxshi;
    }
  return max;
}

unsigned int
GaussianBasisSet::max_ncartesian_in_shell(int aminc) const
{
  unsigned int max = 0;
  for (size_t i=0; i<nshell(); i++) {
      unsigned int maxshi = shells_[i].ncartesian_with_aminc(aminc);
      if (max < maxshi) max = maxshi;
    }
  return max;
}

unsigned int
GaussianBasisSet::max_nprimitive_in_shell() const
{
  unsigned int max = 0;
  for (size_t i=0; i<nshell(); i++) {
      if (max < shells_[i].nprimitive()) max = shells_[i].nprimitive();
    }
  return max;
}

unsigned int
GaussianBasisSet::max_am_for_contraction(int con) const
{
  unsigned int max = 0;
  for (size_t i=0; i<nshell(); i++) {
      if (shells_[i].ncontraction() <= con) continue;
      unsigned int maxshi = shells_[i].am(con);
      if (max < maxshi) max = maxshi;
    }
  return max;
}

int
GaussianBasisSet::function_to_shell(int func) const
{
  return function_to_shell_[func];
}

unsigned int
GaussianBasisSet::ncenter() const
{
  return molecule_->natom();
}

int
GaussianBasisSet::nshell_on_center(int icenter) const
{
  return center_to_nshell_[icenter];
}

unsigned int
GaussianBasisSet::nbasis_on_center(int icenter) const
{
  return center_to_nfunction_[icenter];
}

int
GaussianBasisSet::shell_on_center(int icenter, int ishell) const
{
  return center_to_shell_[icenter] + ishell;
}

const GaussianBasisSet::Shell&
GaussianBasisSet::operator()(int icenter,int ishell) const
{
  return shells_[center_to_shell_[icenter] + ishell];
}

GaussianBasisSet::Shell&
GaussianBasisSet::operator()(int icenter,int ishell)
{
  return shells_[center_to_shell_[icenter] + ishell];
}

#ifdef MPQC_NEW_FEATURES
mpqc::range GaussianBasisSet::range(int s) const {
    int f = this->shell_to_function(s);
    return mpqc::range(f, f + this->shell(s).nfunction());
}
#endif

int
GaussianBasisSet::find(int C, const GaussianShell& S) const {
  const int sfirst = this->center_to_shell_[C];
  const int sfence = sfirst + this->center_to_nshell_[C];
  for(int s=sfirst; s<sfence; ++s)
    if ( static_cast<const GaussianShell&>(shells_[s]).equiv(S) )
      return s;
  return -1;
}

int
GaussianBasisSet::equiv(const Ref<GaussianBasisSet> &b)
{
  if (nshell() != b->nshell()) return 0;
  for (int i=0; i<nshell(); i++) {
      if (not shells_[i].equiv(b->shells_[i]))
        return 0;
    }
  return 1;
}

void
GaussianBasisSet::print_brief(ostream& os) const
{
  os << indent
     << "GaussianBasisSet:" << endl << incindent
     << indent << "nbasis = " << nbasis_ << endl
     << indent << "nshell = " << nshell() << endl
     << indent << "nprim  = " << nprim_ << endl;
  if (!name_.empty()) {
      os << indent
         << "name = \"" << name_ << "\"" << endl;
  }
  if (name_ != label_) {
    os << indent
       << "label = \"" << label_ << "\"" << endl;
  }
  os << decindent;
}

void
GaussianBasisSet::print(ostream& os) const
{
  print_brief(os);
  if (!SCFormIO::getverbose(os)) return;

  os << incindent;

  // Loop over centers
  int icenter = 0;
  int ioshell = 0;
  for (icenter=0; icenter < ncenter(); icenter++) {
      os << endl << indent
         << scprintf(
             "center %d: %12.8f %12.8f %12.8f, nshell = %d, shellnum = %d\n",
             icenter,
             r(icenter,0),
             r(icenter,1),
             r(icenter,2),
             center_to_nshell_[icenter],
             center_to_shell_[icenter]);
      for (int ishell=0; ishell < center_to_nshell_[icenter]; ishell++) {
	  os << indent
             << scprintf("Shell %d: functionnum = %d\n",
                         ishell,shell_to_function_[ioshell]);
          os << incindent;
	  operator()(icenter,ishell).print(os);
          os << decindent;
          ioshell++;
	}
    }

  os << decindent;
}

Ref<GaussianBasisSet>
sc::operator+(const Ref<GaussianBasisSet>& A, const Ref<GaussianBasisSet>& B)
{
  GaussianBasisSet& aref = *(A.pointer());
  return aref + B;
}

void
GaussianBasisSet::init(std::string name,
                       std::string label,
                       const Ref<Molecule> &molecule,
                       const std::vector<Shell>& shell)

{
  name_ = name;
  label_ = label;
  molecule_ = molecule;
  shells_ = shell;

  init2();
}

/////////////////////////////////////////////////////////////////////////////
// GaussianBasisSet::ValueData

GaussianBasisSet::ValueData::ValueData(
    const Ref<GaussianBasisSet> &basis,
    const Ref<Integral> &integral)
{
  maxam_ = basis->max_angular_momentum();

  civec_ = new CartesianIter *[maxam_+1];
  sivec_ = new SphericalTransformIter *[maxam_+1];
  for (int i=0; i<=maxam_; i++) {
      civec_[i] = integral->new_cartesian_iter(i);
      // spherical transforms are necessary even for p-type shells
      if (i>0) sivec_[i] = integral->new_spherical_transform_iter(i);
      else sivec_[i] = 0;
    }
}

GaussianBasisSet::ValueData::~ValueData()
{
  for (int i=0; i<=maxam_; i++) {
      delete civec_[i];
      delete sivec_[i];
    }
  delete[] civec_;
  delete[] sivec_;
}

int
sc::ishell_on_center(int icenter, const Ref<GaussianBasisSet>& bs,
		             const GaussianShell& A)
{
  for (int jshell=0; jshell < bs->nshell_on_center(icenter); jshell++) {
    const GaussianShell& B = bs->shell(icenter, jshell);
    if (A.equiv(B))
      return bs->shell_on_center(icenter, jshell);
  }
  throw ProgrammingError("ishell_on_center() -- did not find the given shell",__FILE__,__LINE__);
}

////

namespace sc {

  std::vector<unsigned int> operator<<(const GaussianBasisSet& B, const GaussianBasisSet& A) {

    typedef std::vector<int> BFMap;
    const bool must_use_same_basis = true;
    BFMap int_map = map(B, A);
    BFMap::iterator i = std::find(int_map.begin(), int_map.end(), -1);
    if (i != int_map.end())
      throw std::logic_error("Cannot construct basis function map");

    std::vector<unsigned int> result(int_map.size());
    std::copy(int_map.begin(), int_map.end(), result.begin());
    return result;
  }

  std::vector<int>
  map(const GaussianBasisSet& B, const GaussianBasisSet& A) {

    // validate input
    // TODO implement Molecule::operator==()
    const bool valid_input = ( A.molecule().pointer() == B.molecule().pointer() // no Molecule::operator== yet, hence compare pointers
                             );
    if (!valid_input)
      throw std::logic_error("Cannot construct basis function map");

    std::vector<int> map(A.nbasis());
    // loop over shells in A
    for(int sa=0; sa<A.nshell(); ++sa) {
      const int center = A.shell_to_center(sa);
      const GaussianShell& a = A[sa];
      bool found_equivalent = false;
      const int fa_start = A.shell_to_function(sa);
      const int fa_fence = fa_start + a.nfunction();

      // loop over all shells in A on the same center
      const int sb_start = B.shell_on_center(center,0);
      const int sb_fence = sb_start + B.nshell_on_center(center);
      for(int sb=sb_start; sb<sb_fence; ++sb) {
        const GaussianShell& b = B[sb];

        if (a.equiv(b)) {
          if (found_equivalent) throw std::logic_error("found more than 1 match for a shell in operator<<(GaussianBasisSet,GaussianBasisSet)");
          found_equivalent = true;

          // map each basis function
          const int fb_start = B.shell_to_function(sb);
          const int fb_fence = fb_start + b.nfunction();
          for(int fa=fa_start, fb=fb_start; fa<fa_fence; ++fa, ++fb)
            map[fa] = fb;

        }

      }

      if (!found_equivalent) {
        for(int fa=fa_start; fa<fa_fence; ++fa)
          map[fa] = -1;
      }

    }

    return map;
  }

}

/////////////////////////////////////////////////////////////////////////////

GaussianBasisSetMap::GaussianBasisSetMap(const Ref<GaussianBasisSet>& source,
                                         const Ref<GaussianBasisSet>& target) :
                                         source_(source),
                                         target_(target)
{
  if (source->molecule().pointer() != target->molecule().pointer())
    throw std::logic_error("GaussianBasisSetMap -- cannot map basis sets, molecules are different");

  Ref<SCMatrixKit> matrixkit = source->matrixkit();
  Ref<Molecule> molecule = source->molecule();
  const int ncenter = source->ncenter();

  // construct shell map
  const int nshell_source = source->nshell();
  smap_.resize(nshell_source);
  for(int s=0; s<nshell_source; ++s) {
    const int c = source->shell_to_center(s);
    const int s_in_target = target->find(c, source->shell(s));
    if (s_in_target != -1) // s2 is not found in bs1
      smap_[s] = s_in_target;
    else
      throw std::logic_error("GaussianBasisSetMap: cannot construct map");
  }

  // construct function map
  const int nbasis_source = source->nbasis();
  fmap_.resize(nbasis_source);
  for(int s=0; s<nshell_source; ++s) {
    const int nf = source->shell(s).nfunction();
    const int f_source_offset = source->shell_to_function(s);
    const int f_target_offset = target->shell_to_function( smap_[s] );
    for(int f=0; f<nf; ++f) {
      fmap_[f + f_source_offset] = f + f_target_offset;
    }
  }

  {
    fblock_to_function_.push_back(0);
    int current_function_target = fmap_[0];
    int nfunction_in_curr_fblock = 1;
    for(int f=1; f<nbasis_source; ++f) {
      const int next_function_target = fmap_[f];
      if (next_function_target == current_function_target + 1) {
        // still in the current fblock
        ++nfunction_in_curr_fblock;
      }
      else {
        // finalize this fblock
        fblock_size_.push_back(nfunction_in_curr_fblock);
        nfunction_in_curr_fblock = 1;
        // move onto the new fblock
        fblock_to_function_.push_back(f);
      }
      current_function_target = next_function_target;
    }
    // finalize the last fblock
    fblock_size_.push_back(nfunction_in_curr_fblock);
  }
}

GaussianBasisSetMap::~GaussianBasisSetMap()
{
}

int
GaussianBasisSetMap::map_shell(int s12) const { return smap_[s12]; }

int
GaussianBasisSetMap::map_function(int f12) const { return fmap_[f12]; }

int
GaussianBasisSetMap::nfblock() const { return fblock_size_.size(); }

int
GaussianBasisSetMap::fblock_to_function(int b) const { return fblock_to_function_[b]; }

int
GaussianBasisSetMap::fblock_size(int b) const { return fblock_size_[b]; }

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:

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

#ifdef __GNUC__
#pragma implementation
#endif

#include <stdio.h>
#include <stdexcept>

#include <scconfig.h>
#ifdef HAVE_SSTREAM
#  include <sstream>
#else
#  include <strstream.h>
#endif

#include <util/keyval/keyval.h>
#include <util/class/scexception.h>
#include <util/misc/formio.h>
#include <util/misc/newstring.h>
#include <util/state/stateio.h>
#include <math/scmat/blocked.h>
#include <chemistry/molecule/molecule.h>
#include <chemistry/qc/basis/gaussshell.h>
#include <chemistry/qc/basis/gaussbas.h>
#include <chemistry/qc/basis/files.h>
#include <chemistry/qc/basis/cartiter.h>
#include <chemistry/qc/basis/transform.h>
#include <chemistry/qc/basis/integral.h>

using namespace std;
using namespace sc;

static ClassDesc GaussianBasisSet_cd(
  typeid(GaussianBasisSet),"GaussianBasisSet",3,
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
      ExEnv::err0() << indent << "GaussianBasisSet: no \"molecule\"\n";
      abort();
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
#ifdef HAVE_SSTREAM
          ostringstream ostrs;
#else
          ostrstream ostrs;
#endif
          // Look at the basisdir and basisfiles variables to find out what
          // basis set files are to be read in.  The files are read on node
          // 0 only.
          ParsedKeyVal::cat_files("basis",topkeyval,ostrs);
#ifdef HAVE_SSTREAM
          int n = 1 + strlen(ostrs.str().c_str());
          in_char_array = strcpy(new char[n],ostrs.str().c_str());
#else
          ostrs << ends;
          in_char_array = ostrs.str();
          int n = ostrs.pcount();
#endif
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

  // if there isn't a matrixkit in the input, init2() will assign the
  // default matrixkit
  matrixkit_ << keyval->describedclassvalue("matrixkit");
  
  // Bases keeps track of what basis set data bases have already
  // been read in.  It also handles the conversion of basis
  // names to file names.
  BasisFileSet bases(keyval);
  init(molecule_,keyval,bases,1,pure);
}

GaussianBasisSet::GaussianBasisSet(UnitType)
{
  molecule_ = new Molecule;
  molecule_->add_atom(molecule()->atominfo()->string_to_Z("Q"),
                      0.0, 0.0, 0.0, // xyz
                      "dummy",       // label
                      0.0,           // no mass
                      1, 0.0         // no charge
      );
  name_ = new_string("Unit");
  label_ = new_string(name_);
  shell_ = new GaussianShell*[1];

  double *exp = new double[1];
  int *am = new int[1];
  int *pure = new int[1];
  double **c = new double*[1];
  *c = new double[1];
  exp[0] = 0.0;
  am[0] = 0;
  pure[0] = 0;
  c[0][0] = 1.0;
  shell_[0] = new GaussianShell(1,1,exp,am,pure,c,
                                GaussianShell::Unnormalized,
                                false);

  ncenter_ = 1;
  nshell_ = 1;
  center_to_nshell_.push_back(1);
  init2(0,1);
}

GaussianBasisSet::GaussianBasisSet(const GaussianBasisSet& gbs) :
  molecule_(gbs.molecule_),
  matrixkit_(gbs.matrixkit_),
  basisdim_(gbs.basisdim_),
  ncenter_(gbs.ncenter_),
  nshell_(gbs.nshell_)
{
  int i,j;
  
  name_ = new_string(gbs.name_);
  label_ = new_string(gbs.label_);

  center_to_nshell_.resize(ncenter_);
  for (i=0; i < ncenter_; i++) {
      center_to_nshell_[i] = gbs.center_to_nshell_[i];
    }
  
  shell_ = new GaussianShell*[nshell_];
  for (i=0; i<nshell_; i++) {
      const GaussianShell& gsi = gbs(i);

      int nc=gsi.ncontraction();
      int np=gsi.nprimitive();
      
      int *ams = new int[nc];
      int *pure = new int[nc];
      double *exps = new double[np];
      double **coefs = new double*[nc];

      for (j=0; j < nc; j++) {
          ams[j] = gsi.am(j);
          pure[j] = gsi.is_pure(j);
          coefs[j] = new double[np];
          for (int k=0; k < np; k++)
              coefs[j][k] = gsi.coefficient_unnorm(j,k);
        }

      for (j=0; j < np; j++)
          exps[j] = gsi.exponent(j);
      
      shell_[i] = new GaussianShell(nc, np, exps, ams, pure, coefs,
                                   GaussianShell::Unnormalized);
    }

  init2();
}

GaussianBasisSet::GaussianBasisSet(const char* name,
                                   const char* label,
                                   const Ref<Molecule>& molecule,
                                   const Ref<SCMatrixKit>& matrixkit,
                                   const RefSCDimension& basisdim,
                                   const int ncenter,
                                   const int nshell,
                                   GaussianShell** shell,
                                   const std::vector<int>& center_to_nshell) :
  molecule_(molecule),
  matrixkit_(matrixkit),
  basisdim_(basisdim),
  ncenter_(ncenter),
  nshell_(nshell),
  shell_(shell),
  center_to_nshell_(center_to_nshell)
{
  name_ = new_string(name);
  label_ = new_string(label);
  
  init2();
}

Ref<GaussianBasisSet>
GaussianBasisSet::operator+(const Ref<GaussianBasisSet>& B)
{
  return (GaussianBasisSetSum(this,B)).bs12();
}

GaussianBasisSet::GaussianBasisSet(StateIn&s):
  SavableState(s)
{
  matrixkit_ = SCMatrixKit::default_matrixkit();

  if (s.version(::class_desc<GaussianBasisSet>()) < 3) {
      // read the duplicate length saved in older versions
      int junk;
      s.get(junk);
    }
  s.get(center_to_nshell_);

  molecule_ << SavableState::restore_state(s);
  basisdim_ << SavableState::restore_state(s);


  ncenter_ = center_to_nshell_.size();
  s.getstring(name_);
  s.getstring(label_);

  nshell_ = 0;
  int i;
  for (i=0; i<ncenter_; i++) {
      nshell_ += center_to_nshell_[i];
    }
  
  shell_ = new GaussianShell*[nshell_];
  for (i=0; i<nshell_; i++) {
      shell_[i] = new GaussianShell(s);
    }

  init2();
}

void
GaussianBasisSet::save_data_state(StateOut&s)
{
  s.put(center_to_nshell_);

  SavableState::save_state(molecule_.pointer(),s);
  SavableState::save_state(basisdim_.pointer(),s);
  
  s.putstring(name_);
  s.putstring(label_);
  for (int i=0; i<nshell_; i++) {
      shell_[i]->save_object_state(s);
    }
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

  name_ = keyval->pcharvalue("name");
  int have_custom, nelement;

  if (keyval->exists("basis")) {
      have_custom = 1;
      nelement = keyval->count("element");
    }
  else {
      have_custom = 0;
      nelement = 0;
      if (!name_) {
          ExEnv::err0() << indent
               << "GaussianBasisSet: No name given for basis set\n";
          abort();
        }
    }

  // Construct label_
  if (name_)
    label_ = new_string(name_);
  else {
    if (have_custom) {
      ostringstream oss;
      Ref<AtomInfo> atominfo = molecule->atominfo();
      // If element is given then construct label_ using element symbol and basis name
      // combinations, e.g. "{ [Fe S1] [Ni S2] [C aug-cc-pVDZ] }"
      if (nelement) {
        oss << "{ ";
        for(int e=0; e<nelement; e++) {
          char* tmpelementname = keyval->pcharvalue("element", e);
          int Z = atominfo->string_to_Z(tmpelementname);
          std::string elemsymbol = atominfo->symbol(Z);
          char* basisname = keyval->pcharvalue("basis", e);
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
          char* basisname = keyval->pcharvalue("basis", a);
          oss << basisname << " ";
        }
        oss << "]";
      }
      std::string label = oss.str();
      label_ = new char[label.size() + 1];
      strcpy(label_,label.c_str());
    }
  }


  // construct prefixes for each atom: :basis:<atom>:<basisname>:<shell#>
  // and read in the shell
  nbasis_ = 0;
  int ishell = 0;
  ncenter_ = molecule->natom();
  int iatom;
  for (iatom=0; iatom < ncenter_; iatom++) {
      if (skip_atom(skip_ghosts,include_q,molecule,iatom)) continue;
      int Z = molecule->Z(iatom);
      // see if there is a specific basisname for this atom
      char* sbasisname = 0;
      if (have_custom && !nelement) {
          sbasisname = keyval->pcharvalue("basis",iatom);
        }
      else if (nelement) {
          int i;
          for (i=0; i<nelement; i++) {
              char *tmpelementname = keyval->pcharvalue("element", i);
              int tmpZ = molecule->atominfo()->string_to_Z(tmpelementname);
              if (tmpZ == Z) {
                  sbasisname = keyval->pcharvalue("basis", i);
                  break;
                }
              delete[] tmpelementname;
            }
        }
      if (!sbasisname) {
          if (!name_) {
              ExEnv::err0()
                  << indent << "GaussianBasisSet: no basis name for atom "
                  << iatom
                  << " (Z=" <<molecule->atominfo()->name(Z) << ")"
                  << std::endl;
              abort();
            }
          sbasisname = strcpy(new char[strlen(name_)+1],name_);
        }
      std::string name(molecule->atominfo()->name(Z));
      ishell += count_shells_(keyval, name.c_str(),
                              sbasisname, bases, havepure, pure, missing_ok);
      delete[] sbasisname;
    }
  nshell_ = ishell;
  shell_ = new GaussianShell*[nshell_];
  ishell = 0;
  center_to_nshell_.resize(ncenter_);
  for (iatom=0; iatom<ncenter_; iatom++) {
      if (skip_atom(skip_ghosts,include_q,molecule,iatom)) {
          center_to_nshell_[iatom] = 0;
          continue;
        }
      int Z = molecule->Z(iatom);
      // see if there is a specific basisname for this atom
      char* sbasisname = 0;
      if (have_custom && !nelement) {
          sbasisname = keyval->pcharvalue("basis",iatom);
        }
      else if (nelement) {
          int i;
          for (i=0; i<nelement; i++) {
              char *tmpelementname = keyval->pcharvalue("element", i);
              int tmpZ = molecule->atominfo()->string_to_Z(tmpelementname);
              if (tmpZ == Z) {
                  sbasisname = keyval->pcharvalue("basis", i);
                  break;
                }
              delete[] tmpelementname;
            }
        }
      if (!sbasisname) {
          if (!name_) {
              ExEnv::err0()
                  << indent << "GaussianBasisSet: no basis name for atom "
                  << iatom
                  << " (Z=" <<molecule->atominfo()->name(Z) << ")"
                  << std::endl;
              abort();
            }
          sbasisname = strcpy(new char[strlen(name_)+1],name_);
        }

      int ishell_old = ishell;
      std::string name(molecule->atominfo()->name(Z));
      get_shells_(ishell, keyval, name.c_str(),
                  sbasisname, bases, havepure, pure, missing_ok);

      center_to_nshell_[iatom] = ishell - ishell_old;

      delete[] sbasisname;
     }

  // delete the name_ if the basis set is customized
  if (have_custom) {
      delete[] name_;
      name_ = 0;
    }

  // finish with the initialization steps that don't require any
  // external information
  init2(skip_ghosts,include_q);
}

double
GaussianBasisSet::r(int icenter, int xyz) const
{
  return molecule_->r(icenter,xyz);
}

void
GaussianBasisSet::init2(int skip_ghosts,bool include_q)
{
  // center_to_shell_ and shell_to_center_
  shell_to_center_.resize(nshell_);
  center_to_shell_.resize(ncenter_);
  center_to_nbasis_.resize(ncenter_);
  int ishell = 0;
  for (int icenter=0; icenter<ncenter_; icenter++) {
      if (skip_atom(skip_ghosts,include_q,molecule(),icenter)) {
          center_to_shell_[icenter] = -1;
          center_to_nbasis_[icenter] = 0;
          continue;
        }
      int j;
      center_to_shell_[icenter] = ishell;
      center_to_nbasis_[icenter] = 0;
      for (j = 0; j<center_to_nshell_[icenter]; j++) {
          center_to_nbasis_[icenter] += shell_[ishell+j]->nfunction();
        }
      ishell += center_to_nshell_[icenter];
      for (j = center_to_shell_[icenter]; j<ishell; j++) {
	  shell_to_center_[j] = icenter;
	}
     }

  // compute nbasis_ and shell_to_function_[]
  shell_to_function_.resize(nshell_);
  shell_to_primitive_.resize(nshell_);
  nbasis_ = 0;
  nprim_ = 0;
  for (ishell=0; ishell<nshell_; ishell++) {
      shell_to_function_[ishell] = nbasis_;
      shell_to_primitive_[ishell] = nprim_;
      nbasis_ += shell_[ishell]->nfunction();
      nprim_ += shell_[ishell]->nprimitive();
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
  for(int i=0; i<nsh; i++)
    has_pure_ = has_pure_ || shell_[i]->has_pure();

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
      ExEnv::err0() << indent
                    << scprintf("GaussianBasisSet::count_shells_ couldn't find \"%s\":\n", keyword);
      keyval->errortrace(ExEnv::err0());
      throw std::runtime_error("GaussianBasisSet::count_shells_ -- couldn't find the basis set");
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
GaussianBasisSet::get_shells_(int& ishell, Ref<KeyVal>& keyval, const char* element, const char* basisname, BasisFileSet& bases,
			      int havepure, int pure, bool missing_ok)
{
  char keyword[KeyVal::MaxKeywordLength];

  sprintf(keyword,":basis:%s:%s:am",
          element,basisname);
  if (keyval->exists(keyword)) {
    get_even_temp_shells_(ishell, keyval, element, basisname, havepure, pure);
  }
  else {
    recursively_get_shell(ishell, keyval, element, basisname,
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
GaussianBasisSet::get_even_temp_shells_(int& ishell, Ref<KeyVal>& keyval, const char* element, const char* basisname,
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
      int* am = new int[1];
      double* exps = new double[1];
      double** coeffs = new double*[1];
      coeffs[0] = new double[1];

      exps[0] = alpha;
      am[0] = l;
      coeffs[0][0] = 1.0;

      if (l < 1)
        shell_[ishell] = new GaussianShell(1,1,exps,am,GaussianShell::Cartesian,coeffs,GaussianShell::Normalized);
      else if (havepure)
        shell_[ishell] = new GaussianShell(1,1,exps,am,GaussianShell::Pure,coeffs,GaussianShell::Normalized);
      else
        shell_[ishell] = new GaussianShell(1,1,exps,am,GaussianShell::Cartesian,coeffs,GaussianShell::Normalized);
      ishell++;

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
  recursively_get_shell(int&ishell,Ref<KeyVal>&keyval,
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
      ExEnv::err0() << indent
           << scprintf("GaussianBasisSet:: couldn't find \"%s\":\n", keyword);
      keyval->errortrace(ExEnv::err0());
      throw std::runtime_error("GaussianBasisSet::recursively_get_shell -- couldn't find the basis set");
    }
  if (!count) return;
  for (int j=0; j<count; j++) {
      sprintf(prefix,":basis:%s:%s",
	      element,basisname);
      Ref<KeyVal> prefixkeyval = new PrefixKeyVal(keyval,prefix,j);
      if (prefixkeyval->exists("get")) {
          char* newbasis = prefixkeyval->pcharvalue("get");
          if (!newbasis) {
	      ExEnv::err0() << indent << "GaussianBasisSet: "
                   << scprintf("error processing get for \"%s\"\n", prefix);
              keyval->errortrace(ExEnv::err0());
	      exit(1);
	    }
	  recursively_get_shell(ishell,keyval,element,newbasis,bases,
                                havepure,pure,get,missing_ok);
          delete[] newbasis;
	}
      else {
          if (get) {
	      if (havepure) shell_[ishell] = new GaussianShell(prefixkeyval,pure);
	      else shell_[ishell] = new GaussianShell(prefixkeyval);
	    }
	  ishell++;
	}
    }
}

GaussianBasisSet::~GaussianBasisSet()
{
  delete[] name_;
  delete[] label_;

  int ii;
  for (ii=0; ii<nshell_; ii++) {
      delete shell_[ii];
    }
  delete[] shell_;
}

int
GaussianBasisSet::max_nfunction_in_shell() const
{
  int i;
  int max = 0;
  for (i=0; i<nshell_; i++) {
      if (max < shell_[i]->nfunction()) max = shell_[i]->nfunction();
    }
  return max;
}

int
GaussianBasisSet::max_ncontraction() const
{
  int i;
  int max = 0;
  for (i=0; i<nshell_; i++) {
      if (max < shell_[i]->ncontraction()) max = shell_[i]->ncontraction();
    }
  return max;
}

int
GaussianBasisSet::max_angular_momentum() const
{
  int i;
  int max = 0;
  for (i=0; i<nshell_; i++) {
      int maxshi = shell_[i]->max_angular_momentum();
      if (max < maxshi) max = maxshi;
    }
  return max;
}

int
GaussianBasisSet::max_cartesian() const
{
  int i;
  int max = 0;
  for (i=0; i<nshell_; i++) {
      int maxshi = shell_[i]->max_cartesian();
      if (max < maxshi) max = maxshi;
    }
  return max;
}

int
GaussianBasisSet::max_ncartesian_in_shell(int aminc) const
{
  int i;
  int max = 0;
  for (i=0; i<nshell_; i++) {
      int maxshi = shell_[i]->ncartesian_with_aminc(aminc);
      if (max < maxshi) max = maxshi;
    }
  return max;
}

int
GaussianBasisSet::max_nprimitive_in_shell() const
{
  int i;
  int max = 0;
  for (i=0; i<nshell_; i++) {
      if (max < shell_[i]->nprimitive()) max = shell_[i]->nprimitive();
    }
  return max;
}

int
GaussianBasisSet::max_am_for_contraction(int con) const
{
  int i;
  int max = 0;
  for (i=0; i<nshell_; i++) {
      if (shell_[i]->ncontraction() <= con) continue;
      int maxshi = shell_[i]->am(con);
      if (max < maxshi) max = maxshi;
    }
  return max;
}

int
GaussianBasisSet::function_to_shell(int func) const
{
  return function_to_shell_[func];
}

int
GaussianBasisSet::ncenter() const
{
  return ncenter_;
}

int
GaussianBasisSet::nshell_on_center(int icenter) const
{
  return center_to_nshell_[icenter];
}

int
GaussianBasisSet::nbasis_on_center(int icenter) const
{
  return center_to_nbasis_[icenter];
}

int
GaussianBasisSet::shell_on_center(int icenter, int ishell) const
{
  return center_to_shell_[icenter] + ishell;
}




const GaussianShell&
GaussianBasisSet::operator()(int icenter,int ishell) const
{
  return *shell_[center_to_shell_[icenter] + ishell];
}

GaussianShell&
GaussianBasisSet::operator()(int icenter,int ishell)
{
  return *shell_[center_to_shell_[icenter] + ishell];
}

int
GaussianBasisSet::equiv(const Ref<GaussianBasisSet> &b)
{
  if (nshell() != b->nshell()) return 0;
  for (int i=0; i<nshell(); i++) {
      if (!shell_[i]->equiv(b->shell_[i])) return 0;
    }
  return 1;
}

void
GaussianBasisSet::print_brief(ostream& os) const
{
  os << indent
     << "GaussianBasisSet:" << endl << incindent
     << indent << "nbasis = " << nbasis_ << endl
     << indent << "nshell = " << nshell_ << endl
     << indent << "nprim  = " << nprim_ << endl;
  if (name_) {
      os << indent
         << "name = \"" << name_ << "\"" << endl;
  }
  else {
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
  for (icenter=0; icenter < ncenter_; icenter++) {
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
             << scprintf("Shell %d: functionnum = %d, primnum = %d\n",
                         ishell,shell_to_function_[ioshell],shell_to_primitive_[ioshell]);
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
GaussianBasisSet::init(
    char *name,
    char *label,
    const Ref<Molecule> &molecule,
    const Ref<SCMatrixKit> &matrixkit,
    const Ref<SCMatrixKit> &so_matrixkit,
    GaussianShell **shell,
    const std::vector<int> shell_to_center)
    
{
  name_ = name;
  label_ = label;
  molecule_ = molecule;
  matrixkit_ = matrixkit;
  so_matrixkit_ = so_matrixkit;
  shell_ = shell;
  shell_to_center_ = shell_to_center;

  ncenter_ = molecule_->natom();
  nshell_ = shell_to_center_.size();

  nbasis_ = 0;
  nprim_  = 0;
  has_pure_ = false;

  shell_to_function_.resize(nshell_);
  shell_to_primitive_.resize(nshell_);
  center_to_nshell_.resize(ncenter_);
  center_to_nbasis_.resize(ncenter_);
  center_to_shell_.resize(nshell_);

  std::fill(center_to_shell_.begin(), center_to_shell_.end(), -1);
  std::fill(center_to_nshell_.begin(), center_to_nshell_.end(), 0);
  std::fill(center_to_nbasis_.begin(), center_to_nbasis_.end(), 0);

  for (int ishell=0; ishell<nshell_; ishell++) {
      int center = shell_to_center_[ishell];
      if (center_to_shell_[center] == -1) {
          center_to_shell_[center] = ishell;
        }
      center_to_nshell_[center]++;
      center_to_nbasis_[center] += shell_[ishell]->nfunction();
      shell_to_function_[ishell] = nbasis_;
      shell_to_primitive_[ishell] = nprim_;
      shell_to_center_[ishell] = center;
      nbasis_ += shell_[ishell]->nfunction();
      nprim_ += shell_[ishell]->nprimitive();
      if (shell_[ishell]->has_pure()) has_pure_ = true;
    }

  if (basisdim_.null()) {
    int *bs = new int[nshell_];
    for (int s=0; s < nshell_; s++)
      bs[s] = shell_[s]->nfunction();
    basisdim_ = new SCDimension(nbasis_, nshell_, bs, "basis set dimension");
    delete[] bs;
  }

  function_to_shell_.resize(nbasis_);
  int ifunc = 0;
  for (int i=0; i<nshell_; i++) {
      int nfun = shell_[i]->nfunction();
      for (int j=0; j<nfun; j++) {
          function_to_shell_[ifunc] = i;
          ifunc++;
        }
    }
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
      if (i>1) sivec_[i] = integral->new_spherical_transform_iter(i);
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

/////////////////////////////////////////////////////////////////////////////

static ClassDesc GaussianBasisSetSum_cd(
  typeid(GaussianBasisSetSum),"GaussianBasisSetSum",1,
  "virtual public SavableState",
  0, 0, create<GaussianBasisSetSum>);


GaussianBasisSetSum::GaussianBasisSetSum(const Ref<GaussianBasisSet>& bs1,
                                         const Ref<GaussianBasisSet>& bs2) :
  bs1_(bs1), bs2_(bs2)
{
  sum(bs1,bs2);
}

GaussianBasisSetSum::GaussianBasisSetSum(StateIn&s):
  SavableState(s)
{
  bs1_ << SavableState::restore_state(s);
  bs2_ << SavableState::restore_state(s);
  bs12_ << SavableState::restore_state(s);

  // maps
  s.get(shell_to_basis_);
  s.get(function_to_basis_);
  s.get(fblock_to_function_);
  s.get(fblock_size_);
}

GaussianBasisSetSum::~GaussianBasisSetSum()
{
}

void
GaussianBasisSetSum::save_data_state(StateOut&s)
{
  SavableState::save_state(bs1_.pointer(),s);
  SavableState::save_state(bs2_.pointer(),s);
  SavableState::save_state(bs12_.pointer(),s);
  
  // maps
  s.put(shell_to_basis_);
  s.put(function_to_basis_);
  s.put(fblock_to_function_);
  s.put(fblock_size_);
}

void
GaussianBasisSetSum::sum(const Ref<GaussianBasisSet>& A,
                         const Ref<GaussianBasisSet>& B)
{
  GaussianBasisSet* a = A.pointer();
  GaussianBasisSet* b = B.pointer();
  if (a->molecule_.pointer() != b->molecule_.pointer())
  throw std::runtime_error("GaussianBasisSetSum::sum -- cannot sum basis sets, molecules are different");

  Ref<Molecule> molecule = a->molecule();
  const int ncenter = a->ncenter();
  Ref<SCMatrixKit> matrixkit = a->matrixkit();
  const int nshell = a->nshell() + b->nshell();
  std::vector<int> center_to_nshell(ncenter);
  shell_to_basis_.resize(nshell);

  GaussianShell** shell = new GaussianShell*[nshell];
  int* func_per_shell = new int[nshell];

  for(int c=0; c<ncenter; c++) {

    int ns1 = a->nshell_on_center(c);
    int ns2 = b->nshell_on_center(c);
    int ns = ns1+ns2;
    int s1off = a->shell_on_center(c,0);
    int s2off = b->shell_on_center(c,0);
    int soff = s1off + s2off;
    center_to_nshell[c] = ns;

    // i is the shell index on this center in the composite basis
    // ii is the absolute index
    for (int i=0; i<ns; i++) {
      const int ii = soff + i;
      const GaussianShell* gsi;
      if (i < ns1) {
        gsi = &a->shell(s1off + i);
        shell_to_basis_[ii] = 1;
      }
      else {
        gsi = &b->shell(s2off + i - ns1);
        shell_to_basis_[ii] = 1;
      }

      int nc=gsi->ncontraction();
      int np=gsi->nprimitive();
      func_per_shell[ii] = gsi->nfunction();

      int *ams = new int[nc];
      int *pure = new int[nc];
      double *exps = new double[np];
      double **coefs = new double*[nc];

      for (int j=0; j < nc; j++) {
        ams[j] = gsi->am(j);
        pure[j] = gsi->is_pure(j);
        coefs[j] = new double[np];
        for (int k=0; k < np; k++)
        coefs[j][k] = gsi->coefficient_unnorm(j,k);
      }

      for (int j=0; j < np; j++)
      exps[j] = gsi->exponent(j);

      shell[soff + i] = new GaussianShell(nc, np, exps, ams, pure, coefs,
          GaussianShell::Unnormalized);
    }
  }

  int nbas = a->nbasis() + b->nbasis();
  RefSCDimension basisdim = new SCDimension(nbas, nshell, func_per_shell, "basis set dimension");

  const char* A_name = a->name();
  const char* B_name = b->name();
  const char* AplusB_name = 0;
  if (!A_name && !B_name) {
    ostringstream oss;
    oss << "[" << A_name << "]+[" << B_name << "]";
    std::string tmpname = oss.str();
    AplusB_name = strcpy(new char[tmpname.size()+1],tmpname.c_str());
  }
  const char* AplusB_label = 0;
  if (AplusB_name) {
    AplusB_label = AplusB_name;
  }
  else {
    ostringstream oss;
    const char* A_label = a->label();
    const char* B_label = b->label();
    oss << "[" << A_label << "]+[" << B_label << "]";
    std::string tmpname = oss.str();
    AplusB_label = strcpy(new char[tmpname.size()+1],tmpname.c_str());
  }

  bs12_ = new GaussianBasisSet(AplusB_name, AplusB_label, molecule,
      matrixkit, basisdim, ncenter,
      nshell, shell, center_to_nshell);

  delete[] func_per_shell;
  delete[] AplusB_name;
  if (AplusB_name != AplusB_label)
  delete[] AplusB_label;

  //
  // compute the rest of the maps using shell_to_basis_
  //
  const int nbasis = a->nbasis() + b->nbasis();
  function_to_basis_.resize(nbasis);
  for(int s=0; s<nshell; ++s) {
    const int bs = shell_to_basis_[s];
    const int nf = bs12_->shell(s).nfunction();
    const int foff = bs12_->shell_to_function(s);
    for(int f=0; f<nf; ++f)
      function_to_basis_[f + foff] = bs;
  }
  {
    int curr_basis = function_to_basis_[0];
    fblock_to_function_.push_back(0);
    int nf_in_curr_fblock = 1;
    for(int f=1; f<nbasis; ++f) {
      const int basis = function_to_basis_[f];
      if (basis == curr_basis) {
        // still in the current fblock
        ++nf_in_curr_fblock;
      }
      else {
        // finalize this fblock
        curr_basis = basis;
        fblock_size_.push_back(nf_in_curr_fblock);
        nf_in_curr_fblock = 0;
        // move onto the new fblock
        fblock_to_function_.push_back(f);
      }
    }
    // finalize the last fblock
    fblock_size_.push_back(nf_in_curr_fblock);
  }
}

const Ref<GaussianBasisSet>&
GaussianBasisSetSum::bs1() const { return bs1_; }

const Ref<GaussianBasisSet>&
GaussianBasisSetSum::bs2() const { return bs2_; }

const Ref<GaussianBasisSet>&
GaussianBasisSetSum::bs12() const { return bs12_; }



/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:

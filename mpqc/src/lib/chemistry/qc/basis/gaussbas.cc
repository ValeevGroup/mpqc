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
#include <strstream.h>

#include <util/keyval/keyval.h>
#include <util/misc/formio.h>
#include <util/misc/newstring.h>
#include <math/scmat/blocked.h>
#include <chemistry/molecule/molecule.h>
#include <chemistry/qc/basis/gaussshell.h>
#include <chemistry/qc/basis/gaussbas.h>
#include <chemistry/qc/basis/files.h>
#include <chemistry/qc/basis/cartiter.h>
#include <chemistry/qc/basis/transform.h>
#include <chemistry/qc/basis/integral.h>

SavableState_REF_def(GaussianBasisSet);

#define CLASSNAME GaussianBasisSet
#define PARENTS public SavableState
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#define VERSION 2
#include <util/state/statei.h>
#include <util/class/classi.h>

void *
GaussianBasisSet::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SavableState::_castdown(cd);
  return do_castdowns(casts,cd);
}

GaussianBasisSet::GaussianBasisSet(const RefKeyVal&topkeyval)
{
  molecule_ = topkeyval->describedclassvalue("molecule");
  if (molecule_.null()) {
      cerr << node0 << indent << "GaussianBasisSet: no \"molecule\"\n";
      abort();
    }

  // see if the user requests pure am or cartesian functions
  int pure;
  pure = topkeyval->booleanvalue("puream");
  if (topkeyval->error() != KeyVal::OK) pure = -1;

  // construct a keyval that contains the basis library

  // this ParsedKeyVal CTOR looks at the basisdir and basisfiles
  // variables to find out what basis set files are to be read in
  // the files are read the on only node 0
  RefMessageGrp grp = MessageGrp::get_default_messagegrp();
  RefParsedKeyVal parsedkv = new ParsedKeyVal();
  char *in_char_array;
  if (grp->me() == 0) {
      ostrstream ostrs;
      ParsedKeyVal::cat_files("basis",topkeyval,ostrs);
      ostrs << ends;
      in_char_array = ostrs.str();
      int n = ostrs.pcount();
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
  RefKeyVal libkeyval = parsedkv.pointer();
  RefKeyVal keyval = new AggregateKeyVal(topkeyval,libkeyval);

  // if there isn't a matrixkit in the input, init2() will assign the
  // default matrixkit
  matrixkit_ = keyval->describedclassvalue("matrixkit");
  
  // Bases keeps track of what basis set data bases have already
  // been read in.  It also handles the conversion of basis
  // names to file names.
  BasisFileSet bases(keyval);
  init(molecule_,keyval,bases,1,pure);
}

GaussianBasisSet::GaussianBasisSet(const GaussianBasisSet& gbs) :
  nshell_(gbs.nshell_),
  ncenter_(gbs.ncenter_),
  molecule_(gbs.molecule_),
  matrixkit_(gbs.matrixkit_),
  basisdim_(gbs.basisdim_)
{
  int i,j;
  
  name_ = new_string(gbs.name_);

  center_to_nshell_.set_length(ncenter_);
  for (i=0; i < ncenter_; i++) {
      center_to_nshell_(i) = gbs.center_to_nshell_(i);
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

GaussianBasisSet::GaussianBasisSet(StateIn&s):
  SavableState(s),
  center_to_nshell_(s)
{
  matrixkit_ = SCMatrixKit::default_matrixkit();

  molecule_.restore_state(s);
  basisdim_.restore_state(s);

  ncenter_ = center_to_nshell_.length();
  s.getstring(name_);

  nshell_ = 0;
  int i;
  for (i=0; i<ncenter_; i++) {
      nshell_ += center_to_nshell_(i);
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
  center_to_nshell_.save_object_state(s);

  molecule_.save_state(s);
  basisdim_.save_state(s);
  
  s.putstring(name_);
  for (int i=0; i<nshell_; i++) {
      shell_[i]->save_object_state(s);
    }
}

void
GaussianBasisSet::init(RefMolecule&molecule,
                       RefKeyVal&keyval,
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
          cerr << node0 << indent
               << "GaussianBasisSet: No name given for basis set\n";
          abort();
        }
    }

  // construct prefixes for each atom: :basis:<atom>:<basisname>:<shell#>
  // and read in the shell
  nbasis_ = 0;
  int ishell = 0;
  ncenter_ = molecule->natom();
  int iatom;
  for (iatom=0; iatom < ncenter_; iatom++) {
      ChemicalElement currentelement(molecule->operator[](iatom).element());
      // see if there is a specific basisname for this atom
      char* sbasisname = 0;
      if (have_custom && !nelement) {
          sbasisname = keyval->pcharvalue("basis",iatom);
        }
      else if (nelement) {
          int i;
          for (i=0; i<nelement; i++) {
              char *tmpelementname = keyval->pcharvalue("element", i);
              ChemicalElement tmpelement(tmpelementname);
              if (tmpelement == currentelement) {
                  sbasisname = keyval->pcharvalue("basis", i);
                  break;
                }
              delete[] tmpelementname;
            }
        }
      if (!sbasisname) {
          if (!name_) {
              cerr << node0 << indent << "GaussianBasisSet: "
                   << scprintf("no basis name for atom %d (%s)\n",
                               iatom, currentelement.name());
              abort();
            }
          sbasisname = strcpy(new char[strlen(name_)+1],name_);
        }
      recursively_get_shell(ishell,keyval,
                            currentelement.name(),
                            sbasisname,bases,havepure,pure,0);
      delete[] sbasisname;
    }
  nshell_ = ishell;
  shell_ = new GaussianShell*[nshell_];
  ishell = 0;
  center_to_nshell_.set_length(ncenter_);
  for (iatom=0; iatom<ncenter_; iatom++) {
      ChemicalElement currentelement(molecule->operator[](iatom).element());
      // see if there is a specific basisname for this atom
      char* sbasisname = 0;
      if (have_custom && !nelement) {
          sbasisname = keyval->pcharvalue("basis",iatom);
        }
      else if (nelement) {
          int i;
          for (i=0; i<nelement; i++) {
              char *tmpelementname = keyval->pcharvalue("element", i);
              ChemicalElement tmpelement(tmpelementname);
              if (tmpelement == currentelement) {
                  sbasisname = keyval->pcharvalue("basis", i);
                  break;
                }
              delete[] tmpelementname;
            }
        }
      if (!sbasisname) {
          if (!name_) {
              cerr << node0 << indent << "GaussianBasisSet: "
                   << scprintf("no basis name for atom %d (%s)\n",
                               iatom, currentelement.name());
              abort();
            }
          sbasisname = strcpy(new char[strlen(name_)+1],name_);
        }

      int ishell_old = ishell;
      recursively_get_shell(ishell,keyval,
                            currentelement.name(),
                            sbasisname,bases,havepure,pure,1);

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
  init2();
}

double
GaussianBasisSet::r(int icenter, int xyz) const
{
  return molecule_->atom(icenter).point()[xyz];
}

void
GaussianBasisSet::init2()
{
  // center_to_shell_ and shell_to_center_
  shell_to_center_.set_length(nshell_);
  center_to_shell_.set_length(ncenter_);
  center_to_nbasis_.set_length(ncenter_);
  int ishell = 0;
  for (int icenter=0; icenter<ncenter_; icenter++) {
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
  shell_to_function_.set_length(nshell_);
  nbasis_ = 0;
  nprim_ = 0;
  for (ishell=0; ishell<nshell_; ishell++) {
      shell_to_function_[ishell] = nbasis_;
      nbasis_ += shell_[ishell]->nfunction();
      nprim_ += shell_[ishell]->nprimitive();
    }

  // would like to do this in function_to_shell(), but it is const
  int n = nbasis();
  int nsh = nshell();
  function_to_shell_.set_length(n);
  int ifunc = 0;
  for (int i=0; i<nsh; i++) {
      int nfun = operator()(i).nfunction();
      for (int j=0; j<nfun; j++) {
          function_to_shell_[ifunc] = i;
          ifunc++;
        }
    }

  if (matrixkit_.null())
    matrixkit_ = SCMatrixKit::default_matrixkit();

  so_matrixkit_ = new BlockedSCMatrixKit(matrixkit_);
  
  if (basisdim_.null()) {
    int nb = nshell();
    int *bs = new int[nb];
    for (int s=0; s < nb; s++)
      bs[s] = shell(s).nfunction();
    basisdim_ = new SCDimension(nbasis(), nb, bs, "basis set dimension");
  }

  civec_ = 0;
  sivec_ = 0;
}

void
GaussianBasisSet::set_matrixkit(const RefSCMatrixKit& mk)
{
  matrixkit_ = mk;
  so_matrixkit_ = new BlockedSCMatrixKit(matrixkit_);
}

void
GaussianBasisSet::
  recursively_get_shell(int&ishell,RefKeyVal&keyval,
			const char*element,
			const char*basisname,
                        BasisFileSet &bases,
			int havepure,int pure,
			int get)
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
      cerr << node0 << indent
           << scprintf("GaussianBasisSet:: couldn't find \"%s\":\n", keyword);
      keyval->errortrace(cerr);
      exit(1);
    }
  if (!count) return;
  for (int j=0; j<count; j++) {
      sprintf(prefix,":basis:%s:%s",
	      element,basisname);
      RefKeyVal prefixkeyval = new PrefixKeyVal(prefix,keyval,j);
      if (prefixkeyval->exists("get")) {
          char* newbasis = prefixkeyval->pcharvalue("get");
          if (!newbasis) {
	      cerr << node0 << indent << "GaussianBasisSet: "
                   << scprintf("error processing get for \"%s\"\n", prefix);
              keyval->errortrace(cerr);
	      exit(1);
	    }
	  recursively_get_shell(ishell,keyval,element,newbasis,bases,
                                havepure,pure,get);
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
  RefIntegral nullint;
  set_integral(nullint);

  delete[] name_;

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
  return center_to_shell_(icenter) + ishell;
}
int
GaussianBasisSet::shell_to_center(int ishell) const
{
  return shell_to_center_(ishell);
}

const GaussianShell&
GaussianBasisSet::operator()(int icenter,int ishell) const
{
  return *shell_[center_to_shell_(icenter) + ishell];
}

GaussianShell&
GaussianBasisSet::operator()(int icenter,int ishell)
{
  return *shell_[center_to_shell_(icenter) + ishell];
}

void
GaussianBasisSet::set_integral(const RefIntegral &integral)
{
  int i;
  int maxam = max_angular_momentum();
  if (civec_) {
      for (i=0; i<=maxam; i++) {
          delete civec_[i];
          delete sivec_[i];
        }
      delete[] civec_;
      delete[] sivec_;
      civec_ = 0;
      sivec_ = 0;
    }
  if (integral.nonnull()) {
      civec_ = new CartesianIter *[maxam+1];
      sivec_ = new SphericalTransformIter *[maxam+1];
      for (i=0; i<=maxam; i++) {
          civec_[i] = integral->new_cartesian_iter(i);
          if (i>1) sivec_[i] = integral->new_spherical_transform_iter(i);
          else sivec_[i] = 0;
        }
    }
}

void
GaussianBasisSet::print_brief(ostream& os) const
{
  os << node0 << indent
     << "GaussianBasisSet:" << endl << incindent
     << indent << "nbasis = " << nbasis_ << endl
     << indent << "nshell = " << nshell_ << endl
     << indent << "nprim  = " << nprim_ << endl;
  if (name_) {
      os << node0 << indent
         << "name = \"" << name_ << "\"" << endl;
    }
  os << decindent;
}

void
GaussianBasisSet::print(ostream& os) const
{
  print_brief(os);
  os << incindent;

  // Loop over centers
  int icenter = 0;
  int ioshell = 0;
  for (icenter=0; icenter < ncenter_; icenter++) {
      os << node0 << endl << indent
         << scprintf(
             "center %d: %12.8f %12.8f %12.8f, nshell = %d, shellnum = %d\n",
             icenter,
             r(icenter,0),
             r(icenter,1),
             r(icenter,2),
             center_to_nshell_[icenter],
             center_to_shell_[icenter]);
      for (int ishell=0; ishell < center_to_nshell_[icenter]; ishell++) {
	  os << node0 << indent
             << scprintf("Shell %d: functionnum = %d\n",
                         ishell,shell_to_function_[ioshell]);
          os << node0 << incindent;
	  operator()(icenter,ishell).print(os);
          os << node0 << decindent;
          ioshell++;
	}
    }

  os << node0 << decindent;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:

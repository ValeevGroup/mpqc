//
// uncontract.cc
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

#include <util/class/scexception.h>
#include <chemistry/qc/basis/uncontract.h>
#include <chemistry/qc/basis/gaussshell.h>

using namespace std;
using namespace sc;

static ClassDesc UncontractedBasisSet_cd(
  typeid(UncontractedBasisSet),"UncontractedBasisSet",1,
  "public GaussianBasisSet",
  0, create<UncontractedBasisSet>, create<UncontractedBasisSet>);

UncontractedBasisSet::UncontractedBasisSet(const Ref<KeyVal>&keyval)
{
  Ref<GaussianBasisSet> basis;
  basis << keyval->describedclassvalue("basis");
  if (basis.null()) {
    basis = new GaussianBasisSet(keyval);
    if (basis.null())
      throw InputError("could not construct a GaussianBasisSet",
                       __FILE__, __LINE__,
                       "basis", 0, class_desc());
    }

  uncontract(basis);
}

UncontractedBasisSet::UncontractedBasisSet(const Ref<GaussianBasisSet>&b)
{
  uncontract(b);
}

UncontractedBasisSet::UncontractedBasisSet(StateIn&s):
  SavableState(s),
  GaussianBasisSet(s)
{
}

void
UncontractedBasisSet::save_data_state(StateOut&s)
{
  GaussianBasisSet::save_data_state(s);
}

static
char *
name_conv(const char *name)
{
  if (name == 0) return 0;
  std::string newname = "Uncontracted(";
  newname += name;
  newname += ")";
  return strcpy(new char[newname.size()+1],newname.c_str());
}

// using to order shells in terms of:
//   increasing center number
//   increasing angular momentum
//   decreasing exponent
struct shelldesc_t {
    int center;
    int am;
    double exponent;

    // we do not sort on pure
    mutable bool pure;

    bool operator<(const shelldesc_t &a) const {
      if (center < a.center) return true;
      else if (center > a.center) return false;
      if (am < a.am) return true;
      else if (am > a.am) return false;
      if (exponent > a.exponent) return true;
      return false;
    }
};

void
UncontractedBasisSet::uncontract(const Ref<GaussianBasisSet>&basis)
{
  std::set<shelldesc_t> shellinfo;

  shelldesc_t current;
  for (int i=0; i<basis->ncenter(); i++) {
      current.center = i;
      for (int j=0; j<basis->nshell_on_center(i); j++) {
          GaussianShell &shell = basis->shell(i,j);
          for (int con=0; con<shell.ncontraction(); con++) {
              current.am = shell.am(con);
              current.pure = shell.is_pure(con);
              for (int prim=0; prim<shell.nprimitive(); prim++) {
                  current.exponent = shell.exponent(prim);
                  if (shell.coefficient_unnorm(con,prim) != 0.0) {
                      // insert the current shell into the set
                      std::pair<std::set<shelldesc_t>::const_iterator,bool>
                          iter_inserted = shellinfo.insert(current);
                      // nonpure shells completely span pure shells,
                      // so if both are given, make sure the nonpure
                      // shell ends up in the shellinfo set.
                      if (!current.pure && iter_inserted.first->pure) {
                          iter_inserted.first->pure = false;
                        }
                    }
                }
            }
        }
    }

  // convert the set shell structures

  int nshell = shellinfo.size();
  GaussianShell **shell = new GaussianShell*[nshell];
  int ishell = 0;

  std::vector<int> shell_to_center(nshell);

  for (std::set<shelldesc_t>::iterator iter = shellinfo.begin();
       iter != shellinfo.end(); iter++) {
      double *exponents = new double[1];
      int *am = new int[1];
      *exponents = iter->exponent;
      *am = iter->am;
      double **c = new double*[1];
      *c = new double[1];
      **c = 1.0;
      int *pure = new int[1];
      if (iter->pure) *pure = 1;
      else            *pure = 0;
      shell[ishell] = new GaussianShell(1, 1, exponents, am,
                                        pure, c);
      shell_to_center[ishell] = iter->center;
      ishell++;
    }

  init(name_conv(basis->name()),
       name_conv(basis->label()),
       basis->molecule(),
       basis->matrixkit(),
       basis->so_matrixkit(),
       shell,
       shell_to_center);
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:

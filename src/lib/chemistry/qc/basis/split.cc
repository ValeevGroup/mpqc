//
// split.cc
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

#include <util/misc/scexception.h>
#include <chemistry/qc/basis/split.h>
#include <chemistry/qc/basis/gaussshell.h>

using namespace std;
using namespace sc;

static ClassDesc SplitBasisSet_cd(
  typeid(SplitBasisSet),"SplitBasisSet",1,
  "public GaussianBasisSet",
  0, create<SplitBasisSet>, create<SplitBasisSet>);

SplitBasisSet::SplitBasisSet(const Ref<KeyVal>&keyval)
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

  split(basis);
}

SplitBasisSet::SplitBasisSet(const Ref<GaussianBasisSet>&b,
                             std::string name)
{
  split(b, name);
}

SplitBasisSet::SplitBasisSet(StateIn&s):
  SavableState(s),
  GaussianBasisSet(s)
{
}

void
SplitBasisSet::save_data_state(StateOut&s)
{
  GaussianBasisSet::save_data_state(s);
}

static
char *
name_conv(const char *name)
{
  if (name == 0) return 0;
  std::string newname = "Split(";
  newname += name;
  newname += ")";
  return strcpy(new char[newname.size()+1],newname.c_str());
}

void
SplitBasisSet::split(const Ref<GaussianBasisSet>&basis,
                     std::string name)
{
  std::vector<int> nshell_on_center(basis->ncenter());
  std::fill(nshell_on_center.begin(), nshell_on_center.end(), 0);

  int nshell = 0;
  for (int icenter=0; icenter<basis->ncenter(); icenter++) {
      for (int ishell=0; ishell<basis->nshell_on_center(icenter);
           ishell++) {
          nshell_on_center[icenter]
              += basis->shell(icenter,ishell).ncontraction();
        }
      nshell += nshell_on_center[icenter];
    }

  GaussianShell **shells = new GaussianShell*[nshell];
  int ishell = 0;

  std::vector<int> shell_to_center(nshell);

  for (int icenter=0, ishellall=0; icenter<basis->ncenter(); icenter++) {
      for (int ishell=0; ishell<basis->nshell_on_center(icenter);
           ishell++) {
          const GaussianShell &shell = basis->shell(icenter,ishell);
          int ncon = shell.ncontraction();
          int nprim = shell.nprimitive();
          for (int icon=0; icon<ncon; icon++, ishellall++) {
              shell_to_center[ishellall] = icenter;
              int *am = new int[1];
              int *pure = new int[1];
              *am = shell.am(icon);
              *pure = shell.is_pure(icon);
              double *exponents = new double[nprim];
              double **c = new double*[1];
              *c = new double[nprim];
              for (int iprim=0; iprim<nprim; iprim++) {
                  exponents[iprim] = shell.exponent(iprim);
                  c[0][iprim] = shell.coefficient_unnorm(icon,iprim);
                }
              shells[ishellall] // is hell all???
                  = new GaussianShell(1, nprim, exponents, am,
                                      pure, c,
                                      GaussianShell::Unnormalized);
            }
        }
    }

  init((name.empty() ? name_conv(basis->name()) : const_cast<char*>(name.c_str())),
       name_conv(basis->label()),
       basis->molecule(),
       basis->matrixkit(),
       basis->so_matrixkit(),
       shells,
       shell_to_center);

//   SCFormIO::setverbose(ExEnv::out0(), 1);
//   basis->print();
//   print();
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:

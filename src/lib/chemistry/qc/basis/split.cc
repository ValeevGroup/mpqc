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
std::string
name_conv(const std::string& name)
{
  if (name.empty()) return name;
  std::string newname = "Split(";
  newname += name;
  newname += ")";
  return newname;
}

namespace {
  struct split_filter {
      split_filter(unsigned int contr) : contr_(contr) {}
      bool operator()(const GaussianShell& shell,
                      unsigned int contr) {
        return contr == contr_;
      }
      unsigned int contr_;
  };
}

void
SplitBasisSet::split(const Ref<GaussianBasisSet>&basis,
                     std::string name)
{
  molecule_ = basis->molecule();

  // create shells
  std::vector<Shell> shells;
  for (int s = 0; s < basis->nshell(); ++s) {
    const GaussianShell &shell = basis->shell(s);

    for(unsigned int c=0; c<shell.ncontraction(); ++c) {
      split_filter f(c);
      shells.push_back(Shell(this, basis->shell_to_center(s), filter(shell, f)));
    }
  }

  init((name.empty() ? name_conv(basis->name()) : const_cast<char*>(name.c_str())),
       name_conv(basis->label()),
       basis->molecule(),
       shells);

  if (debug()) {
    SCFormIO::setverbose(ExEnv::out0(), 1);
    basis->print();
    print();
    SCFormIO::setverbose(ExEnv::out0(), 0);
  }
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:

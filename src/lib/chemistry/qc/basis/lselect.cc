//
// lselect.cc
//
// Copyright (C) 2007 Edward Valeev
//
// Author: Edward Valeev <evaleev@vt.edu>
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

#include <algorithm>
#include <util/state/statein.h>
#include <util/state/stateout.h>
#include <util/misc/scexception.h>
#include <chemistry/qc/basis/lselect.h>
#include <chemistry/qc/basis/split.h>
#include <chemistry/qc/basis/gaussshell.h>

using namespace std;
using namespace sc;

static ClassDesc LSelectBasisSet_cd(
  typeid(LSelectBasisSet),"LSelectBasisSet",1,
  "public GaussianBasisSet",
  0, create<LSelectBasisSet>, create<LSelectBasisSet>);

LSelectBasisSet::LSelectBasisSet(const Ref<KeyVal>&keyval)
{
  Ref<GaussianBasisSet> basis;
  basis << keyval->describedclassvalue("basis");
  if (basis == 0) {
    basis = new GaussianBasisSet(keyval);
    if (basis == 0)
      throw InputError("could not construct a GaussianBasisSet",
                       __FILE__, __LINE__,
                       "basis", 0, class_desc());
  }

  if (keyval->exists("l")) {
      int nl = keyval->count("l");
      if (nl < 1)
	  throw InputError("l must be an array of at least 1 element",__FILE__, __LINE__);
      for(int c=0; c<nl; ++c) {
        const int l = keyval->intvalue("l",c);
        if (l < 0) {
          throw InputError("l < 0", __FILE__, __LINE__);
        }
        l_.push_back(l);
      }
  }
  else {
      int lmin = keyval->intvalue("lmin",KeyValValueint(0));
      if (lmin < 0) lmin = 0;
      int lmax = keyval->intvalue("lmax",KeyValValueint(basis->max_angular_momentum()));
      if (lmin > lmax)
        throw InputError("lmin > lmax",__FILE__, __LINE__);
      if (lmin < 0)
        throw InputError("lmin < 0",__FILE__, __LINE__);
      if (lmax < 0)
        throw InputError("lmax < 0",__FILE__, __LINE__);
      for(int l=lmin; l<=lmax; ++l)
        l_.push_back(l);
  }

  lselect(basis);
}

LSelectBasisSet::LSelectBasisSet(StateIn&s):
  SavableState(s),
  GaussianBasisSet(s)
{
    s.get(l_);
}

void
LSelectBasisSet::save_data_state(StateOut&s)
{
  GaussianBasisSet::save_data_state(s);
  s.put(l_);
}

namespace {
    std::string
    name_conv(const std::string& name, const std::vector<unsigned int>& l_)
    {
	std::ostringstream oss;
	if (name.empty()) return name;
	oss << "LSelect(" << name << ":";
	for(std::vector<unsigned int>::const_iterator am=l_.begin(); am!=l_.end(); ++am)
	    oss << GaussianShell::amtypes[*am];
	oss << ")";
	const std::string& newname = oss.str();
	return newname;
    }

}

namespace {
  struct lselect_filter {
      lselect_filter(const std::vector<unsigned int>& l) : l_(l) {}
      bool operator()(const GaussianShell& shell,
                    unsigned int contr) {
        return std::find(l_.begin(), l_.end(), shell.am(contr)) != l_.end();
      }
      const std::vector<unsigned int>& l_;
  };
}

void
LSelectBasisSet::lselect(const Ref<GaussianBasisSet>&basis)
{
  molecule_ = basis->molecule();

    // create shells
    std::vector<Shell> shells;
    lselect_filter f(l_);
    for (int s = 0; s < basis->nshell(); ++s) {
      const GaussianShell &shell = basis->shell(s);

      shells.push_back(Shell(this, basis->shell_to_center(s), filter(shell, f)));
    }

    init(name_conv(basis->name(),l_),
         name_conv(basis->label(),l_),
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

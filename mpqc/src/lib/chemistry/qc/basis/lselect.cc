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

#ifdef __GNUC__
#pragma implementation
#endif

#include <algorithm>
#include <util/state/statein.h>
#include <util/state/stateout.h>
#include <util/class/scexception.h>
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
  if (basis.null()) {
      throw InputError("missing basis value: must be GaussianBasisSet",
                       __FILE__, __LINE__,
                       "basis", 0, class_desc());
  }

  if (keyval->exists("l")) {
      int nl = keyval->count("l");
      if (nl < 1)
	  throw InputError("l must be an array of at least 1 element",__FILE__, __LINE__);
      for(int c=0; c<nl; ++c) {
	  l_.push_back(keyval->intvalue("l",c));
      }
  }
  else {
      int lmin = keyval->intvalue("lmin",KeyValValueint(0));
      if (lmin < 0) lmin = 0;
      int lmax = keyval->intvalue("lmax",KeyValValueint(basis->max_angular_momentum()));
      if (lmin > lmax)
	  throw InputError("lmin > lmax",__FILE__, __LINE__);
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
    char *
    name_conv(const char *name, const std::vector<int>& l_)
    {
	std::ostringstream oss;
	if (name == 0) return 0;
	oss << "LSelect(" << name << ":";
	for(std::vector<int>::const_iterator am=l_.begin(); am!=l_.end(); ++am)
	    oss << GaussianShell::amtypes[*am];
	oss << ")";
	const std::string& newname = oss.str();
	return strcpy(new char[newname.size()+1],newname.c_str());
    }


    class __shell_selector {
    public:
	__shell_selector(const std::vector<int>& l) : lbegin_(l.begin()), lend_(l.end()) {}

	// return true if this contraction is in l_
	bool selected(const GaussianShell& shell, int contraction) {
	    const int l=shell.am(contraction);
	    const std::vector<int>::const_iterator v = find(lbegin_,lend_,l);
	    if (v != lend_)
	      return true;
	    else
	      return false;
	}

	// return the number of contractions from shell that are found in l_
	int num_selected(const GaussianShell& shell) {
	    const int nc = shell.ncontraction();
	    int nselected = 0;
	    // check every angular momentum, if some are found in l_, return true
	    for(int c=0; c<nc; ++c) {
	      if (selected(shell,c)) ++nselected;
	    }
	    return nselected;
	}

	// return true if only some contractions from shell are found in l_
	bool select_some_only(const GaussianShell& shell) {
	    // only generally contracted shells can return true
	    const int nc = shell.ncontraction();
	    if (nc == 1) return false;
	    const int lmax = shell.max_angular_momentum();
	    const int lmin = shell.min_angular_momentum();
	    if (lmin == lmax) return false;
	    // check every angular momentum, if only some are found in l_, return true
	    bool found = false;
	    bool not_found = false;
	    for(int c=0; c<nc; ++c) {
		const int l=shell.am(c);
		const std::vector<int>::const_iterator v = find(lbegin_,lend_,l);
		if (v != lend_)
		    found = true;
		else
		    not_found = true;
		if (found == true && not_found == true)
		    return true;
	    }
	    return false;
	}
    private:
	const std::vector<int>::const_iterator lbegin_;
	const std::vector<int>::const_iterator lend_;
    };
}

void
LSelectBasisSet::lselect(const Ref<GaussianBasisSet>&basis)
{
    __shell_selector ss(l_);

    // compute the number of shells on each center
    std::vector<int> nshell_on_center(basis->ncenter(),0);
    int nshell = 0;
    for (int icenter=0; icenter<basis->ncenter(); icenter++) {
	for (int ishell=0; ishell<basis->nshell_on_center(icenter);
	     ++ishell) {
	    if (ss.num_selected(basis->shell(icenter,ishell)))
		nshell_on_center[icenter] += 1;
        }
	nshell += nshell_on_center[icenter];
    }

    // create shells
    GaussianShell **shells = new GaussianShell*[nshell];
    int ishell = 0;
    std::vector<int> shell_to_center(nshell);
    for (int icenter=0, ishellall=0; icenter<basis->ncenter(); icenter++) {
	for (int ishell=0; ishell<basis->nshell_on_center(icenter);
	     ishell++) {
	    const GaussianShell &shell = basis->shell(icenter,ishell);
	    int ncon = ss.num_selected(shell);
	    if (ncon == 0)
		continue;
	    int nprim = shell.nprimitive();
	    shell_to_center[ishellall] = icenter;

	    // storage
	    int *am = new int[ncon];
	    int *pure = new int[ncon];
	    double *exponents = new double[nprim*ncon];
	    double **c = new double*[ncon];
	    c[0] = new double[nprim*ncon];
	    for(int i=1; i<ncon; ++i)
		c[i] = c[i-1] + nprim;

	    // loop over all contractions of the original shell
	    const int ncon_orig = shell.ncontraction();
	    int icon = 0;
	    for (int icon_orig=0; icon_orig<ncon_orig; ++icon_orig) {
		if (!ss.selected(shell,icon_orig))
		    continue;
		am[icon] = shell.am(icon_orig);
		pure[icon] = shell.is_pure(icon_orig);
		for (int iprim=0; iprim<nprim; iprim++) {
		    exponents[iprim] = shell.exponent(iprim);
		    c[icon][iprim] = shell.coefficient_unnorm(icon_orig,iprim);
                }
		++icon;
            }

	    shells[ishellall]
		= new GaussianShell(ncon, nprim, exponents, am,
				    pure, c,
				    GaussianShell::Unnormalized);
	    ++ishellall;
        }
    }
    
    init(name_conv(basis->name(),l_),
	 name_conv(basis->label(),l_),
	 basis->molecule(),
	 basis->matrixkit(),
	 basis->so_matrixkit(),
	 shells,
	 shell_to_center);

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

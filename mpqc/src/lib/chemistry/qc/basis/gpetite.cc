//
// gpetite.cc --- implementation of GPetite4 and helpers
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@ca.sandia.gov>
// Maintainer: SNL
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

#include <util/misc/formio.h>

#include <chemistry/qc/basis/basis.h>
#include <chemistry/qc/basis/gpetite.h>

using namespace std;
using namespace sc;

////////////////////////////////////////////////////////////////////////////

canonical_aaaa::canonical_aaaa()
{
}

canonical_aaaa::canonical_aaaa(const Ref<GaussianBasisSet> bi,
                               const Ref<GaussianBasisSet> bj,
                               const Ref<GaussianBasisSet> bk,
                               const Ref<GaussianBasisSet> bl)
{
}

////////////////////////////////////////////////////////////////////////////

canonical_aabc::canonical_aabc(const Ref<GaussianBasisSet> bi,
                               const Ref<GaussianBasisSet> bj,
                               const Ref<GaussianBasisSet> bk,
                               const Ref<GaussianBasisSet> bl)
{
  nk_ = bk->nshell();
  nl_ = bl->nshell();
}

////////////////////////////////////////////////////////////////////////////


canonical_aabb::canonical_aabb(const Ref<GaussianBasisSet> bi,
                               const Ref<GaussianBasisSet> bj,
                               const Ref<GaussianBasisSet> bk,
                               const Ref<GaussianBasisSet> bl)
{
  int ni = bi->nshell();
  nij_ = (ni*long(ni+1))>>1;
}

////////////////////////////////////////////////////////////////////////////

canonical_abcd::canonical_abcd(const Ref<GaussianBasisSet> bi,
                               const Ref<GaussianBasisSet> bj,
                               const Ref<GaussianBasisSet> bk,
                               const Ref<GaussianBasisSet> bl)
{
  ni_ = bi->nshell();
  nj_ = bj->nshell();
  nk_ = bk->nshell();
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:

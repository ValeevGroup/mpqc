//
// linkage.h
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Edward Seidl <seidl@janed.com>
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

#ifndef _chemistry_qc_basis_linkage_h
#define _chemistry_qc_basis_linkage_h

#include <chemistry/qc/basis/uncontract.h>
#include <chemistry/qc/basis/split.h>
#include <chemistry/qc/basis/lselect.h>
#include <chemistry/qc/basis/union.h>
#include <chemistry/qc/basis/gaussbas.h>

namespace sc {

ForceLink<UncontractedBasisSet> basis_force_link_a_;
ForceLink<SplitBasisSet> basis_force_link_b_;
ForceLink<LSelectBasisSet> basis_force_link_c_;
ForceLink<UnionBasisSet> basis_force_link_d_;
ForceLink<WriteBasisGrid> basis_force_link_e_;

}

#endif

//
// linkage.h
//
// Copyright (C) 1997 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@ca.sandia.gov>
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

#ifndef _chemistry_qc_wfn_linkage_h
#define _chemistry_qc_wfn_linkage_h

#include <chemistry/qc/wfn/density.h>
#include <chemistry/qc/wfn/orbital.h>
#include <chemistry/qc/wfn/eht.h>
#include <chemistry/qc/wfn/esp.h>
#include <chemistry/qc/wfn/rdm.h>

namespace sc {

ForceLink<ElectronDensity> wfn_force_link_a_;
ForceLink<Orbital> wfn_force_link_b_;
ForceLink<ExtendedHuckelWfn> wfn_force_link_d_;
ForceLink<WriteElectrostaticPotential> wfn_force_link_e_;
ForceLink<HCoreWfn> wfn_force_link_f_;
ForceLink<OBWfnRDMOne> wfn_force_link_g_;
ForceLink<OBWfnRDMTwo> wfn_force_link_h_;

}

#endif

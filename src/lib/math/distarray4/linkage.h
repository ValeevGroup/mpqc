//
// linkage.h
//
// Copyright (C) 1996 Limit Point Systems, Inc.
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

#ifndef _math_distarray4_linkage_h
#define _math_distarray4_linkage_h

#include <math/distarray4/distarray4.h>
#include <math/distarray4/distarray4_node0file.h>
#include <math/distarray4/distarray4_memgrp.h>
#include <math/distarray4/distarray4_mpiio.h>

namespace sc {

ForceLink<DistArray4_MemoryGrp> math_distarray4_force_link_a_;
ForceLink<DistArray4_Node0File> math_distarray4_force_link_b_;
ForceLink<DistArray4_MPIIO>     math_distarray4_force_link_c_;

}

#endif

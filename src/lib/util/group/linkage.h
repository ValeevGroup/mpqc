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

#ifndef _util_group_linkage_h
#define _util_group_linkage_h

#ifdef HAVE_CONFIG_H
#include <mpqc_config.h>
#endif

#include <util/group/memproc.h>

namespace sc {
static ForceLink<ProcMessageGrp> group_force_link_0_;
static ForceLink<ProcMemoryGrp> group_force_link_1_;
}

# if defined(HAVE_PTHREAD)
#   include <util/group/thpthd.h>
namespace sc {
    static ForceLink<PthreadThreadGrp> group_force_link_c_;
}
# endif

#if defined(HAVE_MPI)
#   include <util/group/memmtmpi.h>
namespace sc {
    static ForceLink<MTMPIMemoryGrp> group_force_link_g_;
}
#endif

#if defined(HAVE_ARMCI)
#   include <util/group/memarmci.h>
namespace sc {
    static ForceLink<ARMCIMemoryGrp> group_force_link_h_;
}
#endif

#endif

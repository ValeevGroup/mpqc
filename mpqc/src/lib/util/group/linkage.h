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

#ifndef __PIC__

#ifdef HAVE_CONFIG_H
#include <scconfig.h>
#endif

static ForceLink<ProcMessageGrp> group_force_link_;

# ifdef HAVE_SYSV_IPC
#   include <util/group/messshm.h>
    static ForceLink<ShmMessageGrp> group_force_link_a_;
# endif

# ifdef HAVE_NX
#   include <util/group/messpgon.h>
#   include <util/group/memipgon.h>
    static ForceLink<ParagonMessageGrp> group_force_link_b_;
    static ForceLink<IParagonMemoryGrp> group_force_link_b1_;
# endif

# if defined(HAVE_PTHREAD)
#   include <util/group/thpthd.h>
    static ForceLink<PthreadThreadGrp> group_force_link_c_;
# endif

#if defined(HAVE_PUMA_MPI2)
#   include <util/group/thpuma.h>
    static ForceLink<PumaThreadGrp> group_force_link_d_;
#endif

#if defined(HAVE_MPI)
#   include <util/group/memmpi.h>
    static ForceLink<MPIMemoryGrp> group_force_link_e_;
#endif

#if defined(HAVE_MPI2_ONE_SIDED)
#   include <util/group/memmpi2.h>
    static ForceLink<MPI2MemoryGrp> group_force_link_f_;
#endif

#if defined(HAVE_MPI)
#   include <util/group/memmtmpi.h>
    static ForceLink<MTMPIMemoryGrp> group_force_link_g_;
#endif

#endif /* __PIC__ */


#endif

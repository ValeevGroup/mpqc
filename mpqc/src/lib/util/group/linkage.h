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

const ClassDesc &group_force_link_ = ProcMessageGrp::class_desc_;

# ifdef HAVE_SYSV_IPC
#   include <util/group/messshm.h>
    const ClassDesc &group_force_link_a_ = ShmMessageGrp::class_desc_;
# endif

# ifdef HAVE_NX
#   include <util/group/messpgon.h>
#   include <util/group/memipgon.h>
    const ClassDesc &group_force_link_b_ = ParagonMessageGrp::class_desc_;
    const ClassDesc &group_force_link_b1_ = IParagonMemoryGrp::class_desc_;
# endif

# if defined(HAVE_LIBPTHREAD) || defined(HAVE_LIBPTHREADS)
#   include <util/group/thpthd.h>
    const ClassDesc &group_force_link_c_ = PthreadThreadGrp::class_desc_;
# endif

#if defined(HAVE_PUMA_MPI2)
#   include <util/group/thpuma.h>
    const ClassDesc &group_force_link_d_ = PumaThreadGrp::class_desc_;
#endif

#if defined(HAVE_MPI)
#   include <util/group/memmpi.h>
    const ClassDesc &group_force_link_e_ = MPIMemoryGrp::class_desc_;
#endif

#endif /* __PIC__ */


#endif

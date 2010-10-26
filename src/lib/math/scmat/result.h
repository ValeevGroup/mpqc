//
// result.h
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

#ifndef _math_scmat_result_h
#define _math_scmat_result_h
#ifdef __GNUC__
#pragma interface
#endif

#include <util/misc/compute.h>
#include <math/scmat/matrix.h>

namespace sc {

typedef AccResult<RefSCMatrix > AccResultRefSCMatrix;
typedef AccResult<RefSymmSCMatrix > AccResultRefSymmSCMatrix;
typedef AccResult<RefDiagSCMatrix > AccResultRefDiagSCMatrix;
typedef AccResult<RefSCVector > AccResultRefSCVector;

typedef AccResult<RefSCMatrix > ResultRefSCMatrix;
typedef AccResult<RefSymmSCMatrix > ResultRefSymmSCMatrix;
typedef AccResult<RefDiagSCMatrix > ResultRefDiagSCMatrix;
typedef AccResult<RefSCVector > ResultRefSCVector;

}

#endif

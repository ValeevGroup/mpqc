//
// array.h
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

#ifdef __GNUC__
#pragma interface
#endif

#ifndef _util_container_array_h
#define _util_container_array_h

#include <iostream.h>
#include <iomanip.h>
#include <stdlib.h>

#include <util/misc/exenv.h>

// Simple array template and macro classes.
#include <util/container/artem.h>
#include <util/container/armac.h>
#define ARRAY_dec(Type) Array_declare(Type)
#define ARRAY_def(Type)

// These are simple arrays that have StateIn CTORS and save_object_state,
// but do not actually inherit from SavableState at this time.  They only
// work if Type is a basic type like int, double, etc.
// The SSB macros and templates cannot be used unless <util/state/state.h>
// is included.  This cannot be done here since state.h grabs this file.
class StateIn;
class StateOut;
#include <util/container/ssartem.h>
#include <util/container/ssarmac.h>
#define SSB_ARRAY_dec(Type) SSBArray_declare(Type)
#define SSB_ARRAY_def(Type)

// Two dimensional versions of the above class.
#include <util/container/ar2tem.h>
#include <util/container/ar2mac.h>
#define ARRAY2_dec(Type) Array2_declare(Type)
#define ARRAY2_def(Type)

#include <util/container/ssar2tem.h>
#include <util/container/ssar2mac.h>
#define SSB_ARRAY2_dec(Type) SSBArray2_declare(Type)
#define SSB_ARRAY2_def(Type)

// Declare arrays of the basic types.
// The SSB versions cannot be declared because state.h cannot be included
// yet.  Use the template versions instead.

ARRAY_dec(int);
ARRAY2_dec(int);

ARRAY_dec(Arrayint);

ARRAY_dec(double);
ARRAY2_dec(double);

ARRAY_dec(Arraydouble);

#endif

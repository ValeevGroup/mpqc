//
// set.h
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

#ifndef _util_container_set_h
#define _util_container_set_h

#include <Pix.h>
#include <stdlib.h>
#include <iomanip.h>
#include <util/container/array.h>

//#include <util/container/settmpl.h> // The template set declaration.
#include <util/container/gnuset.h>
#include <util/container/gnuavlse.h>
#include <util/container/setmacr.h> // The macro set declaration.
#define SET_dec(Type) Set_declare(Type)
#define SET_def(Type)

// This class combines array capabilities with set capabilities.
// It is basically a set with the iseek and operator[](int i) members
// added.  At the moment, it requires that Set use an array internally.
// When Set is improved, Arrayset must be updated.
//#include <util/container/asettmpl.h> // The template set declaration.
#include <util/container/avlaset.h>
#include <util/container/asetmacr.h> // The macro set declaration.
#define ARRAYSET_dec(Type) Arrayset_declare(Type)
#define ARRAYSET_def(Type)

// declare sets for the basic types
SET_dec(int);
ARRAYSET_dec(int);
ARRAY_dec(Arraysetint);
SET_dec(double);
ARRAYSET_dec(double);

#endif

//
// isotest.cc
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

extern "C" {
#include <math.h>
#ifdef SGI
#include <sigfpe.h>
#endif // SGI
}

#include <util/misc/ieee.h>
#include <util/misc/formio.h>

#include <util/keyval/keyval.h>
#include <math/isosurf/volume.h>
#include <math/isosurf/shape.h>
#include <math/isosurf/surf.h>
#include <math/isosurf/isosurf.h>

// Force linkages:
#ifndef __PIC__
const ClassDesc &fl0 = SphereShape::class_desc_;
const ClassDesc &fl1 = TriangulatedImplicitSurface::class_desc_;
#endif

main(int argc,char** argv)
{
  ieee_trap_errors();

  // The first argument is the input filename.  If it doesn't
  // exist or is a '.', then the default ioput file is used.
  const char *defaultinput = SRCDIR "/isotest.in";
  const char *input;
  if (argc == 1 || !strcmp(argv[1], ".")) input = defaultinput;
  else input = argv[1];

  const char *keyword = (argc > 2)? argv[2] : "surf";

  // open keyval input
  RefKeyVal rpkv(new ParsedKeyVal(input));

  RefTriangulatedSurface surf = rpkv->describedclassvalue(keyword);

  cout << scprintf("surf->flat_area() = %f\n", surf->flat_area());
  cout << scprintf("surf->flat_volume() = %f\n", surf->flat_volume());

  cout << scprintf("surf->area() = %f\n", surf->area());
  cout << scprintf("surf->volume() = %f\n", surf->volume());

  return 0;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// eval: (c-set-style "CLJ")
// End:

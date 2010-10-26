//
// isotest.cc
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@limitpt.com>
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

using namespace std;
using namespace sc;

// Force linkages:
#ifndef __PIC__
static ForceLink<SphereShape> fl0;
static ForceLink<TriangulatedImplicitSurface> fl1;
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
  Ref<KeyVal> rpkv(new ParsedKeyVal(input));

  Ref<TriangulatedSurface> surf;
  surf << rpkv->describedclassvalue(keyword);

  cout << scprintf("surf->flat_area() = %f\n", surf->flat_area());
  cout << scprintf("surf->flat_volume() = %f\n", surf->flat_volume());

  cout << scprintf("surf->area() = %f\n", surf->area());
  cout << scprintf("surf->volume() = %f\n", surf->volume());

  for (int order=1; order<4; order++) { 
      cout << scprintf("order = %d", order) << endl;
      int nir = 4;
      for (int ir=0; ir<nir; ir++) {
          double r = double(ir)/(nir-1);
          double s = 1.0 - r;
          cout << scprintf("  r = %6.4f, s = %6.4f", r, s) << endl;
          TriInterpCoefKey key(order, r, s);
          Ref<TriInterpCoef> coef = new TriInterpCoef(key);
          for (int i=0; i<=order; i++) {
              for (int j=0; j+i<=order; j++) {
                  int k = order - i - j;
                  cout << scprintf("    coef(%d,%d,%d) = %12.8f",
                                   i,j,k,coef->coef(i,j,k))
                       << endl;
                }
            }
        }
    }

  return 0;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:


extern "C" {
#include <stdio.h>
#include <math.h>
#ifdef SGI
#include <sigfpe.h>
#endif // SGI
}
#include <util/misc/ieee.h>

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

  printf("surf->flat_area() = %f\n", surf->flat_area());
  printf("surf->flat_volume() = %f\n", surf->flat_volume());

  printf("surf->area() = %f\n", surf->area());
  printf("surf->volume() = %f\n", surf->volume());

  return 0;
}

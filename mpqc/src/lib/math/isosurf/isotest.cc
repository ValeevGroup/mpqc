
extern "C" {
#include <stdio.h>
#include <math.h>
#ifdef SGI
#include <sigfpe.h>
#endif // SGI
}
#include <util/options/GetLongOpt.h>
#include <util/misc/ieee.h>

#include "volume.h"
#include "shape.h"
#include "surf.h"
#include "isosurf.h"

class TestTriangleIntegrator: public TriangleIntegrator {
  public:
    TestTriangleIntegrator(): TriangleIntegrator(1) {
        set_r(0,0.0);
        set_s(0,0.0);
        set_w(0,0.5);
      };
        
    ~TestTriangleIntegrator() {};
};

main(int argc,char** argv)
{
  ieee_trap_errors();

  GetLongOpt option;
  option.enroll("help",GetLongOpt::NoValue,
                "print this option summary",0);
  option.enroll("noclean",GetLongOpt::NoValue,
                "don't clean up the surface",0);
  option.enroll("simple",GetLongOpt::NoValue,
                "use simple triangles","");
  option.enroll("resolution", GetLongOpt::MandatoryValue,
                "set the resolution", "0.5");
  option.enroll("short",  GetLongOpt::MandatoryValue,
                "set the short edge length to resolution * this", "0.1");
  option.enroll("maxval", GetLongOpt::MandatoryValue,
                "set the maximum value used for the bounding box", "1.0");
  int optind = option.parse(argc,argv);
  if (optind < 1) return -1;
  if ( option.retrieve("help") ) {
      option.usage();
      return 0;
    }

  SCVector3 center;
  center[0] = center[1] = center[2] = 0.0;
  double radius = 1.0;
  RefVolume vol(new SphereShape(center,radius));
  double resolution = atof(option.retrieve("resolution"));
  double maxval = atof(option.retrieve("maxval"));

  ImplicitSurfacePolygonizer isogen(vol);
  isogen.set_resolution(resolution);

  double distance_from_surface = 0.5;
  TriangulatedSurface * surfptr;
  if (option.retrieve("simple")) {
      printf("using 10 point triangles\n");
      surfptr = new TriangulatedSurface10(vol,distance_from_surface);
    }
  else {
      surfptr = new TriangulatedSurface();
    }
  TriangulatedSurface& surf = *surfptr;
  isogen.isosurface(distance_from_surface,surf);

  surf.fix_orientation();
  if (!option.retrieve("noclean")) {
      double shortl = atof(option.retrieve("short")) * resolution;

      FILE* fp = fopen("isotestse.off","w");
      surf.print_geomview_format(fp);
      fclose(fp);
      
      surf.remove_short_edges(shortl);

      fp = fopen("isotestst.off","w");
      surf.print_geomview_format(fp);
      fclose(fp);

      surf.remove_slender_triangles(shortl);
    }
  surf.fix_orientation();

  FILE* fp = fopen("isotest.off","w");
  //surf.print_vertices_and_triangles(fp);
  //surf.remove_short_edges();
  //surf.print_vertices_and_triangles(fp);
  surf.print_geomview_format(fp);
  fclose(fp);

  printf("surface is written\n");

  //surf.print();
  printf("surf.area() = %f\n", surf.area());
  printf("surf.volume() = %f\n", surf.volume());

  double area = 0.0;
  for (int i=0; i<surf.ntriangle(); i++) {
      area += surf.triangle(i)->Triangle::area();
    }
  printf("flat triangle area = %f\n",area);

  // if surface is a sphere here are the correct value
  double sr = distance_from_surface + radius;
  double sa = 4.0 * M_PI * sr * sr;
  double sv = sa * sr / 3.0;
  printf("expected area = %f, expected volume = %f\n",sa,sv);

  delete surfptr;
  return 0;
}

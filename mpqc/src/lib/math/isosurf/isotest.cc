
extern "C" {
#include <stdio.h>
#include <math.h>
#ifdef SGI
#include <sigfpe.h>
#endif // SGI
}
#include <util/options/GetLongOpt.h>

#include "volume.h"
#include "shape.h"
//#include "tess.h"
#include "surf.h"

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
  GetLongOpt option;
  option.enroll("help",GetLongOpt::NoValue,
                "print this option summary",0);
  option.enroll("resolution", GetLongOpt::MandatoryValue,
                "set the resolution", "0.5");
  option.enroll("maxval", GetLongOpt::MandatoryValue,
                "set the maximum value used for the bounding box", "1.0");
  int optind = option.parse(argc,argv);
  if (optind < 1) return -1;
  if ( option.retrieve("help") ) {
      option.usage();
      return 0;
    }

#ifdef SGI
  sigfpe_[_OVERFL].abort=1;  
  sigfpe_[_DIVZERO].abort=1;  
  sigfpe_[_INVALID].abort=1;
  handle_sigfpes(_ON,
                 _EN_OVERFL | _EN_DIVZERO | _EN_INVALID,
                 0,
                 _ABORT_ON_ERROR,
                 0);
#endif
  SCVector3 center;
  center[0] = center[1] = center[2] = 0.0;
  double radius = 1.0;
  RefVolume vol(new SphereShape(center,radius));
  //double resolution = 0.5;
  double resolution = atof(option.retrieve("resolution"));
  double minval = 0.0;
  double maxval = atof(option.retrieve("maxval"));

// the tesselation stuff is producing strange results
//   ArraysetRefPoint points;
//   vol->pointset(resolution, minval, maxval, points);
// 
//   DirichletTesselation tess(points);
// 
//   printf("points:\n");
//   for (int i=0; i<points.length(); i++) {
//       fprintf(stderr," %d ",i);
//       points[i]->print(stderr);
//     }
// 
//   tess.print();

  ImplicitUniformLattice lat(vol,resolution,minval,maxval);
  MCubesIsosurfaceGen isogen(lat);

  double distance_from_surface = 0.5;
  //TriangulatedSurface10 surf(vol,distance_from_surface);
  TriangulatedSurface surf;
  surf.set_integrator(new GaussTriangleIntegrator(7));
  //surf.set_integrator(new TestTriangleIntegrator());
  isogen.isosurface(distance_from_surface,surf);

  FILE* fp = fopen("isotest.out","w");
  surf.print_vertices_and_triangles(fp);
  surf.remove_short_edges();
  surf.print_vertices_and_triangles(fp);
  fclose(fp);

  printf("surface is written\n");

  //surf.print();
  printf("surf.area() = %f, surf.volume() = %f\n",
         surf.area(),surf.volume());

  // if surface is a sphere here are the correct value
  double sr = distance_from_surface + radius;
  double sa = 4.0 * M_PI * sr * sr;
  double sv = sa * sr / 3.0;
  printf("expected area = %f, expected volume = %f\n",sa,sv);
}

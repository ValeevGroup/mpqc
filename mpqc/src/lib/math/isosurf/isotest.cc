
extern "C" {
#include <stdio.h>
#include <math.h>
#ifdef SGI
#include <sigfpe.h>
#endif // SGI
}
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

main()
{
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
  Point3 center;
  center[0] = center[1] = center[2] = 0.0;
  double radius = 1.0;
  RefVolume vol(new SphereShape(center,radius));
  //double resolution = 0.5;
  double resolution = 0.5;
  double minval = 0.0;
  double maxval = 1.0;

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
  TriangulatedSurface10 surf(vol,distance_from_surface);
  //TriangulatedSurface surf;
  surf.set_integrator(new GaussTriangleIntegrator(7));
  //surf.set_integrator(new TestTriangleIntegrator());
  isogen.isosurface(distance_from_surface,surf);

  //surf.print();
  printf("surf.area() = %f, surf.volume() = %f\n",
         surf.area(),surf.volume());

  // if surface is a sphere here are the correct value
  double sr = distance_from_surface + radius;
  double sa = 4.0 * M_PI * sr * sr;
  double sv = sa * sr / 3.0;
  printf("expected area = %f, expected volume = %f\n",sa,sv);

  FILE* fp = fopen("isotest.out","w");
  surf.print(fp);
  fclose(fp);
}


#include "tess.h"

extern "C" {
#include "p3dgen.h"
#include "ge_error.h"
#include "dirichlet.h"
}

// this is needed while the creation operator is running for
// DirichletTesselation
static int dim;

static float* access(float*coord,int i,P_Void_ptr*)
{
  return &coord[dim*i];
}

DirichletTesselation::DirichletTesselation(ArraysetRefPoint&points):
     Tesselation(points)
{
  if (points.length() == 0) return;
  // hopefully all points will have the same dimension
  dim = points[0]->dimension();
  float* data = new float[points.length()*dim];

  int ij = 0;
  for (int i=0; i<points.length(); i++) {
      for (int j=0; j<dim; j++) data[ij] = points[i]->operator[](j);
      ij++;
    }

  ger_init("DirichletTesselation");
  _tess = dch_create_dirichlet_tess(data, points.length(), dim, access);

  delete[] data;
}

DirichletTesselation::~DirichletTesselation()
{
  ger_init("DirichletTesselation");
  dch_destroy_tesselation(_tess);  
}

void DirichletTesselation::print()
{
  dch_dump_tesselation(_tess);
}


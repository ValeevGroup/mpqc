/****************************************************************************
 * isosurface.c
 * Author Joel Welling
 * Copyright 1991, Pittsburgh Supercomputing Center, Carnegie Mellon University
 *
 * Permission use, copy, and modify this software and its documentation
 * without fee for personal use or use within your organization is hereby
 * granted, provided that the above copyright notice is preserved in all
 * copies and that that copyright and this permission notice appear in
 * supporting documentation.  Permission to redistribute this software to
 * other organizations or individuals is not granted;  that must be
 * negotiated with the PSC.  Neither the PSC nor Carnegie Mellon
 * University make any representations about the suitability of this
 * software for any purpose.  It is provided "as is" without express or
 * implied warranty.
 *****************************************************************************/
/*
This module provides the capability to extract isosurfaces from 3D gridded
data using a Marching Cubes algorithm (see Lorensen and Cline, Siggraph '87
proceedings;  corrections are needed to the complements of their cases 3,
6, and 7.
*/

/*
 * Cube vertices are numbered as follows:
 *
 *         3                     2
 *        *---------------------*
 *        |\                    :\
 *        | \ 0                 : \ 1
 *        |  *---------------------*
 *        |7 |                  :6 |
 *        *..|..................*  |
 *         \ |                   . |
 *          \|4                   .|5
 *           *---------------------*
 *
 *
 *  All of this is done with reference to "Marching Cubes:  A High 
 *  Resolution 3D Surface Construction Algorithm" by William
 *  Lorensen and Harvey Cline, Computer Graphics Volume 21, number 4,
 *  July 1987, pp. 163-169.  This paper contains a number of errors
 *  discoverable by trying to mate their pattern #3 to its complement
 *  (what you get when you reverse high and low values).  See a letter
 *  in a subsequent Computer Graphics (1987 or 1988) for the correction
 *  and another reference.  The first correction necessitates changes
 *  to the patterns for the complements of their patterns 6 and 7.
 *
 *  The coordinate system within a cube is that the x direction is
 *  along the (0,1) edge, the y direction is along the (1,2) edge,
 *  and the z direction is along the (4,0) edge.  In the loops, the
 *  indices i, j, and k refer to the x, y, and z directions respectively.
 *  The loops run such that i increments fastest and k slowest.
 *
 *  This routine also uses a short array of 12 vertices to access the
 *  vertices needed for a single cube.  The mapping of array indices
 *  in that array to cube edges is as follows:
 *
 *  edge (0,1) : index 0
 *  edge (1,2) : index 1
 *  edge (2,3) : index 2
 *  edge (3,0) : index 3
 *  edge (0,4) : index 4
 *  edge (1,5) : index 5
 *  edge (2,6) : index 6
 *  edge (3,7) : index 7
 *  edge (4,5) : index 8
 *  edge (5,6) : index 9
 *  edge (6,7) : index 10
 *  edge (7,4) : index 11
 *
 *  There are a series of routines to calculate or look up needed vertices,
 *  specialized by where in the cube they fall and thus by which cubes have
 *  been calculated so far.  Their names and (i,j,k) values are:
 *
 *  calc_vertex_main                 i>0, j>0, k>0
 *  calc_vertex_first_cube           i=0, j=0, k=0
 *  calc_vertex_first_row            i>0, j=0, k=0
 *  calc_vertex_begin_first_face     i=0, j>0, k=0
 *  calc_vertex_first_face           i>0, j>0, k=0
 *  calc_vertex_begin_main           i=0, j>0, k>0
 *  calc_vertex_corner_main          i=0, j=0, k>0
 *  calc_vertex_main_first_row       i>0, j=0, k>0
 *
 *  It is assumed that the calling order proceeds through i, j, and k,
 *  with i cycling fastest and k cycling slowest.  Some of the routines
 *  know to call an appropriate special case routine if i or j is zero.
 */

/*
  Notes-
  -whichcase out of new_triangle?
*/

extern "C" {
#include <stdio.h>
#include <math.h>
#include "p3dgen.h"
#include "ge_error.h"
}
#include "surf.h"
#include "isosurf.h"

/* Structure from which to build list of vertices */
typedef struct P_Vertex_struct {
  int type;
  float x;
  float y;
  float z;
  P_Vector *normal;
  float value;
  struct P_Vertex_struct *next;
  int index;
} P_Vertex;

/* Structure to hold a P_Vector, for use in memory management */
typedef struct P_Vector_List_struct {
  P_Vector vec;
  struct P_Vector_List_struct *next;
} P_Vector_List;

/* The following data needs to get saved from cubelet to cubelet */
typedef struct cell_data_struct {
  P_Vertex *right;
  P_Vertex *back;
  int halfcase;
} cell_data;

/* Structure from which to build list of triangular facets */
typedef struct P3DTriangle_struct {
  P_Vertex *v1;
  P_Vertex *v2;
  P_Vertex *v3;
  struct P3DTriangle_struct *next;
} P3DTriangle;

/* Information about the isosurface as a whole */
static float contour_value;
static float deltax, deltay, deltaz;
static P_Point *corner1_save, *corner2_save;
static int flip_normals= 0, left_handed_coords= 0;
int nx, ny, nz;
int ftn_order_flag= 0; /* true if the left array index increments fastest */

/* Macros which access data and test for inside-ness */
#define ACCESS( grid, i, j, k ) \
  (ftn_order_flag ? grid[k][j][i] : grid[i][j][k])
#define IN_CHECK( i, j, k ) (ACCESS( grid, i, j, k ) >= contour_value)

/* handles for the data structures for the data grid and the
 * areas in which data from previous calculations are saved
 */
static float ***grid= (float ***)0;
static float ***valgrid= (float ***)0;
static cell_data **old_plane_saver= (cell_data **)0;
static cell_data **new_plane_saver= (cell_data **)0;
static P_Vertex **old_row_saver= (P_Vertex **)0;
static P_Vertex **new_row_saver= (P_Vertex **)0;

/* Information for the vertex and triangle lists */
static int current_type; /* vertex type */
static P_Vertex *free_vertex_list= (P_Vertex *)0;
static P_Vertex *vertex_list= (P_Vertex *)0;
static P_Vector_List *free_vector_list= (P_Vector_List *)0;
static P3DTriangle *triangle_list= (P3DTriangle *)0;
static P3DTriangle *free_triangle_list= (P3DTriangle *)0;
static int vertex_count= 0, triangle_count= 0;

/* Prototypes for vertex calculation functions */
static void calc_vertex_main( int ind1, int ind2, int i, int j, int k,
			     P_Vertex **varray );
static void calc_vertex_first_cube( int ind1, int ind2, P_Vertex **varray );
static void calc_vertex_first_row( int ind1, int ind2, int i, 
				  P_Vertex **varray );
static void calc_vertex_begin_first_face( int ind1, int ind2, int j,
				  P_Vertex **varray );
static void calc_vertex_first_face( int ind1, int ind2, int i, int j,
				  P_Vertex **varray );
static void calc_vertex_corner_main( int ind1, int ind2, int k,
				    P_Vertex **varray );
static void calc_vertex_main_first_row( int ind1, int ind2, int i, int k,
				  P_Vertex **varray );
static void calc_vertex_begin_main( int ind1, int ind2, int j, int k,
				  P_Vertex **varray );

static void vertex_space_setup( VOIDLIST )
/* This routine allows initialization of the vertex handling
 * routines based on the grid size.
 */
{
  static int initialized= 0;
  ger_debug("isosurf: vertex_space_setup: doing nothing");
}

static P_Vector_List *allocate_vectors( VOIDLIST )
/* This routine allocates a block of vectors */
{
  static int totalvecs= 0;
  int blocksize, i;
  P_Vector_List *result;

  /* If this is the first pass, guess at a number to allocate.  Otherwise,
   * allocate as many as currently exist.
   */
  if (!totalvecs) blocksize= nx*ny;
  else blocksize= totalvecs;

  ger_debug("isosurf: allocate_vectors: blocksize will be %d",blocksize);

  if ( !(result= (P_Vector_List *)malloc(blocksize*sizeof(P_Vector_List))) )
    ger_fatal("isosurf: allocate_vectors: unable to allocate %d bytes!",
	      blocksize*sizeof(P_Vector_List));

  for (i=0; i<blocksize-1; i++) result[i].next= &result[i+1];
  result[blocksize-1].next= (P_Vector_List *)0;

  totalvecs += blocksize;
  return result;
}

static P_Vector *new_vector(VOIDLIST)
{
  P_Vector *result;

  if (!free_vector_list) free_vector_list= allocate_vectors();

  result= &(free_vector_list->vec);
  free_vector_list= free_vector_list->next;

  result->x= result->y= result->z= 0.0;
  return result;
}

static void free_vector( P_Vector *thisvector )
{
  P_Vector_List *thiscell;

  /* Note that we rely on the fact that the address of thisvector is
   * the same as that of the P_Vector_List cell that holds it.
   */
  thiscell= (P_Vector_List *)thisvector;
  thiscell->next= free_vector_list;
  free_vector_list= thiscell;
}

static P_Vertex *allocate_vertices( VOIDLIST )
/* This routine allocates a block of vertices */
{
  static int totalverts= 0;
  int blocksize, i;
  P_Vertex *result;

  /* If this is the first pass, guess at a number to allocate.  Otherwise,
   * allocate as many as currently exist.
   */
  if (!totalverts) blocksize= nx*ny;
  else blocksize= totalverts;

  ger_debug("isosurf: allocate_vertices: blocksize will be %d",blocksize);

  if ( !(result= (P_Vertex *)malloc(blocksize*sizeof(P_Vertex))) )
    ger_fatal("isosurf: allocate_vertices: unable to allocate %d bytes!",
	      blocksize*sizeof(P_Vertex));

  for (i=0; i<blocksize-1; i++) result[i].next= &result[i+1];
  result[blocksize-1].next= (P_Vertex *)0;

  totalverts += blocksize;
  return result;
}

static P_Vertex *new_vertex(int type)
{
  P_Vertex *result;

  if (!free_vertex_list) free_vertex_list= allocate_vertices();

  result= free_vertex_list;
  free_vertex_list= result->next;

  result->type= type;
  result->next= vertex_list;
  vertex_list= result;
  result->x= result->y= result->z= 0.0;
  result->value= 0.0;
  result->index= -1;
  result->normal= (P_Vector *)0;
  vertex_count++;
  return result;
}

static void free_vertex( P_Vertex *thisvertex )
{
  if ( (thisvertex->type==P3D_CNVTX) || (thisvertex->type==P3D_CVNVTX)
      || (thisvertex->type==P3D_CCNVTX) )
    if (thisvertex->normal) free_vector(thisvertex->normal);

  thisvertex->next= free_vertex_list;
  free_vertex_list= thisvertex;
  vertex_count--;
}

static P3DTriangle *allocate_triangles( VOIDLIST )
/* This routine allocates a block of triangles */
{
  static int totaltris= 0;
  int blocksize, i;
  P3DTriangle *result;

  /* If this is the first pass, guess at a number to allocate.  Otherwise,
   * allocate as many as currently exist.
   */
  if (!totaltris) blocksize= nx*ny;
  else blocksize= totaltris;

  ger_debug("isosurf: allocate_triangles: blocksize will be %d",blocksize);

  if ( !(result= (P3DTriangle *)malloc(blocksize*sizeof(P3DTriangle))) )
    ger_fatal("isosurf: allocate_triangles: unable to allocate %d bytes!",
	      blocksize*sizeof(P3DTriangle));

  for (i=0; i<blocksize-1; i++) result[i].next= &result[i+1];
  result[blocksize-1].next= (P3DTriangle *)0;

  totaltris += blocksize;
  return result;
}

static void new_triangle( int whichcase,
			 P_Vertex *v1, P_Vertex *v2, P_Vertex *v3 )
{
  P3DTriangle *result;

  if ( !v1 || !v2 || !v3 )
    ger_fatal(
      "isosurf: new_triangle: algorithm error; null vertex for case %d!",
      whichcase, whichcase);

  if (!free_triangle_list) free_triangle_list= allocate_triangles();

  result= free_triangle_list;
  free_triangle_list= result->next;

  /* If we want the inner surface or the coordinates are left handed,
   * but not both, we need to permute the vertex order to maintain the
   * right hand rule for which face is outward.
   */
  if ( (flip_normals || left_handed_coords) && 
      !(flip_normals && left_handed_coords)) {
    result->v1= v2;
    result->v2= v1;
    result->v3= v3;
  }
  else {
    result->v1= v1;
    result->v2= v2;
    result->v3= v3;
  }

  result->next= triangle_list;
  triangle_list= result;
  triangle_count++;
}

static void free_triangle( P3DTriangle *thistriangle )
{
  thistriangle->next= free_triangle_list;
  free_triangle_list= thistriangle;
  triangle_count--;
}

static void cleanup( VOIDLIST )
/* This function frees and cleans up the vertex and triangle lists */
{
  P_Vertex *thisvtx, *nextvtx;
  P3DTriangle *thistri, *nexttri;

  ger_debug("isosurf: cleanup");

  thisvtx= vertex_list;
  while( thisvtx ) {
    nextvtx= thisvtx->next;
    free_vertex( thisvtx );
    thisvtx= nextvtx;
  }
  vertex_list= (P_Vertex *)0;
  vertex_count= 0;

  thistri= triangle_list;
  while( thistri ) {
    nexttri= thistri->next;
    free_triangle( thistri );
    thistri= nexttri;
  }
  triangle_list= (P3DTriangle *)0;
  triangle_count= 0;
}

static float deriv_forwards( float v1, float v2, float v3, float step )
/* This routine takes a first derivative to second order accuracy by
 * forward differencing.
 */
{
  return( (-0.5*v3 +2.0*v2 -1.5*v1)/step );
}

static float deriv_centered( float v1, float v2, float step )
/* This routine takes a first derivative to second order accuracy by
 * centered differencing.
 */
{
  return( (v2-v1)/(2.0*step) );
}

static void calc_gradient( float *gradx, float *grady, float *gradz,
			  int i, int j, int k )
/* This routine calculates the components of the gradient at the 
 * given point.
 */
{
    if (i==0) *gradx= 
      deriv_forwards( ACCESS(grid,0,j,k), ACCESS(grid,1,j,k), 
		     ACCESS(grid,2,j,k), deltax );
    else if (i==nx-1) *gradx=
      deriv_forwards( ACCESS(grid,nx-1,j,k), ACCESS(grid,nx-2,j,k), 
		     ACCESS(grid,nx-3,j,k), -deltax );
    else *gradx=
      deriv_centered( ACCESS(grid,i-1,j,k), ACCESS(grid,i+1,j,k), deltax );

    if (j==0) *grady= 
      deriv_forwards( ACCESS(grid,i,0,k), ACCESS(grid,i,1,k), 
		     ACCESS(grid,i,2,k), deltay );
    else if (j==ny-1) *grady=
      deriv_forwards( ACCESS(grid,i,ny-1,k), ACCESS(grid,i,ny-2,k), 
		     ACCESS(grid,i,ny-3,k), -deltay );
    else *grady=
      deriv_centered( ACCESS(grid,i,j-1,k), ACCESS(grid,i,j+1,k), deltay );

    if (k==0) *gradz= 
      deriv_forwards( ACCESS(grid,i,j,0), ACCESS(grid,i,j,1), 
		     ACCESS(grid,i,j,2), deltaz );
    else if (k==nz-1) *gradz=
      deriv_forwards( ACCESS(grid,i,j,nz-1), ACCESS(grid,i,j,nz-2), 
		     ACCESS(grid,i,j,nz-3), -deltaz );
    else *gradz=
      deriv_centered( ACCESS(grid,i,j,k-1), ACCESS(grid,i,j,k+1), deltaz );
}

static UniformLattice* lattice;

static P_Vertex *interp_vertex(int i1, int j1, int k1, int i2, int j2, int k2)
/* This routine generates a vertex by interpolation. */
{
  P_Vertex *vtx;

//   float fraction;
// 
//   vtx= new_vertex( current_type );
// 
//   /* Calculate interpolation factor */
//   fraction= (ACCESS(grid,i2,j2,k2) - contour_value) /
//     (ACCESS(grid,i2,j2,k2)-ACCESS(grid,i1,j1,k1));
// 
//   /* Calculate vertex coordinates */
//   vtx->x= corner1_save->x + i1*deltax + (1.0-fraction)*(i2-i1)*deltax;
//   vtx->y= corner1_save->y + j1*deltay + (1.0-fraction)*(j2-j1)*deltay;
//   vtx->z= corner1_save->z + k1*deltaz + (1.0-fraction)*(k2-k1)*deltaz;

  vtx= new_vertex( current_type );

  RefPoint p = lattice->interpolate(i1,j1,k1,i2,j2,k2,contour_value);

  vtx->x = p->operator[](0);
  vtx->y = p->operator[](1);
  vtx->z = p->operator[](2);

  return vtx;
}

static void calc_vertex_main( int ind1, int ind2, int i, int j, int k,
			     P_Vertex **varray )
/* This function calculates a needed vertex for the main body of the
 * data cube, and puts it in the appropriate place to be grabbed
 * by the vertex access functions.  It is assumed that it will be
 * called in order of increasing i, j, and k, with i increasing
 * fastest and k slowest.
 */
{
  if (k==0) calc_vertex_first_face( ind1, ind2, i, j, varray );
  else if (i==0) calc_vertex_begin_main( ind1, ind2, j, k, varray );
  else if (j==0) calc_vertex_main_first_row( ind1, ind2, i, k, varray );
  else switch (ind1) {
  case 0: 
    switch (ind2) {
    case 1: 
      varray[0]= new_plane_saver[i][j].right;  /* saved on last j */
      break;
    case 3: 
      varray[3]= new_plane_saver[i][j].back; /* saved on last i */
      break;
    case 4:
      varray[4]= old_row_saver[i]; /* saved on last j */
      break;
    default: ger_fatal("Forbidden cubes case: vertex on edge %d-%d",
		       ind1,ind2);
      break;
    }
    break;
  case 1: 
    switch (ind2) {
    case 0: 
      varray[0]= new_plane_saver[i][j].right;  /* saved on last j */
      break;
    case 2: 
      varray[1]= new_plane_saver[i+1][j].back= 
	interp_vertex( i+1, j, k+1, i+1, j+1, k+1 );
      break;
    case 5: 
      varray[5]= old_row_saver[i+1]; /* saved on last j */
      break;
    default: ger_fatal("Forbidden cubes case: vertex on edge %d-%d",
		       ind1,ind2);
      break;
    }
    break;
  case 2: 
    switch (ind2) {
    case 1: 
      varray[1]= new_plane_saver[i+1][j].back= 
	interp_vertex( i+1, j, k+1, i+1, j+1, k+1 );
      break;
    case 3: 
      varray[2]= new_plane_saver[i][j+1].right=
	interp_vertex( i, j+1, k+1, i+1, j+1, k+1 );
      break;
    case 6: 
      varray[6]= new_row_saver[i+1]=
	interp_vertex( i+1, j+1, k+1, i+1, j+1, k );
      break;
    default: ger_fatal("Forbidden cubes case: vertex on edge %d-%d",
		       ind1,ind2);
      break;
    }
    break;
  case 3: 
    switch (ind2) {
    case 0: 
      varray[3]= new_plane_saver[i][j].back; /* saved on last i */
      break;
    case 2: 
      varray[2]= new_plane_saver[i][j+1].right=
	interp_vertex( i, j+1, k+1, i+1, j+1, k+1 );
      break;
    case 7: 
      varray[7]= new_row_saver[i]; /* saved last i */
      break;
    default: ger_fatal("Forbidden cubes case: vertex on edge %d-%d",
		       ind1,ind2);
      break;
    }
    break;
  case 4: 
    switch (ind2) {
    case 0: 
      varray[4]= old_row_saver[i]; /* saved on last j */
      break;
    case 5: 
      varray[8]= old_plane_saver[i][j].right; /* saved last k */
      break;
    case 7: 
      varray[11]= old_plane_saver[i][j].back; /* saved last k */
      break;
    default: ger_fatal("Forbidden cubes case: vertex on edge %d-%d",
		       ind1,ind2);
      break;
    }
    break;
  case 5: 
    switch (ind2) {
    case 1: 
      varray[5]= old_row_saver[i+1]; /* saved on last j */
      break;
    case 4: 
      varray[8]= old_plane_saver[i][j].right; /* saved last k */
      break;
    case 6: 
      varray[9]= old_plane_saver[i+1][j].back; /* saved last k */
      break;
    default: ger_fatal("Forbidden cubes case: vertex on edge %d-%d",
		       ind1,ind2);
      break;
    }
    break;
  case 6: 
    switch (ind2) {
    case 2: 
      varray[6]= new_row_saver[i+1]=
	interp_vertex( i+1, j+1, k+1, i+1, j+1, k );
      break;
    case 5: 
      varray[9]= old_plane_saver[i+1][j].back; /* saved last k */
      break;
    case 7: 
      varray[10]= old_plane_saver[i][j+1].right; /* saved last k */
      break;
    default: ger_fatal("Forbidden cubes case: vertex on edge %d-%d",
		       ind1,ind2);
      break;
    }
    break;
  case 7: 
    switch (ind2) {
    case 3: 
      varray[7]= new_row_saver[i]; /* saved on last i */
      break;
    case 4: 
      varray[11]= old_plane_saver[i][j].back; /* saved last k */
      break;
    case 6: 
      varray[10]= old_plane_saver[i][j+1].right; /* saved last k */
      break;
    default: ger_fatal("Forbidden cubes case: vertex on edge %d-%d",
		       ind1,ind2);
      break;
    }
    break;
  }
}

static void calc_vertex_first_cube( int ind1, int ind2, P_Vertex **varray )
/* This function calculates a needed vertex for the first cubelet of the
 * data cube, and puts it in the appropriate place to be grabbed
 * by the vertex access functions.  This is the first cube handled, so
 * i, j, and k are assumed to be 0.
 */
{
  switch (ind1) {
  case 0: 
    switch (ind2) {
    case 1: 
      varray[0]= new_plane_saver[0][0].right=
	interp_vertex( 0, 0, 1, 1, 0, 1 );
      break;
    case 3: 
      varray[3]= new_plane_saver[0][0].back=
	interp_vertex( 0, 0, 1, 0, 1, 1 );
      break;
    case 4:
      varray[4]=
	interp_vertex( 0, 0, 1, 0, 0, 0 );
      break;
    default: ger_fatal("Forbidden cubes case: vertex on edge %d-%d",ind1,ind2);
      break;
    }
    break;
  case 1: 
    switch (ind2) {
    case 0: 
      varray[0]= new_plane_saver[0][0].right=
	interp_vertex( 0, 0, 1, 1, 0, 1 );
      break;
    case 2: 
      varray[1]= new_plane_saver[1][0].back= 
	interp_vertex( 1, 0, 1, 1, 1, 1 );
      break;
    case 5: 
      varray[5]= old_row_saver[1]=
	interp_vertex( 1, 0, 1, 1, 0, 0 );
      break;
    default: ger_fatal("Forbidden cubes case: vertex on edge %d-%d",ind1,ind2);
      break;
    }
    break;
  case 2: 
    switch (ind2) {
    case 1: 
      varray[1]= new_plane_saver[1][0].back= 
	interp_vertex( 1, 0, 1, 1, 1, 1 );
      break;
    case 3: 
      varray[2]= new_plane_saver[0][1].right=
	interp_vertex( 0, 1, 1, 1, 1, 1 );
      break;
    case 6: 
      varray[6]= new_row_saver[1]=
	interp_vertex( 1, 1, 1, 1, 1, 0 );
      break;
    default: ger_fatal("Forbidden cubes case: vertex on edge %d-%d",ind1,ind2);
      break;
    }
    break;
  case 3: 
    switch (ind2) {
    case 0: 
      varray[3]= new_plane_saver[0][0].back=
	interp_vertex( 0, 0, 1, 0, 1, 1 );
      break;
    case 2: 
      varray[2]= new_plane_saver[0][1].right=
	interp_vertex( 0, 1, 1, 1, 1, 1 );
      break;
    case 7: 
      varray[7]= new_row_saver[0]=
	interp_vertex( 0, 1, 1, 0, 1, 0 );
      break;
    default: ger_fatal("Forbidden cubes case: vertex on edge %d-%d",ind1,ind2);
      break;
    }
    break;
  case 4: 
    switch (ind2) {
    case 0: 
      varray[4]=
	interp_vertex( 0, 0, 1, 0, 0, 0 );
      break;
    case 5: 
      varray[8]=
	interp_vertex( 0, 0, 0, 1, 0, 0 );
      break;
    case 7: 
      varray[11]= 
	interp_vertex( 0, 0, 0, 0, 1, 0 );
      break;
    default: ger_fatal("Forbidden cubes case: vertex on edge %d-%d",ind1,ind2);
      break;
    }
    break;
  case 5: 
    switch (ind2) {
    case 1: 
      varray[5]= old_row_saver[1]=
	interp_vertex( 1, 0, 1, 1, 0, 0 );
      break;
    case 4: 
      varray[8]=
	interp_vertex( 0, 0, 0, 1, 0, 0 );
      break;
    case 6: 
      varray[9]= old_plane_saver[1][0].back=
	interp_vertex( 1, 0, 0, 1, 1, 0 );
      break;
    default: ger_fatal("Forbidden cubes case: vertex on edge %d-%d",ind1,ind2);
      break;
    }
    break;
  case 6: 
    switch (ind2) {
    case 2: 
      varray[6]= new_row_saver[1]=
	interp_vertex( 1, 1, 1, 1, 1, 0 );
      break;
    case 5: 
      varray[9]= old_plane_saver[1][0].back=
	interp_vertex( 1, 0, 0, 1, 1, 0 );
      break;
    case 7: 
      varray[10]= old_plane_saver[0][1].right=
	interp_vertex( 1, 1, 0, 0, 1, 0 );
      break;
    default: ger_fatal("Forbidden cubes case: vertex on edge %d-%d",ind1,ind2);
      break;
    }
    break;
  case 7: 
    switch (ind2) {
    case 3: 
      varray[7]= new_row_saver[0]=
	interp_vertex( 0, 1, 1, 0, 1, 0 );
      break;
    case 4: 
      varray[11]=
	interp_vertex( 0, 0, 0, 0, 1, 0 );
      break;
    case 6: 
      varray[10]= old_plane_saver[0][1].right=
	interp_vertex( 1, 1, 0, 0, 1, 0 );
      break;
    default: ger_fatal("Forbidden cubes case: vertex on edge %d-%d",ind1,ind2);
      break;
    }
    break;
  }
}

static void calc_vertex_first_row( int ind1, int ind2, int i, 
				  P_Vertex **varray )
/* This function calculates a needed vertex for the first row of cubelets,
 * and puts it in the appropriate place to be grabbed by the vertex
 * access functions.  It is assumed that calc_vertex_first_cube
 * has been called, and that j and k are 0.  The calling routine
 * guarantees that i>0 here.
 */
{
  switch (ind1) {
  case 0: 
    switch (ind2) {
    case 1: 
      varray[0]= new_plane_saver[i][0].right=
	interp_vertex( i, 0, 1, i+1, 0, 1 );
      break;
    case 3: 
      varray[3]= new_plane_saver[i][0].back; /* saved last i */
      break;
    case 4:
      varray[4]= old_row_saver[i]; /* saved last i */
      break;
    default: ger_fatal("Forbidden cubes case: vertex on edge %d-%d",ind1,ind2);
      break;
    }
    break;
  case 1: 
    switch (ind2) {
    case 0: 
      varray[0]= new_plane_saver[i][0].right=
	interp_vertex( i, 0, 1, i+1, 0, 1 );
      break;
    case 2: 
      varray[1]= new_plane_saver[i+1][0].back= 
	interp_vertex( i+1, 0, 1, i+1, 1, 1 );
      break;
    case 5: 
      varray[5]= old_row_saver[i+1]=
	interp_vertex( i+1, 0, 1, i+1, 0, 0 );
      break;
    default: ger_fatal("Forbidden cubes case: vertex on edge %d-%d",ind1,ind2);
      break;
    }
    break;
  case 2: 
    switch (ind2) {
    case 1: 
      varray[1]= new_plane_saver[i+1][0].back= 
	interp_vertex( i+1, 0, 1, i+1, 1, 1 );
      break;
    case 3: 
      varray[2]= new_plane_saver[i][1].right=
	interp_vertex( i, 1, 1, i+1, 1, 1 );
      break;
    case 6: 
      varray[6]= new_row_saver[i+1]=
	interp_vertex( i+1, 1, 1, i+1, 1, 0 );
      break;
    default: ger_fatal("Forbidden cubes case: vertex on edge %d-%d",ind1,ind2);
      break;
    }
    break;
  case 3: 
    switch (ind2) {
    case 0: 
      varray[3]= new_plane_saver[i][0].back; /* saved last i */
      break;
    case 2: 
      varray[2]= new_plane_saver[i][1].right=
	interp_vertex( i, 1, 1, i+1, 1, 1 );
      break;
    case 7: 
      varray[7]= new_row_saver[i]; /* saved last i */
      break;
    default: ger_fatal("Forbidden cubes case: vertex on edge %d-%d",ind1,ind2);
      break;
    }
    break;
  case 4: 
    switch (ind2) {
    case 0: 
      varray[4]= old_row_saver[i]; /* saved last i */
      break;
    case 5: 
      varray[8]=
	interp_vertex( i, 0, 0, i+1, 0, 0 );
      break;
    case 7: 
      varray[11]= old_plane_saver[i][0].back; /* saved last i */
      break;
    default: ger_fatal("Forbidden cubes case: vertex on edge %d-%d",ind1,ind2);
      break;
    }
    break;
  case 5: 
    switch (ind2) {
    case 1: 
      varray[5]= old_row_saver[i+1]=
	interp_vertex( i+1, 0, 1, i+1, 0, 0 );
      break;
    case 4: 
      varray[8]=
	interp_vertex( 0, 0, 0, 1, 0, 0 );
      break;
    case 6: 
      varray[9]= old_plane_saver[i+1][0].back=
	interp_vertex( i+1, 0, 0, i+1, 1, 0 );
      break;
    default: ger_fatal("Forbidden cubes case: vertex on edge %d-%d",ind1,ind2);
      break;
    }
    break;
  case 6: 
    switch (ind2) {
    case 2: 
      varray[6]= new_row_saver[i+1]=
	interp_vertex( i+1, 1, 1, i+1, 1, 0 );
      break;
    case 5: 
      varray[9]= old_plane_saver[i+1][0].back=
	interp_vertex( i+1, 0, 0, i+1, 1, 0 );
      break;
    case 7: 
      varray[10]= old_plane_saver[i][1].right= 
	interp_vertex( i+1, 1, 0, i, 1, 0 );
      break;
    default: ger_fatal("Forbidden cubes case: vertex on edge %d-%d",ind1,ind2);
      break;
    }
    break;
  case 7: 
    switch (ind2) {
    case 3: 
      varray[7]= new_row_saver[i]; /* saved last i */
      break;
    case 4: 
      varray[11]= old_plane_saver[i][0].back; /* saved last i */
      break;
    case 6: 
      varray[10]= old_plane_saver[i][1].right=
	interp_vertex( i+1, 1, 0, i, 1, 0 );
      break;
    default: ger_fatal("Forbidden cubes case: vertex on edge %d-%d",ind1,ind2);
      break;
    }
    break;
  }
}

static void calc_vertex_begin_first_face( int ind1, int ind2, int j,
				  P_Vertex **varray )
/* This function calculates a needed vertex for the first cell of the 
 * first face of cubelets, and puts it in the appropriate place to be 
 * grabbed by the vertex access functions.  It is assumed that i and
 * k are zero.
 */
{
  if (j==0) calc_vertex_first_cube( ind1, ind2, varray );
  else switch (ind1) {
  case 0: 
    switch (ind2) {
    case 1: 
      varray[0]= new_plane_saver[0][j].right; /* saved last j */
      break;
    case 3: 
      varray[3]= new_plane_saver[0][j].back=
	interp_vertex( 0, j, 1, 0, j+1, 1 );
      break;
    case 4:
      varray[4]= old_row_saver[0]; /* saved last j */
      break;
    default: ger_fatal("Forbidden cubes case: vertex on edge %d-%d",ind1,ind2);
      break;
    }
    break;
  case 1: 
    switch (ind2) {
    case 0: 
      varray[0]= new_plane_saver[0][j].right; /* saved last j */
      break;
    case 2: 
      varray[1]= new_plane_saver[1][j].back= 
	interp_vertex( 1, j, 1, 1, j+1, 1 );
      break;
    case 5: 
      varray[5]= old_row_saver[1]; /* saved last j */
      break;
    default: ger_fatal("Forbidden cubes case: vertex on edge %d-%d",ind1,ind2);
      break;
    }
    break;
  case 2: 
    switch (ind2) {
    case 1: 
      varray[1]= new_plane_saver[1][j].back= 
	interp_vertex( 1, j, 1, 1, j+1, 1 );
      break;
    case 3: 
      varray[2]= new_plane_saver[0][j+1].right=
	interp_vertex( 0, j+1, 1, 1, j+1, 1 );
      break;
    case 6: 
      varray[6]= new_row_saver[1]=
	interp_vertex( 1, j+1, 1, 1, j+1, 0 );
      break;
    default: ger_fatal("Forbidden cubes case: vertex on edge %d-%d",ind1,ind2);
      break;
    }
    break;
  case 3: 
    switch (ind2) {
    case 0: 
      varray[3]= new_plane_saver[0][j].back=
	interp_vertex( 0, j, 1, 0, j+1, 1 );
      break;
    case 2: 
      varray[2]= new_plane_saver[0][j+1].right=
	interp_vertex( 0, j+1, 1, 1, j+1, 1 );
      break;
    case 7: 
      varray[7]= new_row_saver[0]=
	interp_vertex( 0, j+1, 1, 0, j+1, 0 );
      break;
    default: ger_fatal("Forbidden cubes case: vertex on edge %d-%d",ind1,ind2);
      break;
    }
    break;
  case 4: 
    switch (ind2) {
    case 0: 
      varray[4]= old_row_saver[0]; /* saved last j */
      break;
    case 5: 
      varray[8]= old_plane_saver[0][j].right; /* saved last j */
      break;
    case 7: 
      varray[11]= 
	interp_vertex( 0, j+1, 0, 0, j, 0 );
      break;
    default: ger_fatal("Forbidden cubes case: vertex on edge %d-%d",ind1,ind2);
      break;
    }
    break;
  case 5: 
    switch (ind2) {
    case 1: 
      varray[5]= old_row_saver[1]; /* saved last j */
      break;
    case 4: 
      varray[8]= old_plane_saver[0][j].right; /* saved last j */
      break;
    case 6: 
      varray[9]= old_plane_saver[1][j].back=
	interp_vertex( 1, j, 0, 1, j+1, 0 );
      break;
    default: ger_fatal("Forbidden cubes case: vertex on edge %d-%d",ind1,ind2);
      break;
    }
    break;
  case 6: 
    switch (ind2) {
    case 2: 
      varray[6]= new_row_saver[1]=
	interp_vertex( 1, j+1, 1, 1, j+1, 0 );
      break;
    case 5: 
      varray[9]= old_plane_saver[1][j].back=
	interp_vertex( 1, j, 0, 1, j+1, 0 );
      break;
    case 7: 
      varray[10]= old_plane_saver[0][j+1].right=
	interp_vertex( 1, j+1, 0, 0, j+1, 0 );
      break;
    default: ger_fatal("Forbidden cubes case: vertex on edge %d-%d",ind1,ind2);
      break;
    }
    break;
  case 7: 
    switch (ind2) {
    case 3: 
      varray[7]= new_row_saver[0]=
	interp_vertex( 0, j+1, 1, 0, j+1, 0 );
      break;
    case 4: 
      varray[11]=
	interp_vertex( 0, j+1, 0, 0, j, 0 );
      break;
    case 6: 
      varray[10]= old_plane_saver[0][j+1].right=
	interp_vertex( 1, j+1, 0, 0, j+1, 0 );
      break;
    default: ger_fatal("Forbidden cubes case: vertex on edge %d-%d",ind1,ind2);
      break;
    }
    break;
  }
}

static void calc_vertex_first_face( int ind1, int ind2, int i, int j,
				  P_Vertex **varray )
/* This function calculates a needed vertex for the first face of cubelets,
 * and puts it in the appropriate place to be grabbed by the vertex
 * access functions.  It is assumed that k is zero, and that i increments
 * faster than j.
 */
{
  if (i==0) calc_vertex_begin_first_face( ind1, ind2, j, varray );
  else if (j==0) calc_vertex_first_row( ind1, ind2, i, varray );
  else switch (ind1) {
  case 0: 
    switch (ind2) {
    case 1: 
      varray[0]= new_plane_saver[i][j].right; /* saved last j */
      break;
    case 3: 
      varray[3]= new_plane_saver[i][j].back; /* saved last i */
      break;
    case 4:
      varray[4]= old_row_saver[i]; /* saved last j */
      break;
    default: ger_fatal("Forbidden cubes case: vertex on edge %d-%d",ind1,ind2);
      break;
    }
    break;
  case 1: 
    switch (ind2) {
    case 0: 
      varray[0]= new_plane_saver[i][j].right; /* saved last j */
      break;
    case 2: 
      varray[1]= new_plane_saver[i+1][j].back= 
	interp_vertex( i+1, j, 1, i+1, j+1, 1 );
      break;
    case 5: 
      varray[5]= old_row_saver[i+1]; /* saved last j */
      break;
    default: ger_fatal("Forbidden cubes case: vertex on edge %d-%d",ind1,ind2);
      break;
    }
    break;
  case 2: 
    switch (ind2) {
    case 1: 
      varray[1]= new_plane_saver[i+1][j].back= 
	interp_vertex( i+1, j, 1, i+1, j+1, 1 );
      break;
    case 3: 
      varray[2]= new_plane_saver[i][j+1].right=
	interp_vertex( i, j+1, 1, i+1, j+1, 1 );
      break;
    case 6: 
      varray[6]= new_row_saver[i+1]=
	interp_vertex( i+1, j+1, 1, i+1, j+1, 0 );
      break;
    default: ger_fatal("Forbidden cubes case: vertex on edge %d-%d",ind1,ind2);
      break;
    }
    break;
  case 3: 
    switch (ind2) {
    case 0: 
      varray[3]= new_plane_saver[i][j].back; /* saved last i */
      break;
    case 2: 
      varray[2]= new_plane_saver[i][j+1].right=
	interp_vertex( i, j+1, 1, i+1, j+1, 1 );
      break;
    case 7: 
      varray[7]= new_row_saver[i]; /* saved last i */
      break;
    default: ger_fatal("Forbidden cubes case: vertex on edge %d-%d",ind1,ind2);
      break;
    }
    break;
  case 4: 
    switch (ind2) {
    case 0: 
      varray[4]= old_row_saver[i]; /* saved last j */
      break;
    case 5: 
      varray[8]= old_plane_saver[i][j].right; /* saved last j */
      break;
    case 7: 
      varray[11]= old_plane_saver[i][j].back; /* saved last i */
      break;
    default: ger_fatal("Forbidden cubes case: vertex on edge %d-%d",ind1,ind2);
      break;
    }
    break;
  case 5: 
    switch (ind2) {
    case 1: 
      varray[5]= old_row_saver[i+1]; /* saved last j */
      break;
    case 4: 
      varray[8]= old_plane_saver[i][j].right; /* saved last j */
      break;
    case 6: 
      varray[9]= old_plane_saver[i+1][j].back=
	interp_vertex( i+1, j, 0, i+1, j+1, 0 );
      break;
    default: ger_fatal("Forbidden cubes case: vertex on edge %d-%d",ind1,ind2);
      break;
    }
    break;
  case 6: 
    switch (ind2) {
    case 2: 
      varray[6]= new_row_saver[i+1]=
	interp_vertex( i+1, j+1, 1, i+1, j+1, 0 );
      break;
    case 5: 
      varray[9]= old_plane_saver[i+1][j].back=
	interp_vertex( i+1, j, 0, i+1, j+1, 0 );
      break;
    case 7: 
      varray[10]= old_plane_saver[i][j+1].right=
	interp_vertex( i+1, j+1, 0, i, j+1, 0 );
      break;
    default: ger_fatal("Forbidden cubes case: vertex on edge %d-%d",ind1,ind2);
      break;
    }
    break;
  case 7: 
    switch (ind2) {
    case 3: 
      varray[7]= new_row_saver[i]; /* saved last i */
      break;
    case 4: 
      varray[11]= old_plane_saver[i][j].back; /* saved last i */
      break;
    case 6: 
      varray[10]= old_plane_saver[i][j+1].right=
	interp_vertex( i+1, j+1, 0, i, j+1, 0 );
      break;
    default: ger_fatal("Forbidden cubes case: vertex on edge %d-%d",ind1,ind2);
      break;
    }
    break;
  }
}

static void calc_vertex_corner_main( int ind1, int ind2, int k,
				    P_Vertex **varray )
/* This function calculates a needed vertex for the first cubelet of 
 * the a new plane, and puts it in the appropriate place to be grabbed
 * by the vertex access functions.  It is assumed that i and j are 0,
 * and that the last plane has just been finished.
 */
{
  switch (ind1) {
  case 0: 
    switch (ind2) {
    case 1: 
      varray[0]= new_plane_saver[0][0].right=
	interp_vertex( 0, 0, k+1, 1, 0, k+1 );
      break;
    case 3: 
      varray[3]= new_plane_saver[0][0].back=
	interp_vertex( 0, 0, k+1, 0, 1, k+1 );
      break;
    case 4:
      varray[4]=
	interp_vertex( 0, 0, k+1, 0, 0, k );
      break;
    default: ger_fatal("Forbidden cubes case: vertex on edge %d-%d",ind1,ind2);
      break;
    }
    break;
  case 1: 
    switch (ind2) {
    case 0: 
      varray[0]= new_plane_saver[0][0].right=
	interp_vertex( 0, 0, k+1, 1, 0, k+1 );
      break;
    case 2: 
      varray[1]= new_plane_saver[1][0].back= 
	interp_vertex( 1, 0, k+1, 1, 1, k+1 );
      break;
    case 5: 
      varray[5]= old_row_saver[1]=
	interp_vertex( 1, 0, k+1, 1, 0, k );
      break;
    default: ger_fatal("Forbidden cubes case: vertex on edge %d-%d",ind1,ind2);
      break;
    }
    break;
  case 2: 
    switch (ind2) {
    case 1: 
      varray[1]= new_plane_saver[1][0].back= 
	interp_vertex( 1, 0, k+1, 1, 1, k+1 );
      break;
    case 3: 
      varray[2]= new_plane_saver[0][1].right=
	interp_vertex( 0, 1, k+1, 1, 1, k+1 );
      break;
    case 6: 
      varray[6]= new_row_saver[1]=
	interp_vertex( 1, 1, k+1, 1, 1, k );
      break;
    default: ger_fatal("Forbidden cubes case: vertex on edge %d-%d",ind1,ind2);
      break;
    }
    break;
  case 3: 
    switch (ind2) {
    case 0: 
      varray[3]= new_plane_saver[0][0].back=
	interp_vertex( 0, 0, k+1, 0, 1, k+1 );
      break;
    case 2: 
      varray[2]= new_plane_saver[0][1].right=
	interp_vertex( 0, 1, k+1, 1, 1, k+1 );
      break;
    case 7: 
      varray[7]= new_row_saver[0]=
	interp_vertex( 0, 1, k+1, 0, 1, k );
      break;
    default: ger_fatal("Forbidden cubes case: vertex on edge %d-%d",ind1,ind2);
      break;
    }
    break;
  case 4: 
    switch (ind2) {
    case 0: 
      varray[4]=
	interp_vertex( 0, 0, k+1, 0, 0, k );
      break;
    case 5: 
      varray[8]= old_plane_saver[0][0].right; /* saved last k */
      break;
    case 7: 
      varray[11]= old_plane_saver[0][0].back; /* saved last k */
      break;
    default: ger_fatal("Forbidden cubes case: vertex on edge %d-%d",ind1,ind2);
      break;
    }
    break;
  case 5: 
    switch (ind2) {
    case 1: 
      varray[5]= old_row_saver[1]=
	interp_vertex( 1, 0, k+1, 1, 0, k );
      break;
    case 4: 
      varray[8]= old_plane_saver[0][0].right; /* saved last k */
      break;
    case 6: 
      varray[9]= old_plane_saver[1][0].back; /* saved last k */
      break;
    default: ger_fatal("Forbidden cubes case: vertex on edge %d-%d",ind1,ind2);
      break;
    }
    break;
  case 6: 
    switch (ind2) {
    case 2: 
      varray[6]= new_row_saver[1]=
	interp_vertex( 1, 1, k+1, 1, 1, k );
      break;
    case 5: 
      varray[9]= old_plane_saver[1][0].back; /* saved last k */
      break;
    case 7: 
      varray[10]= old_plane_saver[0][1].right; /* saved last k */
      break;
    default: ger_fatal("Forbidden cubes case: vertex on edge %d-%d",ind1,ind2);
      break;
    }
    break;
  case 7: 
    switch (ind2) {
    case 3: 
      varray[7]= new_row_saver[0]=
	interp_vertex( 0, 1, k+1, 0, 1, k );
      break;
    case 4: 
      varray[11]= old_plane_saver[0][0].back; /* saved last k */
      break;
    case 6: 
      varray[10]= old_plane_saver[0][1].right; /* saved last k */
      break;
    default: ger_fatal("Forbidden cubes case: vertex on edge %d-%d",ind1,ind2);
      break;
    }
    break;
  }
}

static void calc_vertex_main_first_row( int ind1, int ind2, int i, int k,
				  P_Vertex **varray )
/* This function calculates a needed vertex for the first row of cubelets
 * of a new plane (after the first plane), and puts it in the appropriate 
 * place to be grabbed by the vertex access functions.  It is assumed that 
 * j is 0.  The calling routine guarantees i>0.
 */
{
  switch (ind1) {
  case 0: 
    switch (ind2) {
    case 1: 
      varray[0]= new_plane_saver[i][0].right=
	interp_vertex( i, 0, k+1, i+1, 0, k+1 );
      break;
    case 3: 
      varray[3]= new_plane_saver[i][0].back; /* saved last i */
      break;
    case 4:
      varray[4]= old_row_saver[i]; /* saved last i */
      break;
    default: ger_fatal("Forbidden cubes case: vertex on edge %d-%d",ind1,ind2);
      break;
    }
    break;
  case 1: 
    switch (ind2) {
    case 0: 
      varray[0]= new_plane_saver[i][0].right=
	interp_vertex( i, 0, k+1, i+1, 0, k+1 );
      break;
    case 2: 
      varray[1]= new_plane_saver[i+1][0].back= 
	interp_vertex( i+1, 0, k+1, i+1, 1, k+1 );
      break;
    case 5: 
      varray[5]= old_row_saver[i+1]=
	interp_vertex( i+1, 0, k+1, i+1, 0, k );
      break;
    default: ger_fatal("Forbidden cubes case: vertex on edge %d-%d",ind1,ind2);
      break;
    }
    break;
  case 2: 
    switch (ind2) {
    case 1: 
      varray[1]= new_plane_saver[i+1][0].back= 
	interp_vertex( i+1, 0, k+1, i+1, 1, k+1 );
      break;
    case 3: 
      varray[2]= new_plane_saver[i][1].right=
	interp_vertex( i, 1, k+1, i+1, 1, k+1 );
      break;
    case 6: 
      varray[6]= new_row_saver[i+1]=
	interp_vertex( i+1, 1, k+1, i+1, 1, k );
      break;
    default: ger_fatal("Forbidden cubes case: vertex on edge %d-%d",ind1,ind2);
      break;
    }
    break;
  case 3: 
    switch (ind2) {
    case 0: 
      varray[3]= new_plane_saver[i][0].back; /* saved last i */
      break;
    case 2: 
      varray[2]= new_plane_saver[i][1].right=
	interp_vertex( i, 1, k+1, i+1, 1, k+1 );
      break;
    case 7: 
      varray[7]= new_row_saver[i]; /* saved last i */
      break;
    default: ger_fatal("Forbidden cubes case: vertex on edge %d-%d",ind1,ind2);
      break;
    }
    break;
  case 4: 
    switch (ind2) {
    case 0: 
      varray[4]= old_row_saver[i]; /* saved last i */
      break;
    case 5: 
      varray[8]= old_plane_saver[i][0].right; /* saved last k */
      break;
    case 7: 
      varray[11]= old_plane_saver[i][0].back; /* saved last k */
      break;
    default: ger_fatal("Forbidden cubes case: vertex on edge %d-%d",ind1,ind2);
      break;
    }
    break;
  case 5: 
    switch (ind2) {
    case 1: 
      varray[5]= old_row_saver[i+1]=
	interp_vertex( i+1, 0, k+1, i+1, 0, k );
      break;
    case 4: 
      varray[8]= old_plane_saver[i][0].right; /* saved last k */
      break;
    case 6: 
      varray[9]= old_plane_saver[i+1][0].back; /* saved last k */
      break;
    default: ger_fatal("Forbidden cubes case: vertex on edge %d-%d",ind1,ind2);
      break;
    }
    break;
  case 6: 
    switch (ind2) {
    case 2: 
      varray[6]= new_row_saver[i+1]=
	interp_vertex( i+1, 1, k+1, i+1, 1, k );
      break;
    case 5: 
      varray[9]= old_plane_saver[i+1][0].back; /* saved last k */
      break;
    case 7: 
      varray[10]= old_plane_saver[i][1].right; /* saved last k */
      break;
    default: ger_fatal("Forbidden cubes case: vertex on edge %d-%d",ind1,ind2);
      break;
    }
    break;
  case 7: 
    switch (ind2) {
    case 3: 
      varray[7]= new_row_saver[i]; /* saved last i */
      break;
    case 4: 
      varray[11]= old_plane_saver[i][0].back; /* saved last k */
      break;
    case 6: 
      varray[10]= old_plane_saver[i][1].right; /* saved last k */
      break;
    default: ger_fatal("Forbidden cubes case: vertex on edge %d-%d",ind1,ind2);
      break;
    }
    break;
  }
}

static void calc_vertex_begin_main( int ind1, int ind2, int j, int k,
				  P_Vertex **varray )
/* This function calculates a needed vertex for the first cell of a 
 * plane after the first plane, and puts it in the appropriate place to be 
 * grabbed by the vertex access functions.  It is assumed that 
 * i is zero, and that previous planes are finished.
 */
{
  if (j==0) calc_vertex_corner_main( ind1, ind2, k, varray );
  else switch (ind1) {
  case 0: 
    switch (ind2) {
    case 1: 
      varray[0]= new_plane_saver[0][j].right; /* saved last j */
      break;
    case 3: 
      varray[3]= new_plane_saver[0][j].back=
	interp_vertex( 0, j, k+1, 0, j+1, k+1 );
      break;
    case 4:
      varray[4]= old_row_saver[0]; /* saved last j */
      break;
    default: ger_fatal("Forbidden cubes case: vertex on edge %d-%d",ind1,ind2);
      break;
    }
    break;
  case 1: 
    switch (ind2) {
    case 0: 
      varray[0]= new_plane_saver[0][j].right; /* saved last j */
      break;
    case 2: 
      varray[1]= new_plane_saver[1][j].back= 
	interp_vertex( 1, j, k+1, 1, j+1, k+1 );
      break;
    case 5: 
      varray[5]= old_row_saver[1]; /* saved last j */
      break;
    default: ger_fatal("Forbidden cubes case: vertex on edge %d-%d",ind1,ind2);
      break;
    }
    break;
  case 2: 
    switch (ind2) {
    case 1: 
      varray[1]= new_plane_saver[1][j].back= 
	interp_vertex( 1, j, k+1, 1, j+1, k+1 );
      break;
    case 3: 
      varray[2]= new_plane_saver[0][j+1].right=
	interp_vertex( 0, j+1, k+1, 1, j+1, k+1 );
      break;
    case 6: 
      varray[6]= new_row_saver[1]=
	interp_vertex( 1, j+1, k+1, 1, j+1, k );
      break;
    default: ger_fatal("Forbidden cubes case: vertex on edge %d-%d",ind1,ind2);
      break;
    }
    break;
  case 3: 
    switch (ind2) {
    case 0: 
      varray[3]= new_plane_saver[0][j].back=
	interp_vertex( 0, j, k+1, 0, j+1, k+1 );
      break;
    case 2: 
      varray[2]= new_plane_saver[0][j+1].right=
	interp_vertex( 0, j+1, k+1, 1, j+1, k+1 );
      break;
    case 7: 
      varray[7]= new_row_saver[0]=
	interp_vertex( 0, j+1, k+1, 0, j+1, k );
      break;
    default: ger_fatal("Forbidden cubes case: vertex on edge %d-%d",ind1,ind2);
      break;
    }
    break;
  case 4: 
    switch (ind2) {
    case 0: 
      varray[4]= old_row_saver[0]; /* saved last j */
      break;
    case 5: 
      varray[8]= old_plane_saver[0][j].right; /* saved last k */
      break;
    case 7: 
      varray[11]= old_plane_saver[0][j].back; /* saved last k */
      break;
    default: ger_fatal("Forbidden cubes case: vertex on edge %d-%d",ind1,ind2);
      break;
    }
    break;
  case 5: 
    switch (ind2) {
    case 1: 
      varray[5]= old_row_saver[1]; /* saved last j */
      break;
    case 4: 
      varray[8]= old_plane_saver[0][j].right; /* saved last k */
      break;
    case 6: 
      varray[9]= old_plane_saver[1][j].back; /* saved last k */
      break;
    default: ger_fatal("Forbidden cubes case: vertex on edge %d-%d",ind1,ind2);
      break;
    }
    break;
  case 6: 
    switch (ind2) {
    case 2: 
      varray[6]= new_row_saver[1]=
	interp_vertex( 1, j+1, k+1, 1, j+1, k );
      break;
    case 5: 
      varray[9]= old_plane_saver[1][j].back; /* saved last k */
      break;
    case 7: 
      varray[10]= old_plane_saver[0][j+1].right; /* saved last k */
      break;
    default: ger_fatal("Forbidden cubes case: vertex on edge %d-%d",ind1,ind2);
      break;
    }
    break;
  case 7: 
    switch (ind2) {
    case 3: 
      varray[7]= new_row_saver[0]=
	interp_vertex( 0, j+1, k+1, 0, j+1, k );
      break;
    case 4: 
    case 7: 
      varray[11]= old_plane_saver[0][j].back; /* saved last k */
      break;
    case 6: 
      varray[10]= old_plane_saver[0][j+1].right; /* saved last k */
      break;
    default: ger_fatal("Forbidden cubes case: vertex on edge %d-%d",ind1,ind2);
      break;
    }
    break;
  }
}

static P_Vertex *get_vertex( int ind1, int ind2, P_Vertex **varray )
/* This function returns needed vertices for the main body of the
 * data cube, computing them if necessary.  It is assumed that all
 * the necessary calc_vertex calls have been made for the current cube.
 */
{
  switch (ind1) {
  case 0: 
    switch (ind2) {
    case 1: return varray[0]; break;
    case 3: return varray[3]; break;
    case 4: return varray[4]; break;
    default: ger_fatal("Forbidden cubes case: vertex on edge %d-%d",ind1,ind2);
      break;
    }
    break;
  case 1: 
    switch (ind2) {
    case 0: return varray[0]; break;
    case 2: return varray[1]; break;
    case 5: return varray[5]; break;
    default: ger_fatal("Forbidden cubes case: vertex on edge %d-%d",ind1,ind2);
      break;
    }
    break;
  case 2:
    switch (ind2) {
    case 1: return varray[1]; break;
    case 3: return varray[2]; break;
    case 6: return varray[6]; break;
    default: ger_fatal("Forbidden cubes case: vertex on edge %d-%d",ind1,ind2);
      break;
    }
    break;
  case 3: 
    switch (ind2) {
    case 0: return varray[3]; break;
    case 2: return varray[2]; break;
    case 7: return varray[7]; break;
    default: ger_fatal("Forbidden cubes case: vertex on edge %d-%d",ind1,ind2);
      break;
    }
    break;
  case 4: 
    switch (ind2) {
    case 0: return varray[4]; break;
    case 5: return varray[8]; break;
    case 7: return varray[11]; break;
    default: ger_fatal("Forbidden cubes case: vertex on edge %d-%d",ind1,ind2);
      break;
    }
    break;
  case 5: 
    switch (ind2) {
    case 1: return varray[5]; break;
    case 4: return varray[8]; break;
    case 6: return varray[9]; break;
    default: ger_fatal("Forbidden cubes case: vertex on edge %d-%d",ind1,ind2);
      break;
    }
    break;
  case 6: 
    switch (ind2) {
    case 2: return varray[6]; break;
    case 5: return varray[9]; break;
    case 7: return varray[10]; break;
    default: ger_fatal("Forbidden cubes case: vertex on edge %d-%d",ind1,ind2);
      break;
    }
    break;
  case 7: 
    switch (ind2) {
    case 3: return varray[7]; break;
    case 4: return varray[11]; break;
    case 6: return varray[10]; break;
    default: ger_fatal("Forbidden cubes case: vertex on edge %d-%d",ind1,ind2);
      break;
    }
    break;
  }
  return 0; // avoid compiler warnings
}

static void do_general_plane( int k )
/* This routine calculates isosurface vertices and triangles 
 * for any plane k>0. It is assumed that the k's increase in order,
 * and that do_first_plane was called first.
 */
{
  int whichcase;
  P_Vertex *varray[12]; /* to hold edge vertices */
  P_Vertex **vholder;
  int i, j;

  /* The following definitions apply to the code included in cube_cases.c */
#define LIVE_CELL /* do nothing */
#define DEAD_CELL /* do nothing */
#define END_LIVE_CELL /* do nothing */
#define VTX( i1, i2 ) get_vertex( i1, i2, varray )
#define TRIANGLE( v1, v2, v3 ) new_triangle( whichcase, v1, v2, v3 )
#define NEEDS_VTX( i1, i2 ) calc_vertex_main( i1, i2, i, j, k, varray )

  for (j=0; j<ny-1; j++) {

    for (i=0; i<nx-1; i++) {
      new_plane_saver[i][j].halfcase= 
	(( IN_CHECK( i,j,k+1 ) ? 1 : 0 ) << 3)
	  | (( IN_CHECK( i+1,j,k+1 ) ? 1 : 0 ) << 2)
	    | (( IN_CHECK( i+1,j+1,k+1 ) ? 1 : 0 ) << 1)
	      | ( IN_CHECK( i,j+1,k+1 ) ? 1 : 0 );

      if (k==0) {
	old_plane_saver[i][j].halfcase= 
	  (( IN_CHECK( i,j,k ) ? 1 : 0 ) << 3)
	    | (( IN_CHECK( i+1,j,k ) ? 1 : 0 ) << 2)
	      | (( IN_CHECK( i+1,j+1,k ) ? 1 : 0 ) << 1)
		| ( IN_CHECK( i,j+1,k ) ? 1 : 0 );
      }

      whichcase= (new_plane_saver[i][j].halfcase << 4) 
	| old_plane_saver[i][j].halfcase;

      /* In this block we actually do the cube, by including the file 
       * containing the cube cases.
       */
      switch (whichcase) {
#include "cube_cases.c"
      }

    }

    /* swap grid row data spaces */
    vholder= new_row_saver;
    new_row_saver= old_row_saver;
    old_row_saver= vholder;
  }

  /* Clean up definitions */
#undef NEEDS_VTX
#undef LIVE_CELL
#undef DEAD_CELL
#undef END_LIVE_CELL
#undef VTX
#undef TRIANGLE
}

static void calc_isosurface(VOIDLIST)
{
  cell_data **holder;
  int k;

  ger_debug("isosurf: calc_isosurface:");

  for (k=0; k<nz-1; k++) {
    /* do a plane */
    do_general_plane(k);

    /* swap grid plane data spaces */
    holder= new_plane_saver;
    new_plane_saver= old_plane_saver;
    old_plane_saver= holder;
  }
}

static P_Void_ptr *create_indexed_2d_array( int nx, int ny, int cellsize )
/* This routine creates a 2 dimensional array, and array of index pointers
 * into it, so that it can be accessed like result[i][j].
 */
{
  P_Void_ptr *result;
  P_Void_ptr storage;
  int i;

  if ( !(storage= (P_Void_ptr)malloc(nx*ny*cellsize)) )
    ger_fatal("isosurf: create_indexed_2d_array: unable to allocate %d bytes!",
	      nx*ny*cellsize);

  if ( !(result= (P_Void_ptr *)malloc(nx*sizeof(P_Void_ptr *))) )
    ger_fatal("isosurf: create_indexed_2d_array: unable to allocate %d bytes!",
	      nx*cellsize);

  for (i=0; i<nx; i++) result[i]= (char *)storage+i*ny*cellsize;

  return( result );
}

static P_Void_ptr **index_3d_array( P_Void_ptr data, int nx, int ny, int nz,
				 int cellsize )
/* This routine creates index arrays for an existing 3 dimensional data
 * grid, so that it can be accessed like result[i][j][k] by successive
 * pointer dereferencing.
 */
{
  P_Void_ptr **result, **runner1;
  P_Void_ptr *block1, *runner2;
  P_Void_ptr runner3;
  int i;

  if (ftn_order_flag) {
    if ( !(block1= (P_Void_ptr *)malloc(nz*ny*sizeof(P_Void_ptr))) )
      ger_fatal("isosurf: index_3d_array: unable to allocate %d bytes!",
		nz*ny*sizeof(P_Void_ptr));
    
    if ( !(result= (P_Void_ptr **)malloc(nz*sizeof(P_Void_ptr *))) )
      ger_fatal("isosurf: index_3d_array: unable to allocate %d bytes!",
		nz*sizeof(P_Void_ptr *));
    
    runner1= result;
    runner2= block1;
    for (i=0; i<nz; i++) {
      *runner1++= runner2;
      runner2 += ny;
    }
    
    runner2= block1;
    runner3= data;
    for (i=0; i<nz*ny; i++) {
      *runner2++= runner3;
      runner3 = (char *)runner3 + nx*cellsize;
    }
  }
  else {
    if ( !(block1= (P_Void_ptr *)malloc(nx*ny*sizeof(P_Void_ptr))) )
      ger_fatal("isosurf: index_3d_array: unable to allocate %d bytes!",
		nx*ny*sizeof(P_Void_ptr));
    
    if ( !(result= (P_Void_ptr **)malloc(nx*sizeof(P_Void_ptr *))) )
      ger_fatal("isosurf: index_3d_array: unable to allocate %d bytes!",
		nx*sizeof(P_Void_ptr *));
    
    runner1= result;
    runner2= block1;
    for (i=0; i<nx; i++) {
      *runner1++= runner2;
      runner2 += ny;
    }
    
    runner2= block1;
    runner3= data;
    for (i=0; i<nx*ny; i++) {
      *runner2++= runner3;
      runner3 = (char *)runner3 + nz*cellsize;
    }
  }

  return result;
}

static void init_storage( float *data, float *valdata )
/* This routine makes sure that the proper static global data structures
 * exist.  It is smart enough to recreate them only when necessary.
 */
{
  static int nx_last= 0, ny_last= 0, nz_last= 0;
  static float *data_last= 0, *valdata_last= 0;

  if ( nx_last != nx ) {
    if (old_row_saver) free( (P_Void_ptr)old_row_saver );
    if ( !(old_row_saver= (P_Vertex **)malloc( nx*sizeof(P_Vertex *))) )
      ger_fatal("isosurf: init_storage: cannot allocate %d pointers!", nx);
    if (new_row_saver) free( (P_Void_ptr)new_row_saver );
    if ( !(new_row_saver= (P_Vertex **)malloc( nx*sizeof(P_Vertex *))) )
      ger_fatal("isosurf: init_storage: cannot allocate %d pointers!", nx);
  }

  if ( (nx_last != nx) || (ny_last != ny) ) {
    if (old_plane_saver) {
      free( (P_Void_ptr)(old_plane_saver[0]) );
      free( (P_Void_ptr)old_plane_saver );
    }
    old_plane_saver= 
      (cell_data **)create_indexed_2d_array( nx, ny, sizeof(cell_data) );
    if (new_plane_saver) {
      free( (P_Void_ptr)(new_plane_saver[0]) );
      free( (P_Void_ptr)new_plane_saver );
    }
    new_plane_saver= 
      (cell_data **)create_indexed_2d_array( nx, ny, sizeof(cell_data) );
  }

  /* Allocate data structures if previously allocated ones won't do */
  if ( (data_last != data) 
      || (nx_last != nx) || (ny_last != ny) || (nz_last != nz) ) {
    if (grid) {
      free( (P_Void_ptr)grid[0] );
      free( (P_Void_ptr)grid );
    }
    grid= (float ***)index_3d_array( (P_Void_ptr)data, nx, ny, nz, 
				    (int)sizeof(float) );

  }
  if ( (current_type==P3D_CVVTX) || (current_type==P3D_CVNVTX) )
    if ( (valdata_last != valdata) 
	|| (nx_last != nx) || (ny_last != ny) || (nz_last != nz) ) {
      if (valgrid) {
	free( (P_Void_ptr)valgrid[0] );
	free( (P_Void_ptr)valgrid );
      }
      valgrid= (float ***)index_3d_array( (P_Void_ptr)valdata, nx, ny, nz, 
					 (int)sizeof(float) );
      valdata_last= valdata;
    }

  nx_last= nx;
  ny_last= ny;
  nz_last= nz;
  data_last= data;
}

static RefEdge find_edge(ArrayArraysetint& edge_point_numbers,
                         ArrayArraysetRefEdge& edges,
                         ArrayRefVertex& points,
                         int i1,int i2)
{
  // canonicalize the indices so i can find the edge
  if (i1 > i2) {
      int tmp=i1; i1=i2; i2=tmp;
    }
  // find the index of the edge
  int index = edge_point_numbers[i1].iseek(i2);
  // if the index isn't valid, then edge didn't exist a new one is needed
  if (index < 0) {
      edge_point_numbers[i1].add(i2);
      edges[i1].add(new Edge(points[i1],points[i2]));
      index = edge_point_numbers[i1].iseek(i2);
    }
  return edges[i1][index];
}

static int convert_to_surface(TriangulatedSurface& surf,UniformLattice&lat)
/* This routine converts the generated vertex and triangle data to
 * a TriangulatedSurface object.
 */
{
  P3DTriangle *thistri;
  P_Vertex *thisvtx;
  int cell_length;
  int i;
  int retcode;

  ger_debug("isosurf: convert_to_p3d:");

  switch (current_type) {
  case P3D_CVTX: cell_length= 3; break;
  case P3D_CNVTX: cell_length= 6; break;
  case P3D_CVVTX: cell_length= 4; break;
  case P3D_CVNVTX: cell_length= 7; break;
  default: 
    ger_fatal("isosurf:convert_to_p3d: algorithm error; unknown type %d!");
  }

  ArrayRefVertex points(vertex_count);
  int nedge = (triangle_count/2)*3;
  ArrayArraysetint edge_point_numbers(vertex_count);
  ArrayArraysetRefEdge edges(vertex_count);

  /* Initialize the surface. */
  surf.initialize_vertices_triangles(vertex_count,triangle_count);

  /* Copy vertex info into the vertex buffer, adding vertex indices
   * as we go.
   */
  thisvtx= vertex_list;
  i= 0;
  DVector gradient(3);
  while (thisvtx) {
      RefPoint newpoint = new Point(thisvtx->x,thisvtx->y,thisvtx->z);
      lat.gradient(newpoint,gradient);
      points[i] = new Vertex(newpoint,gradient);
      surf.add_vertex(points[i]);
      thisvtx->index= i++;
      thisvtx= thisvtx->next;
  }

  /* Put the triangle info into the surface.
   */
  thistri= triangle_list;
  while (thistri) {
      RefEdge e1 = find_edge(edge_point_numbers,edges,points,
                             thistri->v1->index,thistri->v2->index);
      RefEdge e2 = find_edge(edge_point_numbers,edges,points,
                             thistri->v2->index,thistri->v3->index);
      RefEdge e3 = find_edge(edge_point_numbers,edges,points,
                             thistri->v1->index,thistri->v3->index);
      surf.add_triangle(new Triangle(e1,e2,e3),
                        thistri->v1->index,
                        thistri->v2->index,
                        thistri->v3->index);
    thistri= thistri->next;
  }

  return P3D_SUCCESS;
}

int pg_isosurface_ts( int type, float *data, float *valdata, 
                      int nx_in, int ny_in, int nz_in, 
                      double value, P_Point *corner1, P_Point *corner2, 
                      int show_inside, int ftn_order,
                      TriangulatedSurface& surface,
                      UniformLattice& lat)
{
  int retcode;

  ger_debug("pg_isosurface: nx= %d, ny= %d, nz= %d", nx_in, ny_in, nz_in);

  lattice = &lat;

  current_type= type;
  contour_value= value;
  nx= nx_in;
  ny= ny_in;
  nz= nz_in;
  deltax= (corner2->x - corner1->x)/(nx-1);
  deltay= (corner2->y - corner1->y)/(ny-1);
  deltaz= (corner2->z - corner1->z)/(nz-1);
  corner1_save= corner1;
  corner2_save= corner2;
  ftn_order_flag= ftn_order;

  /* Valid input checks */

  if (nx<2) {
    ger_error("pg_isosurface: nx must be at least 2; call ignored.");
    return(P3D_FAILURE);
  }

  if (ny<2) {
    ger_error("pg_isosurface: ny must be at least 2; call ignored.");
    return(P3D_FAILURE);
  }

  if (nz<2) {
    ger_error("pg_isosurface: nz must be at least 2; call ignored.");
    return(P3D_FAILURE);
  }

  if ((corner1->x == corner2->x)
      || (corner1->y == corner2->y)
      || (corner1->z == corner2->z) ) {
    ger_error("pg_isosurface: corners are coplanar; call ignored.");
    return(P3D_FAILURE);
  }

  if (type==P3D_CCVTX) {
    ger_error("pg_isosurface: CCVTX type invalid; using CVTX");
    type= P3D_CVTX;
  }

  if (type==P3D_CCNVTX) {
    ger_error("pg_isosurface: CCNVTX type invalid; using CNVTX");
    type= P3D_CNVTX;
  }

  if (type==P3D_CVVVTX) {
    ger_error("pg_isosurface: CVVVTX type invalid; using CVVTX");
    type= P3D_CVVTX;
  }

  /* Flip normals if requested, and note if the coordinate system
   * is left handed.
   */
  if (show_inside) flip_normals= 1;
  else flip_normals= 0;
  if (deltax*deltay*deltaz < 0) left_handed_coords= 1;
  else left_handed_coords= 0;

  init_storage( data, valdata );
  vertex_space_setup();

  /* Generate the isosurface.  The vertex and facet lists are stored on
   * static global pointers.
   */
  calc_isosurface();

  /* Add the isosurface to the TriangulatedSurface */
  retcode= convert_to_surface(surface,lat);

  /* Restore all the static global data structures to a reasonable
   * approximation of their initial states.
   */
  cleanup();

  return( retcode );
}

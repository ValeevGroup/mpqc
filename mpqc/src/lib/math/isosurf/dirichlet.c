/****************************************************************************
 * dirichlet.c
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

/* This module is an implementation of A. Bowyer's algorithm to calculate
 * the Dirichlet tesselation of a set of randomly distributed data points.
 * See The Computer Journal, Vol. 24, No. 2, 1981, pp. 162-166 for the
 * algorithm.  The algorithm runs in O( k**(1+1/n) ) time, where k is
 * the number of data points and n is the dimensionality of space.  The
 * process used is to create a valid tesselation from contrived data
 * (the 'bogus points' mentioned below), and then add additional points
 * in such a way that at every step the tesselation is valid.
 *
 * This implementation will work in two or more dimensions.  The data
 * structure returned includes not only the Dirichlet tesselation of
 * the data points, but also the Voronoi polyhedra dual to the tesselation.
 * The complete connectivity mesh of data points (called points in the
 * code), Voronoi vertices (called vertices), and their neighbor and
 * forming point relationships is maintained.
 *
 * Boyer's algorithm includes a simplex which completely encloses the
 * given data points.  The points of that simplex are the 'bogus points'
 * mentioned in the code.
 *
 * The function dch_create_dirichlet_tess creates and returns a
 * tesselation, and dch_destroy_tesselation destroys the data structure
 * returned.  
 *
 *
 * Calling Information:
 *
 * The two routines supplied by this package are:
 *
 * dch_Tess *dch_create_dirichlet_tess( float *coorddata, int npts,
 *     int dim, float *(*access_coords)(float *, int, P_Void_ptr *) )
 *
 * void dch_destroy_tesselation ___(( dch_Tess *tess_to_destroy ))

 * Tesselations are created with dch_create_dirichlet_tess, and
 * the structure returned is destroyed with dch_destroy_tesselation.
 * Calling parameters to dch_create_dirichlet_tess are used to
 * create the data points to tesselate, as follows.  dim specifies the
 * dimensionality of the space;  it must be 2 or greater.  npts specifies
 * the number of points to be tesselated.  access_coords should be a
 * pointer to a function defined as follows:
 *
 * float *access_coords( float *coorddata, int ipt, P_Void_ptr *hook )
 *
 * It is called for each ipt from 0 through npts-1, and in each case
 * should return a pointer to an array of dim floats representing 
 * coordinates.  The floats are immediately recopied, so it is not
 * necessary to allocate memory for them.  hook is a pointer to a 
 * P_Void_ptr which the access_coords routine may set;  the value
 * set is stored in the 'user_hook' field of the resulting vertex
 * and may be used to point back to the user's data.  (Note that
 * the 'user_hook' fields of the bogus points are left null).  The
 * vertices of the tesselation also have 'user_hook' fields which are
 * available to the user, though these must be set explicitly in the
 * returned tesselation.  The 'deleted' fields of the vertices of the
 * tesselation will all be 0 when the tesselation is returned, and
 * can also be used by the calling routine.
 *
 * dch_destroy_tesselation deallocates all memory associated with the
 * tesselation.  Note that it does not deallocate any memory which the
 * user may have hung on the various 'user_hook' fields.
 *
 *
 * Implementation Details:
 *
 * Boyer's algorithm includes a simplex which completely encloses the
 * given data points.  The points of that simplex are the 'bogus points'
 * mentioned in the code.
 *
 * Each point has associated with it a list of neighboring points, and
 * a list of vertices of which the point is a forming point.  Each vertex
 * has associated with it a list of forming points, and a list of 
 * neighboring vertices.  Note that the list of neighboring vertices
 * will be null unless MAINTAIN_VTX_NEIGHBOR_LISTS is defined;  this
 * is done to avoid the significant memory use associated with the
 * neighbor lists.
 *
 * The first step of Boyer's algorithm for adding a point to a valid
 * Dirichlet tesselation requires that a vertex which is closer to the
 * new points than to its forming points be found.  This is accomplished
 * by walking the lattice of already-existing points until the point
 * closest to the new point is found, and then searching the vertices of
 * that (existing) point to find a deleted vertex.  (Note that the vertex
 * closest to the new point is not necessarily deleted).  The starting
 * point for the walk is either the most recently added point or the
 * point closest to the center of the tesselation.  Which is used is
 * determined by trying both for the first few points ('few' being set
 * by the constant TRIAL_SEARCHES) and then using the method which works
 * best over that sample.
 *
 * Some of the vertices (those with n bogus points as forming points,
 * where n is the dimensionality of space) are located at infinity.
 * This is denoted by storing null pointers in their vtx->coords slots,
 * but with the vertex's 'degenerate' flag not set.  Degenerate vertices
 * may also form if the data contains cyclic collections of n+1 points
 * all lieing in a hyperplane of the n dimensional space.  Such vertices
 * also have null vtx->coords slots, but their 'degenerate' flag is set.
 * The forming points of such vertices form Delaunay simplices with
 * zero volume.
 *
 * Memory usage:
 *
 * Experience has shown that for n random points in 3D, memory usage
 * is roughly as follows:
 *
 *   points: n
 *   vertices: about 6.5 n
 *   point list cells: about 40 n
 *   vertex list cells: about 60 n if MAINTAIN_VTX_NEIGHBOR_LISTS is
 *                      defined; about 30 n otherwise.
 *   point pair cells: less than 50 (doesn't vary with n)
 *
 * The code provides for vertex pair cells, but doesn't use them.
 * 
 */

#include <math.h>
#include <stdio.h>
#include "p3dgen.h"
#include "ge_error.h"
#include "dirichlet.h"

/* Notes:
   -merging of vertices may never be necessary
*/

/* Uncomment the following to cause vertex neighbor lists to be maintained */
/*
*/
#define MAINTAIN_VTX_NEIGHBOR_LISTS

/* Some parameters */
#define SIMPLEX_MARGIN 0.1 /* buffer around points in outer bounding simplex */
#define TRIAL_SEARCHES 20 /* number of times to try both methods of finding
                               the closest point to a new point */
#define PROXIMITY_TOLERANCE 0.000001 /* factor used to determine if points
                                        or vertices are too close together */

/* Dimensionality, and some globals to help memory management */
static int current_dim= 0;
static int current_vtx_id= 0;
static int current_pt_id= 0;
static int created_vtx_count= 0;
static int created_vtx_highwater= 0;
static int created_pt_count= 0;
static int created_pt_highwater= 0;
static int created_vtx_list_count= 0;
static int created_vtx_list_highwater= 0;
static int created_pt_list_count= 0;
static int created_pt_list_highwater= 0;
static int created_vtx_pair_count= 0;
static int created_vtx_pair_highwater= 0;
static int created_pt_pair_count= 0;
static int created_pt_pair_highwater= 0;
static int dataset_size;

/* Space to perform Gaussian eliminations in */
static float *gauss_matrix= (float *)0;

/* Points from which to hang free cell lists */
static dch_Vtx_list *free_vtx_list_cells= (dch_Vtx_list *)0;
static dch_Pt_list *free_pt_list_cells= (dch_Pt_list *)0;
static dch_Vtx *free_vtx_cells= (dch_Vtx *)0; /* connected via 
						 'neighbors' slots */

/* Some forward definitions */
static void destroy_vtx_list( dch_Vtx_list * );
static void destroy_pt_list( dch_Pt_list * );
static void destroy_vtx( dch_Vtx * );
static void destroy_pt( dch_Pt * );

static void mem_init(int npts, int dim)
{
  ger_debug("mem_init: %d pts of dimension %d",npts,dim);

  dataset_size= npts;
  current_dim= dim;
  current_vtx_id= 0;
  current_pt_id= 0;

  if (gauss_matrix) free( (P_Void_ptr)gauss_matrix );
  if ( !(gauss_matrix= (float *)malloc( dim*(dim+1)*sizeof(float) )) )
    ger_fatal("mem_init: unable to allocate %d floats!",dim*(dim+1));
}

static float distance( float *coords1, float *coords2 )
/* This actually returns the squared distance, since it is only 
 * for comparison.
 */
{
  int i;
  float result= 0.0;

  if(!coords1 || !coords2) ger_fatal("distance: Null coords in distance!");
  for (i=0; i<current_dim; i++) {
    result += ( *coords1 - *coords2 ) * ( *coords1 - *coords2 );
    coords1++;
    coords2++;
  }
  return(result);
}

static float dist_to_center( float *coords, dch_Tess *tess )
/* This actually returns the squared distance to the center of the
 * tesselation's bounding box, since it is only for comparison.
 */
{
  int i;
  float result= 0.0;
  float term;

  for (i=0; i<current_dim; i++) {
    term= *coords - 0.5*( tess->bndbx->corner1[i] + tess->bndbx->corner2[i] );
    result += term*term;
    coords++;
  }
  return(result);
}

static int too_close( float *coords1, float *coords2, dch_Tess *tess )
{
  return( ( distance(coords1,coords2) <= 
           PROXIMITY_TOLERANCE*tess->characteristic_length ) );
}

static dch_Vtx *create_vtx( float *coords )
{
  dch_Vtx *vtx;
  int i;

  if (!free_vtx_cells) {
    if ( !(vtx= (dch_Vtx *)malloc(sizeof(dch_Vtx))) )
      ger_fatal("create_vtx: unable to allocate %d bytes!",sizeof(dch_Vtx));
  }
  else {
    vtx= free_vtx_cells;
    free_vtx_cells= (dch_Vtx *)vtx->neighbors;
  }

  if (coords) { /* don't do this for vertices at infinity */
    if ( !(vtx->coords= (float *)malloc(current_dim * sizeof(float))) )
      ger_fatal("create_vtx: unable to allocate %d floats!",current_dim);
    for (i=0; i<current_dim; i++) vtx->coords[i]= coords[i];
  }
  else vtx->coords= (float *)0;

  vtx->distance= 0.0;
  vtx->deleted= 0;
  vtx->degenerate= 0;
  vtx->forming_pts= (dch_Pt_list *)0;
  vtx->neighbors= (dch_Vtx_list *)0;
  vtx->id= current_vtx_id++;
  vtx->user_hook= (P_Void_ptr)0;
  created_vtx_count++;
  if (created_vtx_count > created_vtx_highwater) 
    created_vtx_highwater= created_vtx_count;

  return (vtx);
}

static void add_list_vtx( dch_Vtx_list **vlist, dch_Vtx *vtx )
/* Add the given vertex to the given list */
{
  dch_Vtx_list *temp;

  if ( !free_vtx_list_cells ) { /* need to create a cell */
    if ( !(temp= (dch_Vtx_list *)malloc(sizeof(dch_Vtx_list))) )
      ger_fatal("add_list_vtx: unable to allocate %d bytes!",
		sizeof(dch_Vtx_list));
  }
  else {
    temp= free_vtx_list_cells;
    free_vtx_list_cells= temp->next;
  }

  temp->next= *vlist;
  temp->vtx= vtx;
  *vlist= temp;
  created_vtx_list_count++;
  if (created_vtx_list_count > created_vtx_list_highwater) 
    created_vtx_list_highwater= created_vtx_list_count;

}

static void add_pair_vtxs( dch_Vtx_pair_list **vplist, 
			  dch_Vtx *vtx1, dch_Vtx *vtx2 )
{
  dch_Vtx_pair_list *temp;

  if ( !(temp= (dch_Vtx_pair_list *)malloc(sizeof(dch_Vtx_pair_list))) )
    ger_fatal("add_pair_vtxs: unable to allocate %d bytes!",
              sizeof(dch_Vtx_pair_list));

  temp->next= *vplist;
  temp->vtx1= vtx1;
  temp->vtx2= vtx2;
  *vplist= temp;
  created_vtx_pair_count++;
  if (created_vtx_pair_count > created_vtx_pair_highwater) 
    created_vtx_pair_highwater= created_vtx_pair_count;

}

static void remove_list_vtx( dch_Vtx_list **vlist, dch_Vtx *vtx )
/* This routine removes the given vertex from the list */
{
  dch_Vtx_list **set_addr;
  dch_Vtx_list *current, *holder;

  set_addr= vlist;
  current= *vlist;
  while (current) {
    if (current->vtx == vtx) {
      holder= current->next;
      current->next= free_vtx_list_cells;
      free_vtx_list_cells= current;
      *set_addr= holder;
      created_vtx_list_count--;
      return;
    }
    set_addr= &(current->next);
    current= current->next;
  }
}

static dch_Vtx *pop_vtx( dch_Vtx_list **vlist )
/* This routine pops the first vertex off the list, destroying its
 * cell but preserving the rest of the list.
 */
{
  dch_Vtx *tmp_vtx;
  dch_Vtx_list *tmp_list;

  tmp_vtx= (*vlist)->vtx;
  tmp_list= *vlist;
  *vlist= (*vlist)->next;
  tmp_list->next= free_vtx_list_cells;
  free_vtx_list_cells= tmp_list;
  created_vtx_list_count--;
  return( tmp_vtx );
}

static void destroy_vtx_list( dch_Vtx_list *vlist )
/* This destroys the vtx list, leaving the vertices intact. */
{
  dch_Vtx_list *temp,*temp2;

  ger_debug("destroy_vtx_list");

  temp= vlist;
  while (temp) {
    temp2= temp->next;
    temp->next= free_vtx_list_cells;
    free_vtx_list_cells= temp;
    temp= temp2;
    created_vtx_list_count--;
  }
}

static void destroy_vtx_pairs( dch_Vtx_pair_list *vplist )
/* This destroys the vertex pair list, leaving the vertices intact. */
{
  dch_Vtx_pair_list *temp,*temp2;

  ger_debug("destroy_vtx_pairs");

  temp= vplist;
  while (temp) {
    temp2= temp->next;
    free( (P_Void_ptr)temp );
    temp= temp2;
    created_vtx_pair_count--;
  }
}

static void destroy_vtx( dch_Vtx *vtx )
{
  if (!vtx) return;

  destroy_pt_list( vtx->forming_pts );
  destroy_vtx_list( vtx->neighbors );
    
  if (vtx->coords) free( (P_Void_ptr)vtx->coords );
  vtx->neighbors= (dch_Vtx_list *)free_vtx_cells;
  free_vtx_cells= vtx;
  created_vtx_count--;
}

static int vtx_in_list( dch_Vtx_list *vlist, dch_Vtx *vtx )
/* Returns true if the given vertex is in the sorted list, false otherwise */
{
  while (vlist) {
    if (vlist->vtx == vtx) return(1);
    vlist= vlist->next;
  }
  return(0);
}

static dch_Pt *create_pt( float *coords )
{
  dch_Pt *pt;
  int i;

  if ( !(pt= (dch_Pt *)malloc(sizeof(dch_Pt))) )
    ger_fatal("create_pt: unable to allocate %d bytes!",sizeof(dch_Vtx));
  if ( !(pt->coords= (float *)malloc(current_dim * sizeof(float))) )
    ger_fatal("create_pt: unable to allocate %d floats!",current_dim);

  for (i=0; i<current_dim; i++) pt->coords[i]= coords[i];

  pt->neighbors= (dch_Pt_list *)0;
  pt->verts= (dch_Vtx_list *)0;
  pt->user_hook= (P_Void_ptr)0;
  pt->id= current_pt_id++;
  created_pt_count++;
  if (created_pt_count > created_pt_highwater) 
    created_pt_highwater= created_pt_count;

  return (pt);
}

static void destroy_pt( dch_Pt *pt )
{
  if (!pt) ger_fatal("Tried to free a null point; algorithm error!");

  destroy_pt_list( pt->neighbors );
  destroy_vtx_list( pt->verts );

  free( (P_Void_ptr)pt );
  created_pt_count--;
}

static void add_list_pt( dch_Pt_list **plist, dch_Pt *pt )
/* Add the given point to the given list */
{
  dch_Pt_list *temp;

  if (!free_pt_list_cells) { /* need to allocate a new cell */
    if ( !(temp= (dch_Pt_list *)malloc(sizeof(dch_Pt_list))) )
      ger_fatal("add_list_pt: unable to allocate %d bytes!",
		sizeof(dch_Pt_list));
  }
  else {
    temp= free_pt_list_cells;
    free_pt_list_cells= temp->next;
  }

  temp->next= *plist;
  temp->pt= pt;
  *plist= temp;
  created_pt_list_count++;
  if (created_pt_list_count > created_pt_list_highwater) 
    created_pt_list_highwater= created_pt_list_count;

}

static void remove_list_pt( dch_Pt_list **plist, dch_Pt *pt )
/* This routine removes the given point from the list */
{
  dch_Pt_list **set_addr;
  dch_Pt_list *current, *holder;

  set_addr= plist;
  current= *plist;
  while (current) {
    if (current->pt == pt) {
      holder= current->next;
      current->next= free_pt_list_cells;
      free_pt_list_cells= current;
      *set_addr= holder;
      created_pt_list_count--;
      return;
    }
    set_addr= &(current->next);
    current= current->next;
  }
}

static void add_pair_pts( dch_Pt_pair_list **pplist, dch_Pt *pt1, dch_Pt *pt2 )
{
  dch_Pt_pair_list *temp;

  if ( !(temp= (dch_Pt_pair_list *)malloc(sizeof(dch_Pt_pair_list))) )
    ger_fatal("add_pair_pts: unable to allocate %d bytes!",
              sizeof(dch_Pt_pair_list));

  temp->next= *pplist;
  temp->pt1= pt1;
  temp->pt2= pt2;
  *pplist= temp;
  created_pt_pair_count++;
  if (created_pt_pair_count > created_pt_pair_highwater) 
    created_pt_pair_highwater= created_pt_pair_count;

}

static dch_Pt *pop_pt( dch_Pt_list **plist )
/* This routine pops the first point off the list, destroying its
 * cell but preserving the rest of the list.
 */
{
  dch_Pt *tmp_pt;
  dch_Pt_list *tmp_list;

  tmp_pt= (*plist)->pt;
  tmp_list= *plist;
  *plist= (*plist)->next;
  tmp_list->next= free_pt_list_cells;
  free_pt_list_cells= tmp_list;
  created_pt_list_count--;
  return( tmp_pt );
}

static void destroy_pt_list( dch_Pt_list *plist )
/* This destroys the point list, leaving the points intact. */
{
  dch_Pt_list *temp,*temp2;

  ger_debug("destroy_pt_list");

  temp= plist;
  while (temp) {
    temp2= temp->next;
    temp->next= free_pt_list_cells;
    free_pt_list_cells= temp;
    temp= temp2;
    created_pt_list_count--;
  }
}

static void destroy_pt_pairs( dch_Pt_pair_list *pplist )
/* This destroys the point pair list, leaving the points intact. */
{
  dch_Pt_pair_list *temp,*temp2;

  ger_debug("destroy_pt_pairs");

  temp= pplist;
  while (temp) {
    temp2= temp->next;
    free( (P_Void_ptr)temp );
    temp= temp2;
    created_pt_pair_count--;
  }
}

static int pt_in_list( dch_Pt_list *plist, dch_Pt *pt )
/* Returns true if the given point is in the sorted list, false otherwise */
{
  while (plist) {
    if (plist->pt == pt) return(1);
    plist= plist->next;
  }
  return(0);
}

static void dump_pt_list( dch_Pt_list *plist )
/* This routine is for debugging only */
{
  fprintf(stderr,"points (");
  while (plist) {
    if (plist->pt) fprintf(stderr," %d",plist->pt->id);
    else ger_fatal("Tried to dump a null point; algorithm error!");
    plist= plist->next;
  }
  fprintf(stderr," )\n");
}

static void dump_vtx_list( dch_Vtx_list *vlist )
/* This routine is for debugging only */
{
  fprintf(stderr,"vertices (");
  while (vlist) {
    fprintf(stderr," %d",vlist->vtx->id);
    vlist= vlist->next;
  }
  fprintf(stderr," )\n");
}

static void dump_vtx( dch_Vtx *vtx )
/* This routine is for debugging only */
{
  int i;
  if (!vtx) fprintf(stderr,"Attempt to dump null vertex!\n");
  else {
    fprintf(stderr,"Vertex %d:\n",vtx->id);
    if (vtx->coords) {
      fprintf(stderr,"  coords ( ");
      for (i=0; i<current_dim; i++) fprintf(stderr,"%f ",vtx->coords[i]);
      fprintf(stderr,")\n");
    }
    else fprintf(stderr,"  coords infinite\n");
    fprintf(stderr,"  distance %f; deleted= %d, degenerate= %d\n",
	    vtx->distance,vtx->deleted,vtx->degenerate);
    fprintf(stderr,"  neighbors: "); dump_vtx_list(vtx->neighbors);
    fprintf(stderr,"  forming points: "); dump_pt_list(vtx->forming_pts);
  }
}

static void dump_pt( dch_Pt *pt )
/* This routine is for debugging only */
{
  int i;
  if (!pt) fprintf(stderr,"Attempt to dump null point!\n");
  else {
    fprintf(stderr,"Point %d:\n",pt->id);
    fprintf(stderr,"  coords ( ");
    for (i=0; i<current_dim; i++) fprintf(stderr,"%f ",pt->coords[i]);
    fprintf(stderr,")\n");
    fprintf(stderr,"  neighbors: "); dump_pt_list(pt->neighbors);
    fprintf(stderr,"  vertices: "); dump_vtx_list(pt->verts);
  }
}

static dch_Bndbx *create_bndbx()
{
  dch_Bndbx *result;

  ger_debug("create_bndbx:");

  if ( !(result= (dch_Bndbx *)malloc(sizeof(dch_Bndbx))) )
    ger_fatal("create_bndbx: unable to allocate %d bytes!",sizeof(dch_Bndbx));
  if ( !(result->corner1= (float *)malloc(current_dim*sizeof(float))) )
    ger_fatal("create_bndbx: unable to allocate %d floats!",current_dim);
  if ( !(result->corner2= (float *)malloc(current_dim*sizeof(float))) )
    ger_fatal("create_bndbx: unable to allocate %d floats!",current_dim);

  result->empty_flag= 1;

  return(result);
}

static void destroy_bndbx( dch_Bndbx *bndbx )
{
  ger_debug("destroy_bndbx:");

  free( (P_Void_ptr)bndbx->corner1 );
  free( (P_Void_ptr)bndbx->corner2 );
  free( (P_Void_ptr)bndbx );
}

static void dump_bndbx( dch_Bndbx *bndbx )
{
  int i;

  ger_debug("dump_bndbx:");

  if (bndbx->empty_flag) fprintf(stderr,"Bounding box: empty\n");
  else {
    fprintf(stderr,"Bounding box: corner1 ( ");
    for (i=0; i<current_dim; i++) fprintf(stderr,"%f ",bndbx->corner1[i]);
    fprintf(stderr,")\n");
    fprintf(stderr,"              corner2 ( ");
    for (i=0; i<current_dim; i++) fprintf(stderr,"%f ",bndbx->corner2[i]);
    fprintf(stderr,")\n");
  }
}

static void add_pt_to_bndbx( dch_Bndbx *bndbx, dch_Pt *pt )
{
  int i;

  if (bndbx->empty_flag) {
    for (i=0; i<current_dim; i++) 
      bndbx->corner1[i]= bndbx->corner2[i]= pt->coords[i];
    bndbx->empty_flag= 0;
  }
  else for (i=0; i<current_dim; i++) {
    if (bndbx->corner1[i] > pt->coords[i]) bndbx->corner1[i]= pt->coords[i];
    if (bndbx->corner2[i] < pt->coords[i]) bndbx->corner2[i]= pt->coords[i];
  }
}

static int pt_in_bndbx( dch_Bndbx *bndbx, dch_Pt *pt )
{
  int i;
  for (i=0; i<current_dim; i++) {
    if (bndbx->corner1[i] > pt->coords[i]) return(0);
    if (bndbx->corner2[i] < pt->coords[i]) return(0);
  }
  return(1);
}

static void dump_tesselation( dch_Tess *tess )
/* This routine is for debugging only */
{
  dch_Vtx_list *vlist;
  dch_Pt_list *plist;

  fprintf(stderr,"***Dumping tesselation***\n");
  fprintf(stderr,"  dimensionality %d, characteristic length %f\n",
	  tess->dimensionality,tess->characteristic_length);
  dump_bndbx( tess->bndbx );
  fprintf(stderr,"Point list contents follow:\n");
  plist= tess->pt_list;
  while (plist) {
    dump_pt( plist->pt );
    plist= plist->next;
  }
  fprintf(stderr,"Vertex list contents follow:\n");
  vlist= tess->vtx_list;
  while (vlist) {
    dump_vtx( vlist->vtx );
    vlist= vlist->next;
  }
  fprintf(stderr,"Bogus point list: ");
  dump_pt_list(tess->bogus_pts);
  fprintf(stderr,"Current center point is %d\n",tess->center_pt->id);
  fprintf(stderr,"Current most recent point is %d\n", 
	  tess->most_recent_pt->id);
  fprintf(stderr,"***end of dump***\n");
}

void dch_dump_tesselation(dch_Tess*tess)
{
  dump_tesselation(tess);
}

static dch_Pt_list *calc_bogus_pts( dch_Bndbx *bndbx )
/* This routine calculates the corners of a simplex completely enclosing
 * the given bounding box.
 */
{
  float edge_length= 0.0;
  float *coords;
  int i,j;
  dch_Pt_list *result= (dch_Pt_list *)0;

  ger_debug("calc_bogus_pts:");

  /* Malloc some space for needed coordinates */
  if ( !(coords= (float *)malloc(current_dim*sizeof(float))) )
    ger_fatal("calc_bogus_pts: unable to allocate %d floats!",current_dim);

  /* Find the maximum edge length; assume an n-cube of that edge length. */
  for (i=0; i<current_dim; i++)
    if ( (bndbx->corner2[i] - bndbx->corner1[i]) > edge_length )
      edge_length= bndbx->corner2[i] - bndbx->corner1[i];

  /* First corner is the minimum bounding box corner minus SIMPLEX_MARGIN
   * times the edge length.
   */
  for (i=0; i<current_dim; i++) 
    coords[i]= bndbx->corner1[i] - SIMPLEX_MARGIN*edge_length;
  add_list_pt( &result, create_pt(coords) );

  /* Other corners are (1 + SIMPLEX_MARGIN)*n*edge length up the coordinate
   * axes, where n is the dimensionality of the space.
   */
  for (i=0; i<current_dim; i++) {
    for (j=0; j<current_dim; j++) 
      coords[j]= bndbx->corner1[j] - SIMPLEX_MARGIN*edge_length;
    coords[i]= 
      bndbx->corner1[i] + (1.0+SIMPLEX_MARGIN)*current_dim*edge_length;
    add_list_pt( &result, create_pt(coords) );
  }

  /* Clean up */
  free( (P_Void_ptr)coords );

  return( result );
}

static void connect_neighboring_pts( dch_Pt_list *pts, dch_Vtx *vtx )
/* This routine takes a set of forming points and their vertex, makes
 * all the points neighbors of eachother, and adds the vertex to each
 * point's vertex list.
 */
{
  dch_Pt_list *thispt,*thatpt;

  thispt= pts;
  while (thispt) {
    add_list_vtx( &(thispt->pt->verts), vtx );
    thatpt= thispt->next;
    while (thatpt) {
      if ( (thispt->pt != thatpt->pt) &&
          !pt_in_list(thispt->pt->neighbors, thatpt->pt) ) {
        add_list_pt( &(thispt->pt->neighbors),thatpt->pt );
        add_list_pt( &(thatpt->pt->neighbors),thispt->pt );
      }
      thatpt= thatpt->next;
    }
    thispt= thispt->next;
  }
}

static float calc_product_term( dch_Pt_list *thispt )
/* This routine is used in find_central_vtx */
{
  float result= 0.0;
  int i;

  for (i=0; i<current_dim; i++) 
    result += (thispt->pt->coords[i]*thispt->pt->coords[i])
             -(thispt->next->pt->coords[i]*thispt->next->pt->coords[i]);
  return( result/2.0 );
}

static dch_Vtx *find_central_vtx( dch_Pt_list *pts )
/* This routine returns a vertex equidistant from all the given points. */
{
  int i,j;
  dch_Pt_list *thispt;
  float *coords;
  dch_Vtx *result;

  ger_debug("find_central_vtx:");

  /* We want a point equidistant from the given points.  This can
   * be found by solving the matrix equation (e.g. in 2D, with
   * forming points A, B, and C and solution P):
   *
   *  / Ax-Bx    Ay-By \  /Px}\     / (A**2 - B**2)/2 \
   * |                  ||     | = |                   |
   *  \ Bx-Cx    By-Cy /  \Py}/     \ (B**2 - C**2)/2 /
   *
   * where A**2 and B**2 are the squared norms of A and B.  This
   * is easy to prove by simply writing the requirement that A and
   * B be equidistant from P, expanding the squares, and simplifying.
   * The algorithm relies on the fact that there will be current_dim+1
   * points in the input point list.
   */

  thispt= pts;
  for (i=0; i<current_dim; i++) {
    for (j=0; j<current_dim; j++)
      *(gauss_matrix + i*(current_dim+1) + j)=
        thispt->pt->coords[j] - thispt->next->pt->coords[j];
    *(gauss_matrix + i*(current_dim+1) + current_dim)=
      calc_product_term(thispt);
    thispt= thispt->next;
  }

  /* Find the vertex coordinates by Gaussian elimination.  If no
   * solution can be found, we create a vertex at infinty and mark
   * it for deletion.
   */
  coords= dch_gauss_elim( gauss_matrix, current_dim );
  result= create_vtx( coords );

  result->forming_pts= pts;
  if (coords) result->distance= distance( coords, pts->pt->coords );
  else result->degenerate= 1;

  /* The forming points of the new vertex are now all neighbors, and
   * need to know about their new vertex */
  connect_neighboring_pts( pts, result );

  free( (P_Void_ptr)coords );

  return(result);
}

static void add_first_vtx_neighbors( dch_Vtx *vtx, dch_Tess *tess )
/* This routine creates neighbor vertices at infinity, and appropriately
 * hooks them up the the newly formed first vertex.
 */
{
  dch_Vtx *inf_vtx;
  dch_Pt **point_array;
  dch_Pt_list *plist;
  int i,j;

  /* Generate an array to hold the current_dim+1 bogus points, to
   * simplify constructing the point lists for forming points.
   */
  if ( !(point_array= (dch_Pt **)malloc( (current_dim+1)*sizeof(dch_Pt *) )) )
    ger_fatal("add_first_vtx_neighbors: unable to allocate %d bytes!",
              (current_dim+1)*sizeof(dch_Pt *));
  
  plist= tess->bogus_pts;
  for (i=0; i<(current_dim+1); i++ ) {
    point_array[i]= plist->pt;
    plist= plist->next;
  }

  /* We create current_dim+1 vertices at infinity, each with the
   * known vertex as a neighbor and current_dim of the bogus
   * points as forming points.  This leaves each forming point
   * list one vertex short.
   */
  for (i=0; i<(current_dim+1); i++) {
    inf_vtx= create_vtx( (float *)0 ); /* create vertex at infinity */
    add_list_vtx( &(tess->infinite_vtxs), inf_vtx );
    for (j=0; j<(current_dim+1); j++) 
      if (i != j) {
	add_list_pt( &(inf_vtx->forming_pts), point_array[j] );
	add_list_vtx( &(point_array[j]->verts), inf_vtx );
      }
#ifdef MAINTAIN_VTX_NEIGHBOR_LISTS
    add_list_vtx( &(inf_vtx->neighbors), vtx );
    add_list_vtx( &(vtx->neighbors), inf_vtx );
#endif
  }

  /* Clean up */
  free( (P_Void_ptr)point_array );
}

static dch_Tess *create_initial_tesselation( dch_Bndbx *bndbx )
{
  dch_Tess *tess;
  dch_Pt_list *plist, *forming_pt_list= (dch_Pt_list *)0;
  dch_Vtx *new_vtx;

  ger_debug("create_initial_tesselation:");

  if ( !(tess= (dch_Tess *)malloc(sizeof(dch_Tess))) )
    ger_fatal("create_initial_tesselation: unable to allocate %d bytes!",
              sizeof(dch_Tess));

  tess->dimensionality= current_dim;
  tess->bndbx= bndbx;
  tess->vtx_list= (dch_Vtx_list *)0;
  tess->infinite_vtxs= (dch_Vtx_list *)0;
  tess->pt_list= (dch_Pt_list *)0;

  /* Calculate the bogus points, which define the outer bound of
   * a simplex guaranteed to hold all the data points.  Copy the
   * resulting point list as the initial pt_list.
   */
  plist= tess->bogus_pts= calc_bogus_pts( tess->bndbx );
  tess->pt_list= (dch_Pt_list *)0;
  while (plist) {
    add_list_pt( &(tess->pt_list), plist->pt );
    plist= plist->next;
  }

  /* Calculate a first vertex based on the bogus points.
   * We need another copy of the list, as it will become the forming
   * point list of the new vertex and we have to keep separate lists
   * for memory management purposes. 
   */
  plist= tess->bogus_pts;
  while (plist) {
    add_list_pt( &(forming_pt_list), plist->pt );
    plist= plist->next;
  }
  new_vtx= find_central_vtx( forming_pt_list );

  /* The new vertex needs current_dim + 1 neighbors at infinity */
  add_first_vtx_neighbors( new_vtx, tess );

  /* The center point we set is not necessarily really closest
   * to the center, but the next point insertion will replace
   * it with an appropriate value.
   */
  tess->most_recent_pt= tess->center_pt= tess->pt_list->pt;

  /* The following data is used to pick an appropriate search method
   * for finding points closest to newly added points.
   */
  tess->best_search_method= UNKNOWN;
  tess->searches_done= 0;
  tess->center_steps= 0;
  tess->most_recent_steps= 0;

  /* The following data is used with PROXIMITY_TOLERANCE to determine 
   * when to drop points or vertices because they are too close together.
   */
  tess->characteristic_length= new_vtx->distance;

  return(tess);
}

void dch_destroy_tesselation( dch_Tess *tess )
{
  dch_Pt_list *plist;
  dch_Vtx_list *vlist;

  ger_debug("destroy_tesselation:");

  destroy_bndbx( tess->bndbx );

  plist= tess->pt_list;
  while (plist) {
    destroy_pt( plist->pt );
    plist= plist->next;
  }

  vlist= tess->vtx_list;
  while (vlist) {
    destroy_vtx( vlist->vtx );
    vlist= vlist->next;
  }

  destroy_pt_list( tess->bogus_pts );
  destroy_pt_list( tess->pt_list );
  destroy_vtx_list( tess->vtx_list );
  destroy_vtx_list( tess->infinite_vtxs );

  free( (P_Void_ptr)tess );
#ifdef never
  fprintf(stderr,"Final counts: %d %d %d %d %d %d %d %d %d %d %d %d\n",
	  created_vtx_count,
	  created_vtx_highwater,
	  created_pt_count,
	  created_pt_highwater,
	  created_vtx_list_count,
	  created_vtx_list_highwater,
	  created_pt_list_count,
	  created_pt_list_highwater,
	  created_vtx_pair_count,
	  created_vtx_pair_highwater,
	  created_pt_pair_count,
	  created_pt_pair_highwater );
#endif

}

static dch_Pt *find_closest_pt( dch_Pt *start, dch_Pt *target, int *steps )
/* This routine walks the vertex list from the given starting point
 * in search of a vertex to be deleted.  steps is a counter, to be
 * used in picking a method for selecting the best starting point.
 */
{
  dch_Pt_list *neighbors;
  dch_Pt *best_pt, *new_pt;
  float best_dist, new_dist;
  int getting_better= 1;

  ger_debug("find_closest_pt: starting pt %d",start->id);

  best_pt= start;
  best_dist= distance( best_pt->coords, target->coords );
  while (getting_better) {
    getting_better= 0;
    neighbors= best_pt->neighbors;
    while (neighbors) {
      new_pt= neighbors->pt;
      if ( (new_dist= distance( new_pt->coords, target->coords ))
          < best_dist ) {
        getting_better= 1;
        best_pt= new_pt;
        best_dist= new_dist;
      }
      neighbors= neighbors->next;
    }    
  *steps += 1;
  }

  return( best_pt );
}

static dch_Vtx *find_vtx_to_delete( dch_Tess *tess, dch_Pt *pt)
{
  dch_Pt *closest_pt;
  dch_Vtx_list *vlist;
  float tol;

  ger_debug("find_vtx_to_delete:");

  /* Find the closest point.  If we don't know what starting point is
   * best, try both and accumulate statistics.  Otherwise, use the 
   * one that's been best so far.
   */
  switch ((int)(tess->best_search_method)) {
  case (int)UNKNOWN: 
    closest_pt= find_closest_pt( tess->center_pt, pt, &(tess->center_steps));
    (void) find_closest_pt( tess->most_recent_pt, pt, 
                           &(tess->most_recent_steps)); /* throw out result */
    if (tess->searches_done >= TRIAL_SEARCHES) {
      if (tess->center_steps < tess->most_recent_steps) 
        tess->best_search_method= CENTER_PT;
      else tess->best_search_method= MOST_RECENT_PT;
    }
    break;
  case (int)MOST_RECENT_PT: 
    closest_pt= find_closest_pt( tess->most_recent_pt, pt, 
                                &(tess->most_recent_steps));
    break;
  case (int)CENTER_PT: 
    closest_pt= find_closest_pt( tess->center_pt, pt, &(tess->center_steps));
    break;
  default: ger_fatal("find_vtx_to_delete: unknown method %d!\n");
  }
  tess->searches_done++;

  /* If the new point is too close to the existing point, return null
   * to signal that this point should be dropped from the tesselation.
   */
  if ( too_close( closest_pt->coords, pt->coords, tess ) ) 
    return( (dch_Vtx *)0 );

  /* At least one of the closest point's vertices will need to be deleted. 
   * We have to watch out for vertices at infinity (which can come from
   * degeneracies);  they will eventually be cleaned up anyway so just
   * skip over them.  We also have to watch out for the possibility of
   * numerical noise making a vertex look like it shouldn't be deleted
   * when in fact it should.
   */
  vlist= closest_pt->verts;
  tol= PROXIMITY_TOLERANCE*tess->characteristic_length;
  while (vlist) {
    if ( vlist->vtx->coords
	&& (distance( vlist->vtx->coords, pt->coords ) 
	    < vlist->vtx->distance + tol) ) {
      vlist->vtx->deleted= 1;
      return( vlist->vtx );
    }
    vlist= vlist->next;
  }

  /* If we get to here, we failed to find a deleted vertex. */
  ger_fatal("find_vtx_to_delete: no deleted vertex found!");
  return( (dch_Vtx *)0 ); /* to satisfy lint */
}

static void find_all_deleted_vtxs( dch_Vtx_list **deleted_vtx_list, 
				  dch_Pt *pt, dch_Tess *tess )
/* Check all neighbors of the single deleted vertex in the input
 * list, and add them to the list.  We must handle the cases
 * where some of the neighbors are at infinity or are degenerate,
 * and the case where numerical noise might make a vertex look like
 * it shouldn't be deleted.
 */
{
  dch_Vtx_list *stack= (dch_Vtx_list *)0, *neighbors;
  dch_Vtx *vtx;
  dch_Pt_list *forming_pts;
  float tol;

  ger_debug("find_all_deleted_vtxs:");

  tol= PROXIMITY_TOLERANCE*tess->characteristic_length;
  add_list_vtx( &stack, (*deleted_vtx_list)->vtx );
  while (stack) {
    vtx= pop_vtx( &stack );
#ifndef MAINTAIN_VTX_NEIGHBOR_LISTS
    forming_pts= vtx->forming_pts;
    while (forming_pts) {
      neighbors= forming_pts->pt->verts;
      while (neighbors) {
	if ( !(neighbors->vtx->deleted)
	    && ( neighbors->vtx->degenerate
		|| (neighbors->vtx->coords
		    && (distance( neighbors->vtx->coords, pt->coords ) <
			neighbors->vtx->distance + tol) ) ) ) {
	  neighbors->vtx->deleted= 1;
	  add_list_vtx( deleted_vtx_list, neighbors->vtx );
	  add_list_vtx( &stack, neighbors->vtx );
	}
	neighbors= neighbors->next;
      }
      forming_pts= forming_pts->next;
    }
#else
    neighbors= vtx->neighbors;
    while (neighbors) {
      if ( !(neighbors->vtx->deleted)
          && ( neighbors->vtx->degenerate
	      || (neighbors->vtx->coords
		  && (distance( neighbors->vtx->coords, pt->coords ) <
		      neighbors->vtx->distance + tol) ) ) ) {
        neighbors->vtx->deleted= 1;
        add_list_vtx( deleted_vtx_list, neighbors->vtx );
        add_list_vtx( &stack, neighbors->vtx );
      }
      neighbors= neighbors->next;
    }
#endif
  }
}

static dch_Pt_list *accumulate_forming_pts( dch_Vtx_list *deleted_vtx_list )
/* Search the deleted vertex list, collecting (uniquely) all the forming
 * points of the deleted vertices.
 */
{
  dch_Pt_list *loose_pts= (dch_Pt_list *)0;
  dch_Vtx *vtx;
  dch_Pt_list *forming_pts;

  ger_debug("accumulate_forming_pts:");

  while (deleted_vtx_list) {
    vtx= deleted_vtx_list->vtx;
    forming_pts= vtx->forming_pts;
    while (forming_pts) {
      if ( !pt_in_list(loose_pts,forming_pts->pt) )
        add_list_pt( &loose_pts, forming_pts->pt );
      forming_pts= forming_pts->next;
    }
    deleted_vtx_list= deleted_vtx_list->next;
  }

  return(loose_pts);
}

static int shared_live_vtx( dch_Pt *pt1, dch_Pt *pt2 )
/* This routine returns true if the two points share an undeleted
 * vertex, and false otherwise.
 */
{
  dch_Vtx_list *verts;

  ger_debug("shared_live_vtx: checking points %d and %d",pt1->id,pt2->id);

  verts= pt1->verts;
  while (verts) {
    if ( !(verts->vtx->deleted) && vtx_in_list( pt2->verts, verts->vtx ) ) 
      return(1);
    verts= verts->next;
  }

  return(0);
}

static void remove_common_contiguities( dch_Pt_list *loose_pt_list )
/* This routine searches the neighbor lists of the given loose point
 * list for neighbor points that are also on the list.  If such a
 * pair of points is found, it is deleted if the two points share
 * no undeleted vertex.
 */
{
  dch_Pt *pt, *neighbor_pt;
  dch_Pt_list *neighbors;
  dch_Pt_pair_list *to_be_deleted= (dch_Pt_pair_list *)0, *pplist;

  ger_debug("remove_common_contiguities:");

  /* First we find all the pairs to be deleted. */
  while (loose_pt_list) {
    pt= loose_pt_list->pt;
    neighbors= pt->neighbors;
    while (neighbors) {
      neighbor_pt= neighbors->pt;
      if ( pt_in_list(loose_pt_list,neighbor_pt)
          && !shared_live_vtx(pt,neighbor_pt) ) {
        add_pair_pts( &to_be_deleted, pt, neighbor_pt );
      }
      neighbors= neighbors->next;
    }
    loose_pt_list= loose_pt_list->next;
  }

  /* Now we delete the contiguities found */
  pplist= to_be_deleted;
  while (pplist) {
    remove_list_pt( &(pplist->pt1->neighbors), pplist->pt2 );
    remove_list_pt( &(pplist->pt2->neighbors), pplist->pt1 );
    pplist= pplist->next;
  }

  /* Clean up */
  destroy_pt_pairs( to_be_deleted );
}

static dch_Pt_list *find_common_forming_pts( dch_Vtx *vtx1, dch_Vtx *vtx2 )
/* This routine searches the forming points lists of the two vertices
 * for forming points which they have in common, returning a list
 * of the points found.
 */
{
  dch_Pt_list *plist= (dch_Pt_list *)0;
  dch_Pt_list *v1_fpts, *v2_fpts;

  ger_debug("find_common_forming_pts: vertices %d and %d",vtx1->id, vtx2->id);

  v1_fpts= vtx1->forming_pts;
  while (v1_fpts) {
    v2_fpts= vtx2->forming_pts;
    while (v2_fpts) {
      if ( v1_fpts->pt == v2_fpts->pt ) add_list_pt( &plist, v1_fpts->pt );
      v2_fpts= v2_fpts->next;
    }
    v1_fpts= v1_fpts->next;
  }

  return(plist);
}

static void remove_all_incoming_connections( dch_Vtx *vtx )
/* Remove the deleted vertex from the vertex lists of all its 
 * forming points, and from the neighbor lists of all its
 * undeleted neighbors.  The vertex will be destroyed shortly,
 * so we don't have to worry about its own forming point or
 * neighbor lists.
 */
{
  dch_Pt_list *forming_pts;
  dch_Vtx_list *neighbors;

  ger_debug("remove_all_incoming_connections: vertex %d",vtx->id);

  forming_pts= vtx->forming_pts;
  while (forming_pts) {
    remove_list_vtx( &(forming_pts->pt->verts), vtx );
    forming_pts= forming_pts->next;
  }

#ifdef MAINTAIN_VTX_NEIGHBOR_LISTS
  neighbors= vtx->neighbors;
  while (neighbors) {
    if ( !(neighbors->vtx->deleted) )
      remove_list_vtx( &(neighbors->vtx->neighbors), vtx );
    neighbors= neighbors->next;
  }
#endif
}

static int count_different_forming_pts( dch_Vtx *vtx1, dch_Vtx *vtx2 )
/* This routine runs through the forming point lists of the two given
 * vertices, looking for points in one which are not in the other.
 */
{
  dch_Pt_list *fp1, *fp2;
  int diffs_found= 0;

  ger_debug("count_different_forming_pts: checking vertices %d and %d",
	    vtx1->id, vtx2->id);

  fp1= vtx1->forming_pts;
  fp2= vtx2->forming_pts;
  while (fp1) {
    if ( !pt_in_list( fp2, fp1->pt ) ) diffs_found++;
    fp1= fp1->next;
  }
  return( diffs_found );
}

static dch_Vtx_list *find_vtx_neighbors( dch_Vtx *vtx )
/* This routine is used to find a list of vertex neighbors.  It is needed
 * only if vertex neighbor lists are not being dynamically maintained.
 */
{
  dch_Vtx_list *neighbors= (dch_Vtx_list *)0, *fp_verts;
  dch_Pt_list *forming_pts;

  ger_debug("find_vtx_neighbors: vertex %d",vtx->id);

  forming_pts= vtx->forming_pts;
  while (forming_pts) {
    fp_verts= forming_pts->pt->verts;
    while (fp_verts) {
      if ( !(vtx_in_list( neighbors, fp_verts->vtx ))
	  && count_different_forming_pts( vtx, fp_verts->vtx )==1 )
	add_list_vtx( &neighbors, fp_verts->vtx );
      fp_verts= fp_verts->next;
    }
    forming_pts= forming_pts->next;
  }

  return( neighbors );
}

static dch_Vtx_list *test_and_generate_new_vtxs( dch_Vtx_list *deleted_vtxs, 
						dch_Pt *pt )
/* Every deleted vertex is now searched for neighbor vertices which are
 * not deleted.  Each such pair of a deleted and an undeleted vertex
 * defines a line, along which lies a new vertex.  That vertex is found
 * by finding the common forming points of the deleted and undeleted
 * vertices, and using them and the new point as the forming points
 * of a new vertex.  The new and undeleted vertices become neighbors.
 * Sometimes this process creates degenerate vertices;  we keep track
 * of them and return them.
 */
{
  dch_Vtx_list *neighbors;
  dch_Vtx_list *created_vtxs= (dch_Vtx_list *)0;
  dch_Pt_list *forming_pts;
  dch_Vtx *vtx;
  dch_Vtx *nbr_vtx;
  dch_Vtx *new_vtx;

  ger_debug("test_and_generate_new_vtx: point %d",pt->id);

  /* Walk the list, doing the given checks to all vertices */
  while (deleted_vtxs) {
    vtx= deleted_vtxs->vtx;

    /* Get the neighbor list of the given vertex, and check it for
     * undeleted vertices.  
     */
#ifdef MAINTAIN_VTX_NEIGHBOR_LISTS
    neighbors= vtx->neighbors;
#else
    neighbors= vtx->neighbors= find_vtx_neighbors( vtx );
#endif
    while (neighbors) {
      nbr_vtx= neighbors->vtx;
      if ( !(nbr_vtx->deleted) ) {
	
	/* Collect the common forming points of the two vertices */
	forming_pts= find_common_forming_pts( vtx, nbr_vtx );
	add_list_pt( &forming_pts, pt );
	
	/* Create a vertex with the given forming points.  The forming
	 * points become the forming points list of the new vertex,
	 * and the new vertex becomes a vertex of all the forming points,
	 * which also become neighbors.
	 */
	new_vtx= find_central_vtx( forming_pts );

	/* Keep a list of the new vertices */
	add_list_vtx( &created_vtxs, new_vtx );

	/* Make the new vertex and the neighbor neighbors */
#ifdef MAINTAIN_VTX_NEIGHBOR_LISTS
	add_list_vtx( &(new_vtx->neighbors), nbr_vtx );
	add_list_vtx( &(nbr_vtx->neighbors), new_vtx );
#endif
      }
      
      /* Continue with the surrounding loop */
      neighbors= neighbors->next;
    }

    /* Disconnect the deleted vertex from its forming points and neighbors */
    remove_all_incoming_connections( vtx );

    deleted_vtxs= deleted_vtxs->next;
  }

  return ( created_vtxs );
}

static void merge_vtxs( dch_Vtx *vtx1, dch_Vtx *vtx2 )
/* This routine merges the two given vertices by making sure that all
 * the neighbors of the second vertex are neighbors of the first and
 * then destroying the second vertex.
 */
{
  dch_Vtx_list *neighbors;

  ger_debug("merge_vtxs: merging vertices %d and %d",vtx1->id,vtx2->id);

#ifdef MAINTAIN_VTX_NEIGHBOR_LISTS
  neighbors= vtx2->neighbors;
  while (neighbors) {
    if ( !vtx_in_list( vtx1->neighbors, neighbors->vtx ) ) {
      add_list_vtx( &(vtx1->neighbors), neighbors->vtx );
      add_list_vtx( &(neighbors->vtx->neighbors), vtx1 );
    }
    neighbors= neighbors->next;
  }
#endif

  remove_all_incoming_connections( vtx2 );
  destroy_vtx( vtx2 );
}

static void connect_neighboring_vtxs( dch_Vtx_list *vlist )
/* This routine takes a vertex list and searches it for pairs of
 * vertices which share all but one common forming point.  Such
 * pairs must lie on the same line, so they are joined as neighbors.
 */
{
  dch_Vtx *vtx1;
  dch_Vtx_list *vlist_rest, *holder;
  int different_pt_count;

  ger_debug("connect_neighboring_vtxs:");

  while (vlist && vlist->next) {
    vtx1= vlist->vtx;
    vlist_rest= vlist->next;
    while (vlist_rest) {
      holder= vlist_rest->next;
      different_pt_count= count_different_forming_pts(vtx1, vlist_rest->vtx);

      if (different_pt_count==0) {
	/* A degeneracy has produced two vertices with the same forming
	 * points; merge them.
	 */
	remove_list_vtx( &vlist, vlist_rest->vtx );
	merge_vtxs( vtx1, vlist_rest->vtx );
      }

      else if (different_pt_count==1) {
	/* These two vertices become neighbors */
#ifdef MAINTAIN_VTX_NEIGHBOR_LISTS
	add_list_vtx( &(vtx1->neighbors), vlist_rest->vtx );
	add_list_vtx( &(vlist_rest->vtx->neighbors), vtx1 );
#endif
      }

      vlist_rest= holder;
    }
    vlist= vlist->next;
  }
}

static dch_Vtx_list *accumulate_vtxs( dch_Pt_list *plist )
/* This routine collects all the vertices of the points.  The vertices'
 * 'deleted' flag is used to avoid double counting. 
 */
{
  dch_Vtx_list *vlist= (dch_Vtx_list *)0, *vl;

  ger_debug("accumulate_vtxs:");

  /* Accumulate vertices by searching the vertices associated with
   * each point, marking each deleted as it is found.  Only unmarked
   * vertices get accumulated, so each only gets put on the list
   * once.
   */
  while (plist) {
    vl= plist->pt->verts;
    while (vl) {
      if ( !(vl->vtx->deleted) ) {
	add_list_vtx( &(vlist), vl->vtx );
	vl->vtx->deleted= 1;
	if ( vl->vtx->degenerate )
	  ger_debug("accumulate_vtxs: found degenerate vertex %d",
		    vl->vtx->id);
	else if ( !(vl->vtx->coords) )
	  ger_debug("accumulate_vtxs: found infinite vertex %d", vl->vtx->id);
      }
      vl= vl->next;
    }
    plist= plist->next;
  }

  /* Run down the new vertex list, turning off all the 'deleted' flags. */
  vl= vlist;
  while (vl) {
    vl->vtx->deleted= 0;
    vl= vl->next;
  }

  return( vlist );
}

static void add_pt_to_tesselation( dch_Tess *tess, dch_Pt *pt )
/* This function adds a point to the tesselation, returning the current
 * list of degenerate vertices.
 */
{
  dch_Pt_list *loose_pt_list= (dch_Pt_list *)0;
  dch_Vtx_list *thisvtx, *deleted_vtx_list= (dch_Vtx_list *)0;
  dch_Vtx_list *created_vtx_list= (dch_Vtx_list *)0;
  dch_Vtx *deleted_vtx;

  ger_debug("add_pt_to_tesselation: adding point %d",pt->id);

  if ( !pt_in_bndbx( tess->bndbx, pt ) ) {
    ger_fatal("add_pt_to_tesselation: point outside bounding box; ignored.\n");
    return;
  }

  /* Find a vertex in the current tesselation that must be deleted. */
  deleted_vtx= find_vtx_to_delete( tess, pt );

  /* A null list means that find_vtx_to_delete discovered a point which
   * was too close to the new point (possibly causing algorithm errors);
   * the fix is to drop the new point.  The solution is to drop the
   * new point.  This constitutes a memory leak, because the new
   * point will never be garbage collected, but that should not be
   * a problem since bad points are (presumably) rare.
   */
  if (!deleted_vtx) {
    ger_error(
      "add_pt_to_tesselation: Point %d dropped; too close to existing point.",
      pt->id);
    return;
  }
  
  /* Add the new point to the tesselation's point list */
  add_list_pt( &(tess->pt_list), pt );

  /* Generate a complete list of deleted vertices by tree search, starting
   * with the known deleted vertex. 
   */
  add_list_vtx( &deleted_vtx_list, deleted_vtx );
  find_all_deleted_vtxs( &deleted_vtx_list, pt, tess );
  
  /* Accumulate a list of the forming points of deleted vertices */
  loose_pt_list= accumulate_forming_pts( deleted_vtx_list );

  /* Remove contiguities shared by points, both of which are in the
   * given point list.
   */
  remove_common_contiguities( loose_pt_list );

  /* We walk the deleted vertex list, looking for links to undeleted
   * vertices.  Each such link defines a new vertex to be created.  The
   * old vertices are disconnected from the mesh of points and vertices.
   */
  created_vtx_list= test_and_generate_new_vtxs( deleted_vtx_list, pt );
  
  /* The vertices of the new point now need to be connected to eachother
   * in appropriate neighbor relationships.
   */
  connect_neighboring_vtxs( created_vtx_list );

  /* Make the new point the tesselation's 'most recent point',
   * and make it the center point if it is closer to the
   * center than the current center point.
   */
  tess->most_recent_pt= pt;
  if ( dist_to_center(pt->coords,tess) 
      < dist_to_center(tess->center_pt->coords,tess) )
    tess->center_pt= pt;

  /* Clean up.  First walk the deleted vertex list, freeing all the
   * vertices there.  Then free the deleted vertex list, created 
   * vertex list and loose point list.
   */
  thisvtx= deleted_vtx_list;
  while( thisvtx ) {
    destroy_vtx( thisvtx->vtx );
    thisvtx= thisvtx->next;
  }
  destroy_pt_list( loose_pt_list );
  destroy_vtx_list( created_vtx_list );
  destroy_vtx_list( deleted_vtx_list );

}

static dch_Pt *create_user_pt( float *coorddata, int i,
  float *(*coord_access_fun)( float *, int, P_Void_ptr *) )
/* This function just invokes the user's coordinate access function */
{
  P_Void_ptr uhook= (P_Void_ptr)0;
  dch_Pt *result;

  result= create_pt( (*coord_access_fun)( coorddata, i, &uhook ) );
  result->user_hook= uhook;
  return( result );
}

dch_Tess *dch_create_dirichlet_tess
  ( float *coorddata, int npts, int dimensionality,
   float *(*coord_access_fun)(float *, int, P_Void_ptr *) )
{
  int i;
  dch_Pt_list *thispt, *input_pt_list= (dch_Pt_list *)0;
  dch_Bndbx *bndbx;
  dch_Tess *tess;

  ger_debug("create_dirichlet_tesselation:");

  if (npts<1) {
    ger_error("create_dirichlet_tesselation: %d is not enough input points!",
	      npts);
    return( (dch_Tess *)0 );
  }

  mem_init( npts, dimensionality );

  /* Convert all the coordinate data to a point list. */
  for (i=0; i<npts; i++)
    add_list_pt( &input_pt_list, 
		create_user_pt( coorddata, i, coord_access_fun ) );

  /* Calculate the bounding box */
  bndbx= create_bndbx();
  thispt= input_pt_list;
  while (thispt) {
    add_pt_to_bndbx( bndbx, thispt->pt );
    thispt= thispt->next;
  }

  /* Create the initial tesselation */
  tess= create_initial_tesselation(bndbx);

  /* Save the access function, for completeness' sake */
  tess->coord_access_fun= coord_access_fun;

  /* Add points, one at a time */
  thispt= input_pt_list;
  while (thispt) {
    add_pt_to_tesselation( tess, thispt->pt );
    thispt= thispt->next;
  }

  /* Accumulate the (non-infinite) vertices */
  tess->vtx_list= accumulate_vtxs( tess->pt_list );

  /* Clean up */
  destroy_pt_list( input_pt_list );

  return(tess);
}

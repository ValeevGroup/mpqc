/****************************************************************************
 * rand_isosurf.c
 * Author Joel Welling
 * Copyright 1990, Pittsburgh Supercomputing Center, Carnegie Mellon University
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
This module provides the ability to generate isosurfaces from randomly
distributed data.
*/
#include <stdio.h>
#include "p3dgen.h"
#include "pgen_objects.h"
#include "ge_error.h"
#include "dirichlet.h"

/* Notes-
   -does ival want to be a double, to please ornery C compilers?
   -I don't think we need to distinguish between cut_links and cuts,
    though it may be a good idea for future flexibility.
*/

#define CVTX_SZ    3
#define CCVTX_SZ   7
#define CNVTX_SZ   6
#define CCNVTX_SZ  10
#define CVVTX_SZ   4
#define CVNVTX_SZ  7

/* A cell to hold information about a triangle, and a counter for them */
typedef struct triangle_struct {
  int i1, i2, i3;
  dch_Pt *outside_pt;
  dch_Pt *inside_pt;
  struct triangle_struct *next;
} Triangle;
static int live_triangles= 0;

/* A cell to hold information about a cut, and a counter for them */
typedef struct cut_struct {
  dch_Pt *p1;
  dch_Pt *p2;
  int id;
  struct cut_struct *next;
} Cut;
static int live_cuts= 0;

/* A cell to hang on a dch_Pt, holding cut information, and a hook to
 * connect all of them to.
 */
typedef struct cut_link_struct {
  Cut *cut;
  struct cut_link_struct *next;
  struct cut_link_struct *all_links_next;
} Cut_link;
static Cut_link *cut_link_list= (Cut_link *)0;

/* A structure to connect to the user_hook field of the dch_Pt's in
 * the Dirichlet tesselation, and a hook for the global array of them.
 */
typedef struct hook_data_struct {
  int index; 
  float value;
  Cut_link *cut_links;
} Hook_data;
static Hook_data *hook_array;

/* Macros to aid point data access */
#define PT_VALUE( pt ) (((Hook_data *)((pt)->user_hook))->value)
#define PT_INDEX( pt ) (((Hook_data *)((pt)->user_hook))->index)
#define PT_CUTS( pt ) (((Hook_data *)((pt)->user_hook))->cut_links)

/* Hook on which to hang the vertex list being tesselated */
static P_Vlist *current_vlist= (P_Vlist *)0;

static void add_triangle( Triangle **list, dch_Pt *in_pt, dch_Pt *out_pt,
			 int i1, int i2, int i3 )
/* This routine adds a triangle to the given list */
{
  Triangle *newcell;

  if ( !(newcell= (Triangle *)malloc( sizeof(Triangle) )) )
    ger_fatal("rand_isosurf: add_triangle: unable to allocate %d bytes!",
	      sizeof(Triangle));

  newcell->outside_pt= out_pt;
  newcell->inside_pt= in_pt;
  newcell->i1= i1;
  newcell->i2= i2;
  newcell->i3= i3;
  newcell->next= *list;
  *list= newcell;
  live_triangles++;
}

static void destroy_triangle_list( Triangle *list )
/* This routine destroys the given triangle list */
{
  Triangle *victim, *next_victim;

  victim= list;
  while (victim) {
    next_victim= victim->next;
    free( (P_Void_ptr)victim );
    live_triangles--;
    victim= next_victim;
  }
}

static Cut *create_and_list_cut( Cut **list, 
				dch_Pt *p1, dch_Pt *p2 )
/* This routine creates a cut and adds it to the given list */
{
  Cut *newcell;

  if ( !(newcell= (Cut *)malloc( sizeof(Cut) )) )
    ger_fatal("rand_isosurf: create_and_list_cut: can't allocate %d bytes!",
	      sizeof(Cut));
  live_cuts++;

  newcell->p1= p1;
  newcell->p2= p2;
  newcell->id= 0;
  newcell->next= *list;
  *list= newcell;
  return( newcell );
}

static void destroy_cut_list( Cut *list )
/* This routine destroys all the cuts in the list */
{
  Cut *next;

  while (list) {
    next= list->next;
    free( (P_Void_ptr)list );
    live_cuts--;
    list= next;
  }
}

static void destroy_cut_link_list( Cut_link *list )
/* This routine destroys all the cut links in the list (traversing
 * the all_links_next connections). 
 */
{
  Cut_link *next;
  
  while (list) {
    next= list->all_links_next;
    free( (P_Void_ptr)list );
    list= next;
  }
}

static Cut *find_cut_edges( dch_Tess *tess, float ival )
/* This routine walks the tesselation, looking for edges which are cut by
 * the new isosurface.
 */
{
  Cut *cut_list= (Cut *)0;
  dch_Pt_list *plist, *neighbors;
  float val1, val2;

  ger_debug("rand_isosurf: find_cut_edges: finding cuts at %f", ival);

  /* Walk the point list, looking for cuts */
  plist= tess->pt_list;
  while (plist) {
    if (plist->pt->user_hook) { /* not a 'bogus point' */
      val1= PT_VALUE( plist->pt );
      neighbors= plist->pt->neighbors;
      while (neighbors) {
	if ( neighbors->pt->user_hook   /* not a 'bogus point' */
	    && (plist->pt->id < neighbors->pt->id) /* avoid double-counting */
	    ) { 
	  val2= PT_VALUE( neighbors->pt );
	  if ( ((val1 >= ival) && (val2 < ival)) 
	      || ((val2 >= ival) && (val1 < ival)) ) /* cut this edge */
	    (void)create_and_list_cut( &cut_list, plist->pt, neighbors->pt );
	}
	neighbors= neighbors->next;
      }
    }
    plist= plist->next;
  }

  return( cut_list );
}

static void interpolate_data( dch_Pt *p1, dch_Pt *p2, float ival, 
			     float *data, int vtxtype )
/* This routine fills the given data slot with appropriate data based
 * an interpolation of vertex data of the given type.
 */
{
  float v1, v2, fraction;
  int v1_index, v2_index;

  v1= PT_VALUE( p1 );
  v2= PT_VALUE( p2 );
  fraction= (ival-v1) / (v2-v1);
  if ( fraction<0.0 || fraction>1.0 )
    ger_fatal("rand_isosurf: interpolate_data: algorithm error; fraction= %f!",
	      fraction);

  /* Get the interpolated point coordinates */
  *data++= (1.0-fraction) * p1->coords[0] + fraction * p2->coords[0];
  *data++= (1.0-fraction) * p1->coords[1] + fraction * p2->coords[1];
  *data++= (1.0-fraction) * p1->coords[2] + fraction * p2->coords[2];

  if (vtxtype==P3D_CVVTX) return; /* no need to do any more */

  /* Recover reference data left in hooks by coord_access, and ready the
   * vertex list object
   */
  v1_index= PT_INDEX(p1);
  v2_index= PT_INDEX(p2); 
  METHOD_RDY( current_vlist );

  /* Get the interpolated point values if necessary */
  if (vtxtype==P3D_CVVVTX) {
    *data++= (1.0-fraction) * (*(current_vlist->v2))(v1_index)
      + fraction * (*(current_vlist->v2))(v2_index);
  }

  /* Get the interpolated point normals if necessary */
  if (vtxtype==P3D_CVNVTX) {
    *data++= (1.0-fraction) * (*(current_vlist->nx))(v1_index)
      + fraction * (*(current_vlist->nx))(v2_index);
    *data++= (1.0-fraction) * (*(current_vlist->ny))(v1_index)
      + fraction * (*(current_vlist->ny))(v2_index);
    *data++= (1.0-fraction) * (*(current_vlist->nz))(v1_index)
      + fraction * (*(current_vlist->nz))(v2_index);
  }
}

static void attach_cut_point( dch_Pt *p1, dch_Pt *p2, Cut *cut )
/* This routine attaches the given cut to the given point pair within
 * the tesselation data structure.  We attach the cut pointer to the
 * Hook_data field of the point with the higher id, where it joins
 * a potentially growing list.
 */
{
  Cut_link *link;
  dch_Pt *target;

  /* Allocate the link */
  if ( !(link= (Cut_link *)malloc(sizeof(Cut_link))) )
    ger_fatal("rand_isosurf: attach_cut_point: unable to allocate %d bytes!",
	      sizeof(Cut_link));
  link->all_links_next= cut_link_list;
  cut_link_list= link;

  /* Establish the connection to the cut */
  link->cut= cut;

  /* Connect it to the proper point */
  if (p1->id == p2->id)
    ger_fatal("rand_isosurf: attach_cut_point: algorithm error!");
  if (p1->id > p2->id) target= p1;
  else target= p2;
  link->next= PT_CUTS(target);
  PT_CUTS(target)= link;
}

static int get_cut_id( dch_Pt *p1, dch_Pt *p2 )
/* This routine recovers the id of the cut associated with the given point
 * pair within the tesselation data structure.  The information is stored
 * in a list of cut links attached to the point with the highest id.
 */
{
  Cut_link *link;
  dch_Pt *target;

  /* Find the point to search, and its list of cuts */
  if (p1->id == p2->id)
    ger_fatal("rand_isosurf: get_cut_id: algorithm error; points match!");
  if (p1->id > p2->id) target= p1;
  else target= p2;
  link= PT_CUTS(target);

  /* Run down the list until we find the right match */
  while (link) {
    if ( ((link->cut->p1 == p1) && (link->cut->p2 == p2))
	|| ((link->cut->p1 == p2) && (link->cut->p2 == p1)) )
      return( link->cut->id );
    link= link->next;
  }

  /* Getting to here means we failed to find the link */
  ger_fatal("rand_isosurf: get_cut_id: algorithm error; \
link between points %d and %d not found!",
	    p1->id, p2->id);

  return(0); /* to satisfy lint */
}

static int generate_cut_edges( Cut *cut_list, float *cuts, int table_step,
			      float ival, int vtxtype )
/* This routine walks the cut list, generating cut point data needed
 * for the isosurface.
 */
{
  Cut *runner;
  int id= 0;

  ger_debug("rand_isosurf: generate_cut_edges: generating cuts at %f", ival);

  /* Walk the cut list, generating vertex information of the correct type */
  runner= cut_list;
  while (runner) {
    interpolate_data( runner->p1, runner->p2, ival, cuts, vtxtype );
    runner->id= id++;
    cuts += table_step;
    runner= runner->next;
  }

  /* Now patch the cut information into the tesselation structure.  We
   * wait until the end to do this, as the storage method makes old
   * point data inaccessible.
   */
  runner= cut_list;
  while (runner) {
    attach_cut_point( runner->p1, runner->p2, runner );
    runner= runner->next;
  }
}

static void handle_one_tet( dch_Vtx *vtx, Triangle **list, float ival )
{
  int vtx_case=0, pow2 = 1;
  dch_Pt_list *fplist;
  dch_Pt *pts[4];
  int i;

  /* Drop uninteresting (degenerate, deleted, or infinite) tets immediately */
  if ( vtx->degenerate || vtx->deleted || !vtx->coords ) return;

  fplist= vtx->forming_pts;
  for (i=0; i<4; i++) { /* guaranteed four forming points in 3D */
    if ( PT_VALUE(fplist->pt) >= ival ) vtx_case += pow2;
    pow2 *= 2;
    pts[i]= fplist->pt;
    fplist= fplist->next;
  }

  /* Always build triangles with a vertex which is greater than the
   * isosurface value first, so check_polarity can know the proper
   * orientation.
   */
  switch (vtx_case) {
  case 0: break;
  case 1: 
    add_triangle( list, pts[0], pts[1],
		 get_cut_id( pts[0], pts[1] ),
		 get_cut_id( pts[0], pts[2] ),
		 get_cut_id( pts[0], pts[3] ) );
    break;
  case 2: 
    add_triangle( list, pts[1], pts[0],
		 get_cut_id( pts[1], pts[0] ),
		 get_cut_id( pts[1], pts[3] ),
		 get_cut_id( pts[1], pts[2] ) );
    break;
  case 3: 
    add_triangle( list, pts[0], pts[3],
		 get_cut_id( pts[0], pts[3] ),
		 get_cut_id( pts[1], pts[3] ),
		 get_cut_id( pts[0], pts[2] ) );
    add_triangle( list, pts[1], pts[2],
		 get_cut_id( pts[1], pts[2] ),
		 get_cut_id( pts[0], pts[2] ),
		 get_cut_id( pts[1], pts[3] ) );
    break;
  case 4: 
    add_triangle( list, pts[2], pts[0],
		 get_cut_id( pts[2], pts[0] ),
		 get_cut_id( pts[2], pts[1] ),
		 get_cut_id( pts[2], pts[3] ) );
    break;
  case 5: 
    add_triangle( list, pts[0], pts[1],
		 get_cut_id( pts[0], pts[1] ),
		 get_cut_id( pts[2], pts[1] ),
		 get_cut_id( pts[0], pts[3] ) );
    add_triangle( list, pts[2], pts[3],
		 get_cut_id( pts[2], pts[3] ),
		 get_cut_id( pts[0], pts[3] ),
		 get_cut_id( pts[2], pts[1] ) );
    break;
  case 6: 
    add_triangle( list, pts[1], pts[3],
		 get_cut_id( pts[1], pts[3] ),
		 get_cut_id( pts[2], pts[3] ),
		 get_cut_id( pts[1], pts[0] ) );
    add_triangle( list, pts[2], pts[0],
		 get_cut_id( pts[2], pts[0] ),
		 get_cut_id( pts[1], pts[0] ),
		 get_cut_id( pts[2], pts[3] ) );
    break;
  case 7: 
    add_triangle( list, pts[1], pts[3],
		 get_cut_id( pts[1], pts[3] ),
		 get_cut_id( pts[2], pts[3] ),
		 get_cut_id( pts[0], pts[3] ) );
    break;
  case 8: 
    add_triangle( list, pts[3], pts[0],
		 get_cut_id( pts[3], pts[0] ),
		 get_cut_id( pts[3], pts[2] ),
		 get_cut_id( pts[3], pts[1] ) );
    break;
  case 9: 
    add_triangle( list, pts[0], pts[2],
		 get_cut_id( pts[0], pts[2] ),
		 get_cut_id( pts[3], pts[2] ),
		 get_cut_id( pts[0], pts[1] ) );
    add_triangle( list, pts[3], pts[1],
		 get_cut_id( pts[3], pts[1] ),
		 get_cut_id( pts[0], pts[1] ),
		 get_cut_id( pts[3], pts[2] ) );
    break;
  case 10: 
    add_triangle( list, pts[1], pts[0],
		 get_cut_id( pts[1], pts[0] ),
		 get_cut_id( pts[3], pts[0] ),
		 get_cut_id( pts[1], pts[2] ) );
    add_triangle( list, pts[3], pts[2],
		 get_cut_id( pts[3], pts[2] ),
		 get_cut_id( pts[1], pts[2] ),
		 get_cut_id( pts[3], pts[0] ) );
    break;
  case 11: 
    add_triangle( list, pts[1], pts[2],
		 get_cut_id( pts[0], pts[2] ),
		 get_cut_id( pts[3], pts[2] ),
		 get_cut_id( pts[1], pts[2] ) );
    break;
  case 12: 
    add_triangle( list, pts[2], pts[1],
		 get_cut_id( pts[2], pts[1] ),
		 get_cut_id( pts[3], pts[1] ),
		 get_cut_id( pts[2], pts[0] ) );
    add_triangle( list, pts[3], pts[0],
		 get_cut_id( pts[3], pts[0] ),
		 get_cut_id( pts[2], pts[0] ),
		 get_cut_id( pts[3], pts[1] ) );
    break;
  case 13: 
    add_triangle( list, pts[0], pts[1],
		 get_cut_id( pts[0], pts[1] ),
		 get_cut_id( pts[3], pts[1] ),
		 get_cut_id( pts[2], pts[1] ) );
    break;
  case 14: 
    add_triangle( list, pts[1], pts[0],
		 get_cut_id( pts[2], pts[0] ),
		 get_cut_id( pts[1], pts[0] ),
		 get_cut_id( pts[3], pts[0] ) );
    break;
  case 15: break;
  }

}

static void check_orientation( Triangle *tri, float *cuts, int step, 
			      int show_inside )
/* This routine flips the orientation of the triangle if its back
 * face is pointing the wrong way.
 */
{
  float gx, gy, gz, e1x, e1y, e1z, e2x, e2y, e2z, cx, cy, cz, dot;
  int temp;

  /* From the inside and outside points of the triangle, generate
   * a vector pointing up the gradient (i.e. toward the inside of
   * the isosurface).
   */
  gx= tri->inside_pt->coords[0] - tri->outside_pt->coords[0];
  gy= tri->inside_pt->coords[1] - tri->outside_pt->coords[1];
  gz= tri->inside_pt->coords[2] - tri->outside_pt->coords[2];

  /* Generate the right-hand-rule normal of the triangle by taking
   * the cross product of the edges.
   */
  e1x= *( cuts + step * tri->i2 ) - *( cuts + step * tri->i1 );
  e1y= *( cuts + step * tri->i2 + 1 ) - *( cuts + step * tri->i1 + 1 );
  e1z= *( cuts + step * tri->i2 + 2 ) - *( cuts + step * tri->i1 + 2 );
  e2x= *( cuts + step * tri->i3 ) - *( cuts + step * tri->i2 );
  e2y= *( cuts + step * tri->i3 + 1 ) - *( cuts + step * tri->i2 + 1 );
  e2z= *( cuts + step * tri->i3 + 2 ) - *( cuts + step * tri->i2 + 2 );
  cx= e1y*e2z - e1z*e2y;
  cy= e1z*e2x - e1x*e2z;
  cz= e1x*e2y - e1y*e2x;

  /* The dot product will be positive if the normal points toward
   * the inside.  Then we flip the indices to invert the normal.
   */
  dot= gx*cx + gy*cy + gz*cz;
  if ( (show_inside && dot<0.0) || (!show_inside && dot>0.0)) {
    temp= tri->i2;
    tri->i2= tri->i3;
    tri->i3= temp;
  }
}

static Triangle *march_tets( dch_Tess *tess, float *cuts, int step,
			    float ival, int show_inside )
/* This routine actually implements the marching tets algorithm */
{
  Triangle *triangle_list= (Triangle *)0, *runner;
  dch_Vtx_list *vlist;

  ger_debug("rand_isosurf: march_tets: generating isosurface at value %f",
	    ival);

  vlist= tess->vtx_list;
  while (vlist) {
    handle_one_tet( vlist->vtx, &triangle_list, ival );
    vlist= vlist->next;
  }

  runner= triangle_list;
  while (runner) {
    check_orientation( runner, cuts, step, show_inside );
    runner= runner->next;
  }
  return( triangle_list );
}

static float get_value( P_Vlist *vlist, int i )
/* This function returns the value associated with the given point */
{
  METHOD_RDY( vlist );
  switch( current_vlist->type ) {
  case P3D_CVVTX:
  case P3D_CVNVTX:
  case P3D_CVVVTX:
    return( (*(vlist->v))( i ) );
    break;
  default:
    ger_fatal(
       "rand_isosurf: get_value: vertex type %d has no value to contour!",
	      vlist->type);
    break;
  }
  return(0.0); /* to please lint */
}

static float *coord_access( float *dummy, int i, P_Void_ptr *user_hook )
/* This routine is called by the Dirichlet tesselation package to get
 * coordinate data.  Note that the coordinate data is recopied, so we
 * can use static storage for it.
 */
{
  static float coords[3];

  METHOD_RDY( current_vlist );
  coords[0]= (*(current_vlist->x))(i);
  coords[1]= (*(current_vlist->y))(i);
  coords[2]= (*(current_vlist->z))(i);

  /* Connect and fill out the hook data cell */
  *user_hook= (P_Void_ptr)&hook_array[i];
  hook_array[i].index= i;
  hook_array[i].value= get_value( current_vlist, i );
  hook_array[i].cut_links= (Cut_link *)0;

  return(coords);
}

static int do_mesh( int vtxtype, float *cuts, Triangle *triangle_list )
/* This routine actually emits the isosurface mesh */
{
  int retval;
  int *vertices, *facet_lengths, *vptr, *flptr;

  ger_debug("do_mesh:");

  /* Allocate space for the vertex and facet length tables */
  if ( !(facet_lengths= (int *)malloc( live_triangles*sizeof(int) )) )
    ger_fatal("rand_isosurf: do_mesh: unable to allocate %d ints for lengths!",
	      live_triangles);
  if ( !(vertices= (int *)malloc( 3*live_triangles*sizeof(int) )) )
    ger_fatal("rand_isosurf: do_mesh: unable to allocate %d ints for indices!",
	      3*live_triangles);

  /* Transcribe the facet information */
  vptr= vertices;
  flptr= facet_lengths;
  while (triangle_list) {
    *flptr++= 3;
    *vptr++= triangle_list->i1;
    *vptr++= triangle_list->i2;
    *vptr++= triangle_list->i3;
    triangle_list= triangle_list->next;
  }

  /* Generate the mesh */
  retval= dp_mesh( vtxtype, P3D_RGB, cuts, live_cuts,
		  vertices, facet_lengths, live_triangles );

  /* Clean up */
  free( (P_Void_ptr)facet_lengths );
  free( (P_Void_ptr)vertices );

  return( retval );
}

static void mark_outside_vtxs( dch_Tess *tess )
/* This routine marks all tesselation vertices with at least one 'bogus'
 * point as deleted.  This amounts to marking all the Delaunay simplices
 * outside the convex hull, so they can be skipped as we accumulate
 * triangles.
 */
{
  dch_Pt_list *plist;
  dch_Vtx_list *vlist;

  ger_debug("rand_isosurf: mark_outside_vtxs:");

  plist= tess->bogus_pts; /* the list of outside forming points */
  while (plist) {
    vlist= plist->pt->verts; /* the vertices this point helps form */
    while (vlist) {
      vlist->vtx->deleted= 1;
      vlist= vlist->next;
    }
    plist= plist->next;
  }
}

int pg_rand_isosurf( P_Vlist *vlist, double ival, int show_inside )
/* Creates a random isosurface from the points in data. */
{
  dch_Tess *tess;
  Cut *cut_list;
  float *cuts;
  Triangle *triangle_list;
  int in_vtxtype, out_vtxtype;
  int r_isosurf_flag;
  int float_per_vtx;

  ger_debug(
       "p3dgen: pg_rand_isosurf: \
adding random isosurface at %f from %d vertices",
            ival, vlist->length);

  /* Quit if no gob is open */
  if (!pg_gob_open()) {
    ger_error("pg_rand_isosurf: No gob is currently open; call ignored.");
    return(P3D_FAILURE);
  }

  in_vtxtype= vlist->type;
  /* calculates the necessary size for each vertex type */
  switch( in_vtxtype ) {
  case P3D_CVVTX:
    out_vtxtype= P3D_CVTX;
    float_per_vtx = CVTX_SZ;
    break;
  case P3D_CVNVTX:
    out_vtxtype= P3D_CNVTX;
    float_per_vtx = CNVTX_SZ;
    break;
  case P3D_CVVVTX:
    out_vtxtype= P3D_CVVTX;
    float_per_vtx= CVVTX_SZ;
    break;
  default:
    ger_error(
       "rand_isosurf: vertex type %d has no value to contour; call ignored.",
	      (int)in_vtxtype);
    return(P3D_FAILURE);
    break;
  }

  /* Allocate cells to be hung from the user_hook fields of the
   * Dirichlet tesselation, one per data point.
   */
  if ( !(hook_array= (Hook_data *)malloc( vlist->length*sizeof(Hook_data) )) )
    ger_fatal("p3dgen:pg_rand_isosurf: unable to allocate %d hook data cells!",
	      vlist->length);

  /* Calculate the Dirichlet tesselation, and mark all vertices outside
   * the convex hull as deleted.
   */
  current_vlist= vlist;
  tess= dch_create_dirichlet_tess( (float *)0, vlist->length, 3, 
					 coord_access );
  mark_outside_vtxs( tess );

  /* Generate the cut list */
  cut_list= find_cut_edges( tess, ival );

  /* Allocate space for the table of cut vertices */
  if ( !(cuts = 
	 ( float * ) malloc( live_cuts*float_per_vtx*sizeof( float ) ) ) )
    ger_fatal("p3dgen: rand_zsurf: unable to allocate %d floats!",
	      float_per_vtx-1);

  /* Fill out the cut table, and in the process edit cut information
   * into the tesselation data structure.
   */
  generate_cut_edges( cut_list, cuts, float_per_vtx, ival, in_vtxtype );

  /* Uncomment this to put a marker at every interpolated vertex */
  /*
    dp_polymarker( out_vtxtype, P3D_RGB, cuts, live_cuts );
  */

  /* Run the marching tets algorithm, using the known cuts (pointers
   * to which have been plugged into tess)
   */
  triangle_list= march_tets( tess, cuts, float_per_vtx, ival, show_inside );

  /* Generate the isosurface mesh */
  r_isosurf_flag= do_mesh( out_vtxtype, cuts, triangle_list );

  /* Clean up */
  dch_destroy_tesselation(tess);
  free( (P_Void_ptr)cuts );
  free( (P_Void_ptr)hook_array );
  destroy_cut_list( cut_list );
  destroy_triangle_list( triangle_list );
  destroy_cut_link_list( cut_link_list );
  cut_link_list= (Cut_link *)0;

  return( r_isosurf_flag );
}

/****************************************************************************
 * dirichlet.h
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
 * This module provides an interface to the Dirichlet tesselation package.
 */

/* Include prototypes unless requested not to.  NO_SUB_PROTO causes
 * all function prototypes not explicitly associated with an 'extern'ed
 * function to be omitted.  NO_PROTO causes all prototypes to be omitted.
 */
#ifdef NO_PROTO
#define NO_SUB_PROTO
#endif

#ifdef NO_PROTO
#define ___(prototype) ()
#else
#define ___(prototype) prototype
#endif
#ifdef NO_SUB_PROTO
#define __(prototype) ()
#else
#define __(prototype) prototype
#endif

enum dch_search_method { MOST_RECENT_PT, CENTER_PT, UNKNOWN };

struct dch_vtx_struct;
struct dch_pt_struct;
struct dch_vtx_list_struct;
struct dch_pt_list_struct;

typedef struct dch_pt_struct {
  struct dch_pt_list_struct *neighbors;
  struct dch_vtx_list_struct *verts;
  int id;
  float *coords;
  P_Void_ptr user_hook;
} dch_Pt;

typedef struct dch_pt_list_struct {
  dch_Pt *pt;
  struct dch_pt_list_struct *next;
} dch_Pt_list;

typedef struct dch_pt_pair_struct {
  dch_Pt *pt1, *pt2;
  struct dch_pt_pair_struct *next;
} dch_Pt_pair_list;

typedef struct dch_vtx_struct {
  struct dch_pt_list_struct *forming_pts;
  struct dch_vtx_list_struct *neighbors;
  int id;
  float distance; /* square of distance, actually */
  int deleted;
  int degenerate;
  float *coords;
  P_Void_ptr user_hook;
} dch_Vtx;

typedef struct dch_vtx_list_struct {
  dch_Vtx *vtx;
  struct dch_vtx_list_struct *next;
} dch_Vtx_list;

typedef struct dch_vtx_pair_struct {
  dch_Vtx *vtx1, *vtx2;
  struct dch_vtx_pair_struct *next;
} dch_Vtx_pair_list;

typedef struct dch_bndbx_struct {
  float *corner1;
  float *corner2;
  int empty_flag;
} dch_Bndbx;

typedef struct dch_tess_struct {
  int dimensionality;
  dch_Vtx_list *vtx_list;
  dch_Pt_list *pt_list;
  dch_Pt_list *bogus_pts;
  dch_Vtx_list *infinite_vtxs;
  dch_Pt *center_pt;
  dch_Pt *most_recent_pt;
  dch_Bndbx *bndbx;
  enum dch_search_method best_search_method;
  int searches_done;
  int center_steps;
  int most_recent_steps;
  float characteristic_length;
  float *(*coord_access_fun) __((float *, int, P_Void_ptr *));
} dch_Tess;

float *dch_gauss_elim ___(( float *, int ));
dch_Tess *dch_create_dirichlet_tess
  ___(( float *, int npts, int dimensionality,
       float *(*access_coords) __((float *, int, P_Void_ptr *)) ));
void dch_destroy_tesselation ___(( dch_Tess * ));
void dch_dump_tesselation ___((dch_Tess * ));

/* clean up the prototyping macros */
#undef __
#undef ___

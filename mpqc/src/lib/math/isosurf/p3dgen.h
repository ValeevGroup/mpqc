/****************************************************************************
 * p3dgen.h
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
This include file defines the entry points provided by p3dgen.
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

/* Provide a generic pointer even if this is not standard C.  HP Unix
 * defines __STDC__ but doesn't seem to handle the void pointer correctly.
 */
#if (__STDC__ && ! __hpux)
typedef void *P_Void_ptr;
#else
typedef char *P_Void_ptr;
#endif

/* The Stardent C compiler doesn't properly handle (void) parameter lists. */
#ifdef stardent
#define VOIDLIST /* nothing */
#else
#define VOIDLIST void
#endif

/* P3dgen version */
#define P3D_P3DGEN_VERSION 2.2

/* Error codes */
#define P3D_SUCCESS 1
#define P3D_FAILURE 0

/* Booleans */
#define P3D_TRUE 1
#define P3D_FALSE 0

/* Vertex list types */
#define P3D_CVTX 0
#define P3D_CCVTX 1
#define P3D_CCNVTX 2
#define P3D_CNVTX 3
#define P3D_CVVTX 4
#define P3D_CVNVTX 5
#define P3D_CVVVTX 6

/* Color specification types */
#define P3D_RGB 0

/* Material specification values */
#define P3D_DEFAULT_MATERIAL 0
#define P3D_DULL_MATERIAL 1
#define P3D_SHINY_MATERIAL 2
#define P3D_METALLIC_MATERIAL 3
#define P3D_MATTE_MATERIAL 4
#define P3D_ALUMINUM_MATERIAL 5

/* Maximum object name length */
#define P3D_NAMELENGTH 64

/* types */
typedef struct P_Color_struct { int ctype; float r, g, b, a; } P_Color;
typedef struct P_Point_struct { float x, y, z; } P_Point;
typedef struct P_Vector_struct { float x, y, z; } P_Vector;
typedef struct P_Transform_struct { float d[16]; } P_Transform;
typedef struct P_Material_struct { int type; } P_Material;
typedef P_Void_ptr P_Symbol;

/* predefined materials */
extern P_Material *p3d_default_material;
extern P_Material *p3d_dull_material;
extern P_Material *p3d_shiny_material;
extern P_Material *p3d_metallic_material;
extern P_Material *p3d_matte_material;
extern P_Material *p3d_aluminum_material;

/* vlist objects */
typedef struct P_Vlist_struct {
  int type;                   /* one of the vertex list types */
  int length;                 /* number of vertices */
  int retained;               /* true if non-volatile copy of data exists */
  int data_valid;             /* true if there is no risk anyone has freed
				 the vertex data */
  double (*x) __((int));      /* returns x[index] */
  double (*y) __((int));      /* returns y[index] */
  double (*z) __((int));      /* returns z[index] */
  double (*nx) __((int));     /* returns normal_x[index] if present */
  double (*ny) __((int));     /* returns normal_y[index] if present */
  double (*nz) __((int));     /* returns normal_z[index] if present */
  double (*r) __((int));      /* returns r[index] if present */
  double (*g) __((int));      /* returns g[index] if present */
  double (*b) __((int));      /* returns b[index] if present */
  double (*a) __((int));      /* returns a[index] if present */
  double (*v) __((int));      /* returns v[index] if present */
  double (*v2) __((int));     /* returns v2[index] if present */
  void (*print) __(( void )); /* print method */
  void (*destroy_self) __((void)); /* destroy method */
  P_Void_ptr object_data;     /* object data */
} P_Vlist;

#ifdef __cplusplus

extern "C" P_Vlist *po_create_cvlist( int, int, const float * );

extern "C" P_Vlist *po_create_fvlist( int, int, 
                                 float *, float *, float *, 
                                 float *, float *, float *,
                                 float *, float *, float *, float * );

extern "C" P_Vlist *po_create_mvlist( int , int, const float *, 
				     const float *, const float * );

/* Overall control routines */
extern "C" int pg_initialize( void );
extern "C" int pg_shutdown( void );

/* Renderer control routines */
extern "C" int pg_init_ren( char *, char *, char *, char *);
extern "C" int pg_open_ren( char * );
extern "C" int pg_close_ren( char * );
extern "C" int pg_shutdown_ren( char * );
extern "C" int pg_print_ren( char * );

/* Gob routines */
extern "C" int pg_open( char * );
extern "C" int pg_close( void );
extern "C" int pg_free( char * );
extern "C" int pg_gob_open( void );
extern "C" int pg_int_attr( char *, int );
extern "C" int pg_bool_attr( char *, int );
extern "C" int pg_float_attr( char *, double );
extern "C" int pg_string_attr( char *, char * );
extern "C" int pg_color_attr( char *, P_Color * );
extern "C" int pg_point_attr( char *, P_Point * );
extern "C" int pg_vector_attr( char *, P_Vector * );
extern "C" int pg_trans_attr( char *, P_Transform * );
extern "C" int pg_material_attr( char *, P_Material * );
extern "C" int pg_gobcolor( P_Color * );
extern "C" int pg_textheight( double );
extern "C" int pg_backcull( int );
extern "C" int pg_gobmaterial( P_Material * );
extern "C" int pg_transform( P_Transform * );
extern "C" int pg_translate( double, double, double );
extern "C" int pg_rotate( P_Vector *, double );
extern "C" int pg_scale( double );
extern "C" int pg_ascale( double, double, double );
extern "C" int pg_child( char * );
extern "C" int pg_print_gob( char * );

/* Primitive routines */
extern "C" int pg_cylinder( void );
extern "C" int pg_sphere( void );
extern "C" int pg_torus( double, double );
extern "C" int pg_polymarker( P_Vlist * );
extern "C" int pg_polyline( P_Vlist * );
extern "C" int pg_polygon( P_Vlist * );
extern "C" int pg_tristrip( P_Vlist * );
extern "C" int pg_mesh( P_Vlist *, int *, int *, int );
extern "C" int pg_bezier( P_Vlist * ); /* always 16 */
extern "C" int pg_text( char *, P_Point *, P_Vector *, P_Vector * );
extern "C" int pg_light( P_Point *, P_Color * );
extern "C" int pg_ambient( P_Color * );

/* Composite routines */
extern "C"  int pg_axis( P_Point *, P_Point *, P_Vector *, double, 
                        double, int, char *, double, int );
extern "C" int pg_boundbox( P_Point *, P_Point * );
extern "C" int pg_isosurface( int type, float *data, float *valdata,
		  int nx, int ny, int nz, double value,
		  P_Point *corner1, P_Point *corner2,
		  int show_inside, int ftn_order );
extern "C" int pg_zsurface( int, float *, float *, 
                  int, int, P_Point *, P_Point *, 
                  void (*)(int *, float *, int *, int *), int );
extern "C" int pg_rand_zsurf(P_Vlist *vlist, 
			   void (*)(int *, float *, float *, float *, int *));
extern "C" int pg_rand_isosurf(P_Vlist *vlist, double value, int show_inside);

/* Camera */
extern "C" int pg_camera( char *, P_Point *, P_Point *,
                     P_Vector *, double, double, double );
extern "C" int pg_print_camera( char * );
extern "C" int pg_camera_background( char *, P_Color * );

/* Snap */
extern "C" int pg_snap( char *, char *, char * );

/* Color map */
extern "C" int pg_set_cmap(double, double, void (*)( float *,
						    float *, float *,
						    float *, float * ) );

extern "C" int pg_std_cmap(double, double, int);

extern "C" int pg_cmap_color(const float *, float *, float *, float *, 
				  float *);

#else /* __cplusplus not defined */

extern P_Vlist *po_create_cvlist ___(( int, int, float * ));

extern P_Vlist *po_create_fvlist ___(( int, int, 
                                 float *, float *, float *, 
                                 float *, float *, float *,
                                 float *, float *, float *, float * ));

extern P_Vlist *po_create_mvlist ___(( int , int, float *, float *, float * ));

/* Overall control routines */
extern int pg_initialize ___(( void ));
extern int pg_shutdown ___(( void ));

/* Renderer control routines */
extern int pg_init_ren ___(( char *, char *, char *, char *));
extern int pg_open_ren ___(( char * ));
extern int pg_close_ren ___(( char * ));
extern int pg_shutdown_ren ___(( char * ));
extern int pg_print_ren ___(( char * ));

/* Gob routines */
extern int pg_open ___(( char * ));
extern int pg_close ___(( void ));
extern int pg_free ___(( char * ));
extern int pg_gob_open ___(( void ));
extern int pg_int_attr ___(( char *, int ));
extern int pg_bool_attr ___(( char *, int ));
extern int pg_float_attr ___(( char *, double ));
extern int pg_string_attr ___(( char *, char * ));
extern int pg_color_attr ___(( char *, P_Color * ));
extern int pg_point_attr ___(( char *, P_Point * ));
extern int pg_vector_attr ___(( char *, P_Vector * ));
extern int pg_trans_attr ___(( char *, P_Transform * ));
extern int pg_material_attr ___(( char *, P_Material * ));
extern int pg_gobcolor ___(( P_Color * ));
extern int pg_textheight ___(( double ));
extern int pg_backcull ___(( int ));
extern int pg_gobmaterial ___(( P_Material * ));
extern int pg_transform ___(( P_Transform * ));
extern int pg_translate ___(( double, double, double ));
extern int pg_rotate ___(( P_Vector *, double ));
extern int pg_scale ___(( double ));
extern int pg_ascale ___(( double, double, double ));
extern int pg_child ___(( char * ));
extern int pg_print_gob ___(( char * ));

/* Primitive routines */
extern int pg_cylinder ___(( void ));
extern int pg_sphere ___(( void ));
extern int pg_torus ___(( double, double ));
extern int pg_polymarker ___(( P_Vlist * ));
extern int pg_polyline ___(( P_Vlist * ));
extern int pg_polygon ___(( P_Vlist * ));
extern int pg_tristrip ___(( P_Vlist * ));
extern int pg_mesh ___(( P_Vlist *, int *, int *, int ));
extern int pg_bezier ___(( P_Vlist * )); /* always 16 */
extern int pg_text ___(( char *, P_Point *, P_Vector *, P_Vector * ));
extern int pg_light ___(( P_Point *, P_Color * ));
extern int pg_ambient ___(( P_Color * ));

/* Composite routines */
extern int pg_axis ___(( P_Point *, P_Point *, P_Vector *, double, 
                        double, int, char *, double, int ));
extern int pg_boundbox ___(( P_Point *, P_Point * ));
extern int pg_isosurface ___(( int type, float *data, float *valdata,
		  int nx, int ny, int nz, double value,
		  P_Point *corner1, P_Point *corner2,
		  int show_inside, int ftn_order ));
extern int pg_zsurface ___(( int, float *, float *, 
                  int, int, P_Point *, P_Point *, 
                  void (*)(int *, float *, int *, int * ), int ));
extern int pg_rand_zsurf ___(( P_Vlist *vlist, 
			 void (*)(int *, float *, float *, float *, int *) ));
extern int pg_rand_isosurf ___((P_Vlist *vlist, double value, 
				int show_inside));

/* Camera */
extern int pg_camera ___(( char *, P_Point *, P_Point *,
                     P_Vector *, double, double, double ));
extern int pg_print_camera ___(( char * ));
extern int pg_camera_background ___(( char *, P_Color * ));

/* Snap */
extern int pg_snap ___(( char *, char *, char * ));

/* Color map */
extern int pg_set_cmap ___((double, double, void (*)( float *,
						float *, float *,
						float *, float * ) ));

extern int pg_std_cmap ___((double, double, int));

extern int pg_cmap_color ___((float *, float *, float *, float *, float *));

#endif /* __cplusplus */

/* Clean up the prototyping macros */
#undef __
#undef ___



/* $Log$
 * Revision 1.2  1993/12/30 13:32:55  etseidl
 * mostly rcs id stuff
 *
 * Revision 1.3  1992/05/26  20:25:21  jannsen
 * make derivative bounds checking optional
 * add code to allow bound intermediates computable in shell blocks
 *
 * Revision 1.2  1992/05/13  18:29:42  jannsen
 * added bounds checking for derivative integrals
 *
 * Revision 1.1.1.1  1992/03/17  16:33:13  seidl
 * DOE-NIH Quantum Chemistry Library 0.0
 *
 * Revision 1.1  1992/03/17  16:33:12  seidl
 * Initial revision
 *
 * Revision 1.1  1992/02/03  15:49:49  seidl
 * Initial revision
 *
 * Revision 1.2  1992/01/10  18:01:37  cljanss
 * added globals needed to deal with integral storage
 *
 * Revision 1.1  1991/12/14  00:19:36  cljanss
 * Initial revision
 * */
/*
 * Revision 1.5  91/11/22  17:49:25  cljanss
 * bound matrix generated is now handled by a separate flag
 * 
 * Revision 1.4  91/09/28  18:17:19  cljanss
 * Intermediates are no longer stored, if requested with flags.
 * 
 * Revision 1.3  91/06/19  23:18:19  janssen
 * added computation and checking of upper bounds to shell quartet integrals
 * 
 * Revision 1.2  1991/06/19  15:19:38  janssen
 * First stab at two electron derivative integrals--might work might not
 *
 * Revision 1.1  1991/06/16  16:40:07  janssen
 * Initial revision
 * */

#ifndef _INTERMEDIATES_H

#ifdef ALLOC_INTERMEDIATES
# define EXTERN
#define INITIALIZE(x,y) x=y
#else
# define EXTERN extern
#define INITIALIZE(x,y) x
#endif

EXTERN int int_derivative_bounds;
EXTERN int int_storage_threshold;
EXTERN int INITIALIZE(int_maxsize,0);
EXTERN int int_integral_storage;
EXTERN int int_used_integral_storage;
EXTERN int int_store1;
EXTERN int int_store2;
EXTERN int int_store_bounds;
EXTERN int_vector_t int_shell_to_prim;
EXTERN doublep_vector_t int_shell_r;
EXTERN double_matrix_t int_shell_Q;
EXTERN double_matrix_t int_prim_zeta;
EXTERN double_matrix_t int_prim_k;
EXTERN double_matrix_t int_prim_oo2zeta;
EXTERN double_array3_t int_prim_p;
EXTERN doublep_array4_t *int_con_ints;
EXTERN doublep_array4_t ****int_con_ints_array;  /* The contr. int. inter. */

EXTERN centers_t *int_cs1;
EXTERN centers_t *int_cs2;
EXTERN centers_t *int_cs3;
EXTERN centers_t *int_cs4;

EXTERN shell_t *int_shell1;
EXTERN shell_t *int_shell2;
EXTERN shell_t *int_shell3;
EXTERN shell_t *int_shell4;

EXTERN double *int_buffer;
EXTERN double *int_derint_buffer;

EXTERN int int_expweight1; /* For exponent weighted contractions. */
EXTERN int int_expweight2; /* For exponent weighted contractions. */
EXTERN int int_expweight3; /* For exponent weighted contractions. */
EXTERN int int_expweight4; /* For exponent weighted contractions. */

#undef EXTERN
#define _INTERMEDIATES_H
#endif /* _INTERMEDIATES_H */

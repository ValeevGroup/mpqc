
/* $Log$
 * Revision 1.1  1993/12/29 12:53:01  etseidl
 * Initial revision
 *
 * Revision 1.1.1.1  1992/03/17  16:32:32  seidl
 * DOE-NIH Quantum Chemistry Library 0.0
 *
 * Revision 1.1  1992/03/17  16:32:31  seidl
 * Initial revision
 *
 * Revision 1.1  1992/02/03  15:49:49  seidl
 * Initial revision
 *
 * Revision 1.1  1991/12/14  00:19:36  cljanss
 * Initial revision
 * */
/*
 * Revision 1.2  91/09/28  18:17:13  cljanss
 * Intermediates are no longer stored, if requested with flags.
 * 
 * Revision 1.1  91/06/16  16:40:07  janssen
 * Initial revision
 *  */

#ifdef ALLOC_BUILDINTER
# define EXTERN
#else
# define EXTERN extern
#endif

EXTERN double int_v_ooze;
EXTERN double int_v_zeta12;
EXTERN double int_v_zeta34;
EXTERN double int_v_oo2zeta12;
EXTERN double int_v_oo2zeta34;
EXTERN double int_v_W0;
EXTERN double int_v_W1;
EXTERN double int_v_W2;
EXTERN double int_v_p120;
EXTERN double int_v_p121;
EXTERN double int_v_p122;
EXTERN double int_v_p340;
EXTERN double int_v_p341;
EXTERN double int_v_p342;
EXTERN double int_v_r10;
EXTERN double int_v_r11;
EXTERN double int_v_r12;
EXTERN double int_v_r20;
EXTERN double int_v_r21;
EXTERN double int_v_r22;
EXTERN double int_v_r30;
EXTERN double int_v_r31;
EXTERN double int_v_r32;
EXTERN double int_v_r40;
EXTERN double int_v_r41;
EXTERN double int_v_r42;
EXTERN double int_v_k12;
EXTERN double int_v_k34;
EXTERN doublep_array3_t int_v_list;

#undef EXTERN


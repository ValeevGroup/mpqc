
/* $Log$
 * Revision 1.2  1993/12/30 13:32:49  etseidl
 * mostly rcs id stuff
 *
 * Revision 1.1.1.1  1992/03/17  16:32:47  seidl
 * DOE-NIH Quantum Chemistry Library 0.0
 *
 * Revision 1.1  1992/03/17  16:32:46  seidl
 * Initial revision
 *
 * Revision 1.1  1992/02/03  15:49:49  seidl
 * Initial revision
 *
 * Revision 1.1  1991/12/14  00:19:36  cljanss
 * Initial revision
 * */
/*
 * Revision 1.1  91/06/16  16:40:07  janssen
 * Initial revision
 *  */

#ifdef ALLOC_FJTTABLE
# define EXTERN
#else
# define EXTERN extern
#endif

EXTERN double_vector_t int_fjttable;

#undef EXTERN


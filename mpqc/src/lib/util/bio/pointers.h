
/* $Id$ */
/* $Log$
 * Revision 1.1  1993/12/29 12:53:41  etseidl
 * Initial revision
 *
 * Revision 1.2  1992/03/30  22:44:17  jannsen
 * Merged in Sandia non CVS codes.
 *
 * Revision 1.3  91/11/18  17:34:38  cljanss
 * put the bio_error FILE pointer in.
 * 
 * Revision 1.2  91/09/28  20:50:04  cljanss
 * changed the global name ptr to bio_ptr
 * 
 * Revision 1.1  1991/06/16  20:28:50  seidl
 * Initial revision
 * */

#ifdef ALLOC_GLOBALS
#define EXTERN
#else
#define EXTERN extern
#endif

struct pointer {
     int *wptr;
     };

EXTERN struct pointer bio_ptr;

#ifdef ALLOC_GLOBALS
EXTERN FILE *bio_error=stderr;
#else
EXTERN FILE *bio_error;
#endif

#undef EXTERN
#undef INITIALIZE


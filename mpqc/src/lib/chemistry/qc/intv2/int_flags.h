
/* $Log$
 * Revision 1.1  1993/12/29 12:53:02  etseidl
 * Initial revision
 *
 * Revision 1.2  1992/05/26  20:25:19  jannsen
 * make derivative bounds checking optional
 * add code to allow bound intermediates computable in shell blocks
 *
 * Revision 1.1.1.1  1992/03/17  16:33:01  seidl
 * DOE-NIH Quantum Chemistry Library 0.0
 *
 * Revision 1.1  1992/03/17  16:33:00  seidl
 * Initial revision
 *
 * Revision 1.1  1992/02/03  15:49:49  seidl
 * Initial revision
 *
 * Revision 1.1  1991/12/14  00:19:36  cljanss
 * Initial revision
 * */
/*
 * Revision 1.5  91/11/22  17:49:24  cljanss
 * bound matrix generated is now handled by a separate flag
 * 
 * Revision 1.4  91/09/28  18:17:18  cljanss
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

#define INT_EREP     1
#define INT_REDUND   2  /* Do not remove redundant integrals. */
#define INT_NOPERM   4  /* Do not return integrals with permuted the centers. */
#define INT_NOBCHK   8  /* Do not check the shell quartet's upper bound. */
#define INT_NOSTR1  16 /* Do not store the O(n) intermediates, recompute them as needed. */
#define INT_NOSTR2  32 /* Do not store the O(n^2) intermediates, recompute them as needed. */
#define INT_NOSTRB  64 /* Do not store the bounds checking array. */
#define INT_NODERB 128 /* Do not prepare the routines for derivative bounds. */

#define INT_TYPE_SP 1 /* This is an sp shell. */


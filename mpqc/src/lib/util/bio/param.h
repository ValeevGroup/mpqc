
/* $Id$ */
/* $Log$
 * Revision 1.1  1993/12/29 12:53:41  etseidl
 * Initial revision
 *
 * Revision 1.3  1992/04/06  12:28:18  seidl
 * merge in sandia changes
 *
 * Revision 1.2  1992/03/30  22:43:58  jannsen
 * Merged in Sandia non CVS codes.
 *
 * Revision 1.2  91/06/17  16:07:06  seidl
 * add defs of BIO_KEEP and BIO_UNLINK
 * 
 * Revision 1.1  1991/06/16  20:28:50  seidl
 * Initial revision
 * */

#define BIO_KEEP 3
#define BIO_UNLINK 4

#define MAX_UNIT 129
#define MAX_VOLUME 8
#define MAX_STRING 512

#define BIO_UNINITED 0
#define BIO_UNOPENED 1
#define BIO_NOMETHOD 0
#define BIO_SEQUENTIAL 1
#define BIO_S_ASYNC 2
#define BIO_R_ASYNC 3
#define BIO_RAM 4
#define BIO_FORTRAN 5

#define IOOP_READ 1
#define IOOP_WRITE 2


#if (defined(DEC)||defined(SUN))
#define SEEK_SET 0
#define SEEK_CUR 1
#define SEEK_END 2
#endif


/* $Log$
 * Revision 1.1  1993/12/29 12:53:39  etseidl
 * Initial revision
 *
 * Revision 1.4  1992/07/20  18:35:40  seidl
 * add code to make sure a string is non-null
 *
 * Revision 1.3  1992/06/17  22:16:46  jannsen
 * cleaned up for saber-c
 *
 * Revision 1.2  1992/03/30  23:16:28  seidl
 * Merged in Sandia non CVS codes.
 *
 * Revision 1.2  91/12/02  17:34:56  cljanss
 * clean up use of void for old compilers
 * 
 * Revision 1.1  91/11/18  18:16:26  cljanss
 * Initial revision
 *  */
static char *rcsid = "$Id$";

#define NO_TEMPLATES
#include <stdio.h>
#include <stdlib.h>
#include <tmpl.h>
#include <util/bio/libbio.h>
#include "sgen.h"

/* bread_boolean.c,v
 * Revision 1.2  91/09/30  13:49:13  cljanss
 * renamed the global ptr to bio_ptr.
 * 
 * Revision 1.1  1991/06/17  15:57:19  seidl
 * Initial revision
 * */

VOID
bread_boolean(unit,buff,offset,size)
int unit;
int *buff;
int *offset;
int size;
{
  int fword = *offset;

  biordr_(&unit,buff,&fword,&size);
  *offset = fword+size;
  bio_ptr.wptr[unit] = *offset;
  }

/* bread_char.c,v
 * Revision 1.2  91/09/30  13:49:27  cljanss
 * renamed the global ptr to bio_ptr.
 * 
 * Revision 1.1  1991/06/17  15:57:19  seidl
 * Initial revision
 * */

VOID
bread_char(unit,buff,offset,size)
int unit;
char *buff;
int *offset;
int size;
{
  int fword = *offset;

  biordr_(&unit,buff,&fword,&size);
  *offset = fword+size;
  bio_ptr.wptr[unit] = *offset;
  }

/* bread_double.c,v
 * Revision 1.2  91/09/30  13:49:29  cljanss
 * renamed the global ptr to bio_ptr.
 * 
 * Revision 1.1  1991/06/17  15:57:19  seidl
 * Initial revision
 * */

VOID
bread_double(unit,buff,offset,size)
int unit;
double *buff;
int *offset;
int size;
{
  int fword = *offset;

  biordr_(&unit,buff,&fword,&size);
  *offset = fword+size;
  bio_ptr.wptr[unit] = *offset;
  }

/* bread_float.c,v
 * Revision 1.2  91/09/30  13:49:30  cljanss
 * renamed the global ptr to bio_ptr.
 * 
 * Revision 1.1  1991/06/17  15:57:19  seidl
 * Initial revision
 * */

VOID
bread_float(unit,buff,offset,size)
int unit;
float *buff;
int *offset;
int size;
{
  int fword = *offset;

  biordr_(&unit,buff,&fword,&size);
  *offset = fword+size;
  bio_ptr.wptr[unit] = *offset;
  }

/* bread_int.c,v
 * Revision 1.2  91/09/30  13:49:31  cljanss
 * renamed the global ptr to bio_ptr.
 * 
 * Revision 1.1  1991/06/17  15:57:19  seidl
 * Initial revision
 * */

VOID
bread_int(unit,buff,offset,size)
int unit;
int *buff;
int *offset;
int size;
{
  int fword = *offset;

  biordr_(&unit,buff,&fword,&size);
  *offset = fword+size;
  bio_ptr.wptr[unit] = *offset;
  }

/* bread_long.c,v
 * Revision 1.2  91/09/30  13:49:32  cljanss
 * renamed the global ptr to bio_ptr.
 * 
 * Revision 1.1  1991/06/17  15:57:19  seidl
 * Initial revision
 * */

VOID
bread_long(unit,buff,offset,size)
int unit;
long *buff;
int *offset;
int size;
{
  int fword = *offset;

  biordr_(&unit,buff,&fword,&size);
  *offset = fword+size;
  bio_ptr.wptr[unit] = *offset;
  }

/* bread_pointer.c,v
 * Revision 1.4  91/09/30  13:49:33  cljanss
 * renamed the global ptr to bio_ptr.
 * 
 * Revision 1.3  1991/07/19  20:55:39  cljanss
 * changed name to be compatible with the new version of sgen.
 *
 * Revision 1.2  1991/07/19  16:17:30  cljanss
 * These are Ed's latest changes from CCQC.
 *
 * Revision 1.1  1991/06/19  23:06:50  seidl
 * Initial revision
 * */

VOID
bread_pointer(unit,buff,offset,size)
int unit;
VOID_PTR buff;
int *offset;
int size;
{
  int fword = *offset;

  biordr_(&unit,buff,&fword,&size);
  *offset = fword+size;
  bio_ptr.wptr[unit] = *offset;
  }

int *
bread_test_pointer(unit,offset,size)
int unit;
int *offset;
int size;
{
  int *buff;
  int fword = *offset;

  biordr_(&unit,&buff,&fword,&size);
  *offset = fword+size;
  bio_ptr.wptr[unit] = *offset;

  return buff;
  }

/* bread_string.c,v
 * Revision 1.3  91/09/30  13:49:34  cljanss
 * renamed the global ptr to bio_ptr.
 * 
 * Revision 1.2  1991/07/19  20:44:29  cljanss
 * Fixed substantial bugs in the binary string read and write routines.
 *
 * Revision 1.1  1991/06/17  15:57:19  seidl
 * Initial revision
 * */

VOID
bread_string(unit,buff,offset,size)
int unit;
char **buff;
int *offset;
int size;
{
  int ilength = sizeof(int);
  int length;
  int fword;
  int i;

  for (i=0; i<size; i+=sizeof(*buff)) {
    /* Read in the length of the string. */
    fword = *offset;
    biordr_(&unit,&length,&fword,&ilength);
    *offset += ilength;
    bio_ptr.wptr[unit] = *offset;

    /* Allocate storage for the string. */
    if(length) {
      *buff = (char *) malloc(length);

      /* Read in the string. */
      fword = *offset;
      biordr_(&unit,*buff,&fword,&length);
      *offset += length;
      bio_ptr.wptr[unit] = *offset;
      }
    buff++;
    }
  }


/* $Log$
 * Revision 1.1  1993/12/29 12:53:39  etseidl
 * Initial revision
 *
 * Revision 1.4  1992/07/20  18:35:44  seidl
 * add code to make sure a string is non-null
 *
 * Revision 1.3  1992/06/17  22:16:48  jannsen
 * cleaned up for saber-c
 *
 * Revision 1.2  1992/03/30  23:16:33  seidl
 * Merged in Sandia non CVS codes.
 *
 * Revision 1.2  91/12/02  17:35:08  cljanss
 * clean up use of void for old compilers
 * 
 * Revision 1.1  91/11/18  18:16:27  cljanss
 * Initial revision
 *  */
static char *rcsid = "$Id$";

#define NO_TEMPLATES
#include <stdio.h>
#include <tmpl.h>
#include <util/bio/libbio.h>
#include "sgen.h"

/* bwrite_boolean.c,v
 * Revision 1.2  91/09/30  13:49:35  cljanss
 * renamed the global ptr to bio_ptr.
 * 
 * Revision 1.1  1991/06/17  15:57:19  seidl
 * Initial revision
 * */

VOID
bwrite_boolean(unit,buff,offset,size)
int unit;
int *buff;
int *offset;
int size;
{
  int fword = *offset;

  biowrr_(&unit,buff,&fword,&size);
  *offset = fword+size;
  bio_ptr.wptr[unit] = *offset;
  }

/* bwrite_char.c,v
 * Revision 1.2  91/09/30  13:49:36  cljanss
 * renamed the global ptr to bio_ptr.
 * 
 * Revision 1.1  1991/06/17  15:57:19  seidl
 * Initial revision
 * */

VOID
bwrite_char(unit,buff,offset,size)
int unit;
char *buff;
int *offset;
int size;
{
  int fword = *offset;

  biowrr_(&unit,buff,&fword,&size);
  *offset = fword+size;
  bio_ptr.wptr[unit] = *offset;
  }

/* bwrite_double.c,v
 * Revision 1.2  91/09/30  13:49:37  cljanss
 * renamed the global ptr to bio_ptr.
 * 
 * Revision 1.1  1991/06/17  15:57:19  seidl
 * Initial revision
 * */

VOID
bwrite_double(unit,buff,offset,size)
int unit;
double *buff;
int *offset;
int size;
{
  int fword = *offset;

  biowrr_(&unit,buff,&fword,&size);
  *offset = fword+size;
  bio_ptr.wptr[unit] = *offset;
  }

/* bwrite_float.c,v
 * Revision 1.2  91/09/30  13:49:38  cljanss
 * renamed the global ptr to bio_ptr.
 * 
 * Revision 1.1  1991/06/17  15:57:19  seidl
 * Initial revision
 * */

VOID
bwrite_float(unit,buff,offset,size)
int unit;
float *buff;
int *offset;
int size;
{
  int fword = *offset;

  biowrr_(&unit,buff,&fword,&size);
  *offset = fword+size;
  bio_ptr.wptr[unit] = *offset;
  }

/* bwrite_int.c,v
 * Revision 1.2  91/09/30  13:49:39  cljanss
 * renamed the global ptr to bio_ptr.
 * 
 * Revision 1.1  1991/06/17  15:57:19  seidl
 * Initial revision
 * */

VOID
bwrite_int(unit,buff,offset,size)
int unit;
int *buff;
int *offset;
int size;
{
  int fword = *offset;

  biowrr_(&unit,buff,&fword,&size);
  *offset = fword+size;
  bio_ptr.wptr[unit] = *offset;
  }

/* bwrite_long.c,v
 * Revision 1.2  91/09/30  13:49:39  cljanss
 * renamed the global ptr to bio_ptr.
 * 
 * Revision 1.1  1991/06/17  15:57:19  seidl
 * Initial revision
 * */

VOID
bwrite_long(unit,buff,offset,size)
int unit;
long *buff;
int *offset;
int size;
{
  int fword = *offset;

  biowrr_(&unit,buff,&fword,&size);
  *offset = fword+size;
  bio_ptr.wptr[unit] = *offset;
  }

/* bwrite_pointer.c,v
 * Revision 1.3  91/09/30  13:49:40  cljanss
 * renamed the global ptr to bio_ptr.
 * 
 * Revision 1.2  1991/07/19  16:17:30  cljanss
 * These are Ed's latest changes from CCQC.
 *
 * Revision 1.1  1991/06/19  23:06:57  seidl
 * Initial revision
 * */

VOID
bwrite_pointer(unit,buff,offset,size)
int unit;
VOID_PTR buff;
int *offset;
int size;
{
  int fword = *offset;

  biowrr_(&unit,buff,&fword,&size);

  *offset = fword+size;
  bio_ptr.wptr[unit] = *offset;
  }

/* bwrite_string.c,v
 * Revision 1.3  91/09/30  13:49:41  cljanss
 * renamed the global ptr to bio_ptr.
 * 
 * Revision 1.2  1991/07/19  20:44:29  cljanss
 * Fixed substantial bugs in the binary string read and write routines.
 *
 * Revision 1.1  1991/06/17  15:57:19  seidl
 * Initial revision
 * */

VOID
bwrite_string(unit,buff,offset,size)
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
    /* Write out the length of the string. */
    if(*buff==NULL) length=0;
    else length = strlen(*buff)+1;
    fword = *offset;
    biowrr_(&unit,&length,&fword,&ilength);
    *offset += ilength;
    bio_ptr.wptr[unit] = *offset;

    /* Write out the string. */
    if(length) {
      fword = *offset;
      biowrr_(&unit,*buff,&fword,&length);
      *offset += length;
      bio_ptr.wptr[unit] = *offset;
      }
    buff++;
    }
  }


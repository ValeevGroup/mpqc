
/* $Log$
 * Revision 1.1  1993/12/29 12:53:39  etseidl
 * Initial revision
 *
 * Revision 1.4  1992/07/20  18:35:46  seidl
 * add code to make sure a string is non-null
 *
 * Revision 1.3  1992/06/17  22:16:54  jannsen
 * cleaned up for saber-c
 *
 * Revision 1.2  1992/03/30  23:16:42  seidl
 * Merged in Sandia non CVS codes.
 *
 * Revision 1.2  91/12/02  17:35:44  cljanss
 * clean up use of void for old compilers
 * 
 * Revision 1.1  91/11/18  18:16:29  cljanss
 * Initial revision
 *  */
static char *rcsid = "$Id$";

#define NO_TEMPLATES
#include <stdio.h>
#include <tmpl.h>
#include <util/bio/libbio.h>
#include "sgen.h"

VOID
fwrite_boolean(fp,buff,offset,size)
FILE *fp;
int *buff;
int *offset;
int size;
{
  int fword = *offset;

  fseek(fp,(long)*offset,SEEK_SET);
  *offset = fword+fwrite(buff,1,size,fp);
  }

VOID
fwrite_char(fp,buff,offset,size)
FILE *fp;
char *buff;
int *offset;
int size;
{
  int fword = *offset;

  fseek(fp,(long)*offset,SEEK_SET);
  *offset = fword+fwrite(buff,1,size,fp);
  }

VOID
fwrite_double(fp,buff,offset,size)
FILE *fp;
double *buff;
int *offset;
int size;
{
  int fword = *offset;

  fseek(fp,(long)*offset,SEEK_SET);
  *offset = fword+fwrite(buff,1,size,fp);
  }

VOID
fwrite_float(fp,buff,offset,size)
FILE *fp;
float *buff;
int *offset;
int size;
{
  int fword = *offset;

  fseek(fp,(long)*offset,SEEK_SET);
  *offset = fword+fwrite(buff,1,size,fp);
  }

VOID
fwrite_int(fp,buff,offset,size)
FILE *fp;
int *buff;
int *offset;
int size;
{
  int fword = *offset;

  fseek(fp,(long)*offset,SEEK_SET);
  *offset = fword+fwrite(buff,1,size,fp);
  }

VOID
fwrite_long(fp,buff,offset,size)
FILE *fp;
long *buff;
int *offset;
int size;
{
  int fword = *offset;

  fseek(fp,(long)*offset,SEEK_SET);
  *offset = fword+fwrite(buff,1,size,fp);
  }

VOID
fwrite_pointer(fp,buff,offset,size)
FILE *fp;
void *buff;
int *offset;
int size;
{
  int fword = *offset;


  fseek(fp,(long)*offset,SEEK_SET);
  *offset = fword+fwrite(buff,1,size,fp);
  }

VOID
fwrite_string(fp,buff,offset,size)
FILE *fp;
char **buff;
int *offset;
int size;
{
  CONST int ilength = sizeof(int);
  int length;
  int i;

  fseek(fp,(long)*offset,SEEK_SET);
  for (i=0; i<size; i+=sizeof(*buff)) {
    /* Write out the length of the string. */
    if(*buff==NULL) length=0;
    else length = strlen(*buff)+1;
    *offset += fwrite(&length,1,ilength,fp);

    /* Write out the string. */
    if(length) *offset += fwrite(*buff,1,length,fp);
    buff++;
    }
  }

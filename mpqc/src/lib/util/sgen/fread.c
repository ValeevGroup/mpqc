
/* $Log$
 * Revision 1.1  1993/12/29 12:53:39  etseidl
 * Initial revision
 *
 * Revision 1.4  1992/07/20  18:35:45  seidl
 * add code to make sure a string is non-null
 *
 * Revision 1.3  1992/06/17  22:16:52  jannsen
 * cleaned up for saber-c
 *
 * Revision 1.2  1992/03/30  23:16:37  seidl
 * Merged in Sandia non CVS codes.
 *
 * Revision 1.2  91/12/02  17:35:36  cljanss
 * clean up use of void for old compilers
 * 
 * Revision 1.1  91/11/18  18:16:28  cljanss
 * Initial revision
 *  */
static char *rcsid = "$Id$";

#define NO_TEMPLATES
#include <stdio.h>
#include <stdlib.h>
#include <tmpl.h>
#include <util/bio/libbio.h>
#include "sgen.h"

VOID
fread_boolean(fp,buff,offset,size)
FILE *fp;
int *buff;
int *offset;
int size;
{
  int fword = *offset;

  fseek(fp,(long)*offset,SEEK_SET);
  *offset = fword + fread(buff,1,size,fp);
  }

VOID
fread_char(fp,buff,offset,size)
FILE *fp;
char *buff;
int *offset;
int size;
{
  int fword = *offset;

  fseek(fp,(long)*offset,SEEK_SET);
  *offset = fword+fread(buff,1,size,fp);
  }

VOID
fread_double(fp,buff,offset,size)
FILE *fp;
double *buff;
int *offset;
int size;
{
  int fword = *offset;

  fseek(fp,(long)*offset,SEEK_SET);
  *offset = fword+fread(buff,1,size,fp);
  }

VOID
fread_float(fp,buff,offset,size)
FILE *fp;
float *buff;
int *offset;
int size;
{
  int fword = *offset;

  fseek(fp,(long)*offset,SEEK_SET);
  *offset = fword+fread(buff,1,size,fp);
  }

VOID
fread_int(fp,buff,offset,size)
FILE *fp;
int *buff;
int *offset;
int size;
{
  int fword = *offset;

  fseek(fp,(long)*offset,SEEK_SET);
  *offset = fword+fread(buff,1,size,fp);
  }

VOID
fread_long(fp,buff,offset,size)
FILE *fp;
long *buff;
int *offset;
int size;
{
  int fword = *offset;

  fseek(fp,(long)*offset,SEEK_SET);
  *offset = fword+fread(buff,1,size,fp);
  }

VOID
fread_pointer(fp,buff,offset,size)
FILE *fp;
VOID_PTR buff;
int *offset;
int size;
{
  int fword = *offset;

  fseek(fp,(long)*offset,SEEK_SET);
  *offset = fword+fread(buff,1,size,fp);
  }

int *
fread_test_pointer(fp,offset,size)
FILE *fp;
int *offset;
int size;
{
  int *buff;
  int fword = *offset;

  fseek(fp,(long)*offset,SEEK_SET);
  *offset = fword+fread(&buff,1,size,fp);

  return buff;
  }

VOID
fread_string(fp,buff,offset,size)
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
    /* Read in the length of the string. */
    *offset += fread(&length,1,ilength,fp);

    /* Allocate storage for the string. */
    if(length) {
      *buff = (char *) malloc(length);

      /* Read in the string. */
      *offset += fread(*buff,1,length,fp);
      }
    buff++;
    }
  }

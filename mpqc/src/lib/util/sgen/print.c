
/* $Log$
 * Revision 1.2  1994/08/25 22:48:23  etseidl
 * remove rcsids and fix some warnings
 *
 * Revision 1.1.1.1  1993/12/29  12:53:41  etseidl
 * SC source tree 0.1
 *
 * Revision 1.4  1992/10/09  19:34:02  seidl
 * print strings in quotes
 *
 * Revision 1.3  1992/06/17  22:17:03  jannsen
 * cleaned up for saber-c
 *
 * Revision 1.2  1992/03/30  23:16:51  seidl
 * Merged in Sandia non CVS codes.
 *
 * Revision 1.1  91/11/18  18:17:28  cljanss
 * Initial revision
 *  */

#include <stdio.h>
#include "sgen.h"

int sgen_print_nindent=0;

static int suppress_indent = 0;

/* print.c,v
 * Revision 1.1  91/06/16  02:52:39  janssen
 * Initial revision
 *  */

void
sgen_print_suppress_indent()
{
  suppress_indent = 1;
  }

void
sgen_print_indent(fp)
FILE *fp;
{
  int i;

  if (suppress_indent) {
    suppress_indent = 0;
    return;
    }

  for (i=0; i<sgen_print_nindent; i++) {
    fprintf(fp," ");
    }
  }


/* print_boolean.c,v
 * Revision 1.1  91/06/16  02:52:39  janssen
 * Initial revision
 *  */

void
print_boolean(fp,value)
FILE *fp;
int *value;
{
  if (*value) fprintf(fp," yes");
  else fprintf(fp," no");
  }


/* print_char.c,v
 * Revision 1.1  91/06/16  02:52:39  janssen
 * Initial revision
 *  */

void
print_char(fp,value)
FILE *fp;
char *value;
{
  fprintf(fp," %c",*value);
  }


/* print_double.c,v
 * Revision 1.2  91/09/30  13:49:59  cljanss
 * the format of the double datum can now be adjusted.
 * 
 * Revision 1.1  1991/06/16  02:52:39  janssen
 * Initial revision
 * */

static char *double_fmt=" %lf";

void
print_double(fp,value)
FILE *fp;
double *value;
{
  fprintf(fp,double_fmt,*value);
  }

void
print_double_set(fmt)
char *fmt;
{
  double_fmt = fmt;
  }


/* print_float.c,v
 * Revision 1.1  91/06/16  02:52:39  janssen
 * Initial revision
 *  */

void
print_float(fp,value)
FILE *fp;
float *value;
{
  fprintf(fp," %f",*value);
  }


/* print_int.c,v
 * Revision 1.1  91/06/16  02:52:39  janssen
 * Initial revision
 *  */

void
print_int(fp,value)
FILE *fp;
int *value;
{
  fprintf(fp," %d",*value);
  }


/* print_long.c,v
 * Revision 1.1  91/06/16  02:52:39  janssen
 * Initial revision
 *  */

void
print_long(fp,value)
FILE *fp;
long *value;
{
  fprintf(fp,"%ld",*value);
  }


/* print_string.c,v
 * Revision 1.1  91/06/16  02:52:39  janssen
 * Initial revision
 *  */

void
print_string(fp,value)
FILE *fp;
char **value;
{
  fprintf(fp," \"%s\"",*value);
  }


/* print_unsigned_int.c,v
 * Revision 1.1  91/06/16  02:52:39  janssen
 * Initial revision
 *  */

void
print_unsigned(fp,value)
FILE *fp;
unsigned int *value;
{
  fprintf(fp,"%u",*value);
  }


/* print_unsigned_long.c,v
 * Revision 1.1  91/06/16  02:52:39  janssen
 * Initial revision
 *  */

void
print_unsigned_long(fp,value)
FILE *fp;
unsigned long *value;
{
  fprintf(fp,"%lu",*value);
  }


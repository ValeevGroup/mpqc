
/* $Log$
 * Revision 1.1  1993/12/29 12:53:40  etseidl
 * Initial revision
 *
 * Revision 1.3  1992/06/17  22:16:55  jannsen
 * cleaned up for saber-c
 *
 * Revision 1.2  1992/03/30  23:16:46  seidl
 * Merged in Sandia non CVS codes.
 *
 * Revision 1.1  91/11/18  18:16:29  cljanss
 * Initial revision
 *  */
static char *rcsid = "$Id$";

#include <stdio.h>
#include <util/ipv2/ip_libv2.h>
#include "sgen.h"

/* ip_read_boolean_v.c,v
 * Revision 1.1  91/06/16  02:52:39  janssen
 * Initial revision
 *  */

int
ip_read_boolean_v(keyword,value,n,v)
char *keyword;
int *value;
int n;
int *v;
{
  return ip_boolean_v(keyword,value,n,v);
  }


/* ip_read_char_v.c,v
 * Revision 1.1  91/06/16  02:52:39  janssen
 * Initial revision
 *  */

int
ip_read_char_v(keyword,value,n,v)
char *keyword;
char *value;
int n;
int *v;
{
  char *val;
  int errcod;

  if ((errcod = ip_value_v(keyword,&val,n,v))!=IPE_OK) return errcod;

  if (sscanf(val,"%c",value) != 1) return IPE_TYPE;

  return IPE_OK;
  }


/* ip_read_double_v.c,v
 * Revision 1.1  91/06/16  02:52:39  janssen
 * Initial revision
 *  */

int
ip_read_double_v(keyword,value,n,v)
char *keyword;
double *value;
int n;
int *v;
{
  char *val;
  int errcod;

  if ((errcod = ip_value_v(keyword,&val,n,v))!=IPE_OK) return errcod;

  if (sscanf(val,"%lf",value) != 1) return IPE_TYPE;

  return IPE_OK;
  }


/* ip_read_float_v.c,v
 * Revision 1.1  91/06/16  02:52:39  janssen
 * Initial revision
 *  */

int
ip_read_float_v(keyword,value,n,v)
char *keyword;
float *value;
int n;
int *v;
{
  char *val;
  int errcod;

  if ((errcod = ip_value_v(keyword,&val,n,v))!=IPE_OK) return errcod;

  if (sscanf(val,"%f",value) != 1) return IPE_TYPE;

  return IPE_OK;
  }


/* ip_read_int_v.c,v
 * Revision 1.1  91/06/16  02:52:39  janssen
 * Initial revision
 *  */

int
ip_read_int_v(keyword,value,n,v)
char *keyword;
int *value;
int n;
int *v;
{
  char *val;
  int errcod;

  if ((errcod = ip_value_v(keyword,&val,n,v))!=IPE_OK) return errcod;

  if (sscanf(val,"%d",value) != 1) return IPE_TYPE;

  return IPE_OK;
  }


/* ip_read_long_v.c,v
 * Revision 1.1  91/06/16  02:52:39  janssen
 * Initial revision
 *  */

int
ip_read_long_v(keyword,value,n,v)
char *keyword;
long *value;
int n;
int *v;
{
  char *val;
  int errcod;

  if ((errcod = ip_value_v(keyword,&val,n,v))!=IPE_OK) return errcod;

  if (sscanf(val,"%ld",value) != 1) return IPE_TYPE;

  return IPE_OK;
  }


/* ip_read_string_v.c,v
 * Revision 1.1  91/06/16  02:52:39  janssen
 * Initial revision
 *  */

int
ip_read_string_v(keyword,value,n,v)
char *keyword;
char **value;
int n;
int *v;
{
  return ip_string_v(keyword,value,n,v);
  }


/* ip_read_unsigned_int_v.c,v
 * Revision 1.1  91/06/16  02:52:39  janssen
 * Initial revision
 *  */

int
ip_read_unsigned_int_v(keyword,value,n,v)
char *keyword;
unsigned int *value;
int n;
int *v;
{
  char *val;
  int errcod;

  if ((errcod = ip_value_v(keyword,&val,n,v))!=IPE_OK) return errcod;

  if (sscanf(val,"%u",value) != 1) return IPE_TYPE;

  return IPE_OK;
  }


/* ip_read_unsigned_long_v.c,v
 * Revision 1.1  91/06/16  02:52:39  janssen
 * Initial revision
 *  */

int
ip_read_unsigned_long_v(keyword,value,n,v)
char *keyword;
unsigned long *value;
int n;
int *v;
{
  char *val;
  int errcod;

  if ((errcod = ip_value_v(keyword,&val,n,v))!=IPE_OK) return errcod;

  if (sscanf(val,"%lu",value) != 1) return IPE_TYPE;

  return IPE_OK;
  }








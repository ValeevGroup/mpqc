
/* $Log$
 * Revision 1.1  1993/12/29 12:53:40  etseidl
 * Initial revision
 *
 * Revision 1.3  1992/07/20  18:35:38  seidl
 * add code to make sure a string is non-null
 *
 * Revision 1.2  1992/06/17  22:16:58  jannsen
 * cleaned up for saber-c
 *
 * Revision 1.1.1.1  1992/03/17  17:10:00  seidl
 * DOE-NIH Quantum Chemistry Library 0.0
 *
 * Revision 1.1  1992/03/17  17:09:59  seidl
 * Initial revision
 *
 * Revision 1.1  1992/01/09  12:50:06  seidl
 * Initial revision
 * */

static char rcsid[] = "$Id$";

#define NO_TEMPLATES
#include <stdio.h>
#include <math.h>
#include "sgen.h"


int
iseq_boolean(bool1,bool2)
int bool1;
int bool2;
{
  return (bool1==bool2);
  }

int
iseq_char(c1,c2)
int c1;
int c2;
{
  return(c1==c2);
  }

int
iseq_double(d1,d2)
double d1;
double d2;
{
  double tol=1.0e-15;

  return (fabs(d1-d2) < tol);
  }

int
iseq_float(f1,f2)
double f1;
double f2;
{
  float tol=1.0e-15;

  return (fabs(f1-f2) < tol);
  }

int
iseq_int(i1,i2)
int i1;
int i2;
{
  return (i1==i2);
  }

int
iseq_long(l1,l2)
long l1;
long l2;
{
  return (l1==l2);
  }

int
iseq_string(s1,s2)
char *s1;
char *s2;
{
  if(s1==s2) return(1);
  else if(s1==NULL || s2==NULL) return(0);
  return(!strcmp(s1,s2));
  }

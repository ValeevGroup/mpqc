
/* really just a copy of clj's rbcast0 routines */

/* $Log$
 * Revision 1.1  1993/12/29 12:53:40  etseidl
 * Initial revision
 *
 * Revision 1.3  1992/07/20  18:35:51  seidl
 * add code to make sure a string is non-null
 *
 * Revision 1.2  1992/06/17  22:17:07  jannsen
 * cleaned up for saber-c
 *
 * Revision 1.1.1.1  1992/03/17  17:10:09  seidl
 * DOE-NIH Quantum Chemistry Library 0.0
 *
 * Revision 1.1  1992/03/17  17:10:08  seidl
 * Initial revision
 *
 * Revision 1.1  1992/01/09  12:23:05  seidl
 * Initial revision
 * */

static char rcsid[] = "$Id$";

#define NO_TEMPLATES
#include <stdio.h>
#include <stdlib.h>
#include <tmpl.h>
#include <comm/picl/picl.h>
#include "sgen.h"

#include "sndrcv0.h"

/* rbcast0_boolean.c,v
 * Revision 1.2  91/09/30  13:50:43  cljanss
 * messed around with the way types work
 * 
 * Revision 1.1  1991/07/22  14:38:09  cljanss
 * Initial revision
 * */

void
recv0_boolean(buff,type,from,size)
int *buff;
int type;
int from;
int size;
{
  PRINT('r',TYPENOINC(),from,size);
  recv0(buff,size*sizeof(int),TYPE());
#if 0
#if defined(SUN) && defined(NIH)
  CTOHL(buff,size/sizeof(int));
#endif
#endif
  PRINT_DATA('r',"%d\n",*buff);
  }

/* rbcast0_double.c,v
 * Revision 1.2  91/09/30  13:51:11  cljanss
 * messed around with the way types work
 * 
 * Revision 1.1  1991/07/22  14:38:09  cljanss
 * Initial revision
 * */

void
recv0_double(buff,type,from,size)
double *buff;
int type;
int from;
int size;
{
  PRINT('r',TYPENOINC(),from,size);
  recv0(buff,size*sizeof(double),TYPE());
#if 0
#if defined(SUN) && defined(NIH)
  CTOHD(buff,size/sizeof(double));
#endif
#endif
  PRINT_DATA('r',"%lf\n",*buff);
  }

/* rbcast0_int.c,v
 * Revision 1.2  91/09/30  13:51:20  cljanss
 * messed around with the way types work
 * 
 * Revision 1.1  1991/07/22  14:38:09  cljanss
 * Initial revision
 * */

void
recv0_char(buff,type,from,size)
char *buff;
int type;
int from;
int size;
{
  PRINT('r',TYPENOINC(),from,size);
  recv0(buff,size,TYPE());
  PRINT_DATA('r',"%c\n",*buff);
  }

void
recv0_int(buff,type,from,size)
int *buff;
int type;
int from;
int size;
{
  PRINT('r',TYPENOINC(),from,size);
  recv0(buff,size*sizeof(int),TYPE());
#if 0
#if defined(SUN) && defined(NIH)
  CTOHL(buff,size/sizeof(int));
#endif
#endif
  PRINT_DATA('r',"%d\n",*buff);
  }

/* rbcast0_pointer.c,v
 * Revision 1.2  91/09/30  13:51:21  cljanss
 * messed around with the way types work
 * 
 * Revision 1.1  1991/07/22  14:38:09  cljanss
 * Initial revision
 * */

void
recv0_pointer(buff,type,from,size)
void *buff;
int type;
int from;
int size;
{
  PRINT('r',TYPENOINC(),from,size);
  recv0(buff,size*sizeof(void*),TYPE());
#if 0
#if defined(SUN) && defined(NIH)
  CTOHL(buff,size/sizeof(int *));
#endif
#endif
  PRINT_DATA('r',"0x%x\n",*(void**)buff);
  }

int *
recv0_test_pointer(type,from,size)
int type;
int from;
int size;
{
  int *buff;
  PRINT('r',TYPENOINC(),from,size);
  recv0(&buff,size*sizeof(void*),TYPE());
#if 0
#if defined(SUN) && defined(NIH)
  CTOHL(&buff,size/sizeof(int *));
#endif
#endif
  PRINT_DATA('r',"0x%x\n",buff);
  return buff;
  }

/* rbcast0_string.c,v
 * Revision 1.2  91/09/30  13:51:28  cljanss
 * messed around with the way types work
 * 
 * Revision 1.1  1991/07/22  14:38:09  cljanss
 * Initial revision
 * */

void
recv0_string(buff,type,from,size)
char **buff;
int type;
int from;
int size;
{
  int ilength = sizeof(int);
  int length;
  int i;

  for (i=0; i<size; i+=sizeof(*buff)) {
    /* Get the length of the string. */
    PRINT('r',TYPENOINC(),from,ilength);
    recv0(&length,ilength*sizeof(int),TYPE());
#if 0
#if defined(SUN) && defined(NIH)
  CTOHL(&length,1);
#endif
#endif
    PRINT_DATA('r',"%d\n",length);

    /* Allocate storage for the string. */
    if(length) {
      *buff = (char *) malloc(length);

     /* Read in the string. */
      PRINT('r',TYPENOINC(),from,length);
      recv0(*buff,length,TYPE());
      PRINT_DATA('r',"%s\n",*buff);
      }
    buff++;
    }

  }

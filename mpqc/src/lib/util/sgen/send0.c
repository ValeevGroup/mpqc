
/* $Log$
 * Revision 1.1  1993/12/29 12:53:40  etseidl
 * Initial revision
 *
 * Revision 1.3  1992/07/20  18:35:50  seidl
 * add code to make sure a string is non-null
 *
 * Revision 1.2  1992/06/17  22:17:11  jannsen
 * cleaned up for saber-c
 *
 * Revision 1.1.1.1  1992/03/17  17:10:16  seidl
 * DOE-NIH Quantum Chemistry Library 0.0
 *
 * Revision 1.1  1992/03/17  17:10:14  seidl
 * Initial revision
 *
 * Revision 1.1  1992/01/09  12:23:11  seidl
 * Initial revision
 * */

static char rcsid[] = "$Id$";

#define NO_TEMPLATES
#include <stdio.h>
#include <tmpl.h>
#include <comm/picl/picl.h>
#include "sgen.h"

#include "sndrcv0.h"

/* sbcast0_boolean.c,v
 * Revision 1.2  91/09/30  13:51:29  cljanss
 * messed around with the way types work
 * 
 * Revision 1.1  1991/07/22  14:38:09  cljanss
 * Initial revision
 * */

void
send0_boolean(buff,type,dest,size)
int *buff;
int type;
int dest;
int size;
{
  PRINT('s',TYPENOINC(),dest,size);
  PRINT_DATA('s',"%d\n",*buff);
#if 0
#if defined(SUN) && defined(NIH)
  HTOCL(buff,size/sizeof(int));
#endif
#endif
  send0(buff,size*sizeof(int),TYPE(),dest);
#if 0
#if defined(SUN) && defined(NIH)
  CTOHL(buff,size/sizeof(int));
#endif
#endif
  }

/* sbcast0_double.c,v
 * Revision 1.2  91/09/30  13:51:29  cljanss
 * messed around with the way types work
 * 
 * Revision 1.1  1991/07/22  14:38:09  cljanss
 * Initial revision
 * */

void
send0_double(buff,type,dest,size)
double *buff;
int type;
int dest;
int size;
{
  PRINT('s',TYPENOINC(),dest,size);
  PRINT_DATA('s',"%lf\n",*buff);
#if 0
#if defined(SUN) && defined(NIH)
  HTOCD(buff,size/sizeof(double));
#endif
#endif
  send0(buff,size*sizeof(double),TYPE(),dest);
#if 0
#if defined(SUN) && defined(NIH)
  CTOHD(buff,size/sizeof(double));
#endif
#endif
  }

/* sbcast0_int.c,v
 * Revision 1.2  91/09/30  13:51:30  cljanss
 * messed around with the way types work
 * 
 * Revision 1.1  1991/07/22  14:38:09  cljanss
 * Initial revision
 * */

void
send0_char(buff,type,dest,size)
char *buff;
int type;
int dest;
int size;
{
  PRINT('s',TYPENOINC(),dest,size);
  PRINT_DATA('s',"%c\n",*buff);
  send0(buff,size,TYPE(),dest);
  }

void
send0_int(buff,type,dest,size)
int *buff;
int type;
int dest;
int size;
{
  PRINT('s',TYPENOINC(),dest,size);
  PRINT_DATA('s',"%d\n",*buff);
#if 0
#if defined(SUN) && defined(NIH)
  HTOCL(buff,size/sizeof(int));
#endif
#endif
  send0(buff,size*sizeof(int),TYPE(),dest);
#if 0
#if defined(SUN) && defined(NIH)
  CTOHL(buff,size/sizeof(int));
#endif
#endif
  }

/* sbcast0_pointer.c,v
 * Revision 1.2  91/09/30  13:51:31  cljanss
 * messed around with the way types work
 * 
 * Revision 1.1  1991/07/22  14:38:09  cljanss
 * Initial revision
 * */

void
send0_pointer(buff,type,dest,size)
void *buff;
int type;
int dest;
int size;
{
  PRINT('s',TYPENOINC(),dest,size);
  PRINT_DATA('s',"0x%x\n",*(void**)buff);
#if 0
#if defined(SUN) && defined(NIH)
  HTOCL(buff,size/sizeof(int *));
#endif
#endif
  send0(buff,size*sizeof(void*),TYPE(),dest);
#if 0
#if defined(SUN) && defined(NIH)
  CTOHL(buff,size/sizeof(int *));
#endif
#endif
  }

/* sbcast0_string.c,v
 * Revision 1.2  91/09/30  13:51:32  cljanss
 * messed around with the way types work
 * 
 * Revision 1.1  1991/07/22  14:38:09  cljanss
 * Initial revision
 * */

void
send0_string(buff,type,dest,size)
char **buff;
int type;
int dest;
int size;
{
  int ilength = sizeof(int);
  int length;
  int i;

  for (i=0; i<size; i+=sizeof(*buff)) {
    /* Broadcast the length of the string. */
    if(*buff==NULL) length=0;
    else length = strlen(*buff)+1;
    PRINT('s',TYPENOINC(),dest,ilength);
    PRINT_DATA('s',"%d\n",length);
#if 0
#if defined(SUN) && defined(NIH)
  HTOCL(&length,1);
#endif
#endif
    send0(&length,ilength*sizeof(int),TYPE(),dest);
#if 0
#if defined(SUN) && defined(NIH)
  CTOHL(&length,1);
#endif
#endif

    /* Broadcast the string. */
    if(length) {
      PRINT('s',TYPENOINC(),dest,length);
      PRINT_DATA('s',"%s\n",*buff);
      send0(*buff,length,TYPE(),dest);
      }
    buff++;
    }

  }

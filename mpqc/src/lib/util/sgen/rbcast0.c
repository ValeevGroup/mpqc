
/* $Log$
 * Revision 1.1  1993/12/29 12:53:39  etseidl
 * Initial revision
 *
 * Revision 1.4  1992/07/20  18:35:47  seidl
 * add code to make sure a string is non-null
 *
 * Revision 1.3  1992/06/17  22:17:05  jannsen
 * cleaned up for saber-c
 *
 * Revision 1.2  1992/04/06  12:51:21  seidl
 * add _char function
 *
 * Revision 1.1.1.1  1992/03/17  17:10:07  seidl
 * DOE-NIH Quantum Chemistry Library 0.0
 *
 * Revision 1.1  1992/03/17  17:10:04  seidl
 * Initial revision
 *
 * Revision 1.2  1992/01/09  12:20:35  seidl
 * add byte swapping calls
 *
 * Revision 1.1  1991/12/20  16:10:32  seidl
 * Initial revision
 *
 * Revision 1.1  91/11/18  18:16:30  cljanss
 * Initial revision
 *  */
static char *rcsid = "$Id$";

#define NO_TEMPLATES
#include <stdio.h>
#include <stdlib.h>
#include <tmpl.h>
#include <comm/picl/picl.h>
#include "sgen.h"

#include "bcast0.h"

/* rbcast0_boolean.c,v
 * Revision 1.2  91/09/30  13:50:43  cljanss
 * messed around with the way types work
 * 
 * Revision 1.1  1991/07/22  14:38:09  cljanss
 * Initial revision
 * */

void
rbcast0_boolean(buff,type,root,size)
int *buff;
int type;
int root;
int size;
{
  PRINT('r',TYPENOINC(),root,size);
  bcast0(buff,size,TYPE(),root);
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
rbcast0_double(buff,type,root,size)
double *buff;
int type;
int root;
int size;
{
  PRINT('r',TYPENOINC(),root,size);
  bcast0(buff,size,TYPE(),root);
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
rbcast0_char(buff,type,root,size)
char *buff;
int type;
int root;
int size;
{
  PRINT('r',TYPENOINC(),root,size);
  bcast0(buff,size,TYPE(),root);
  PRINT_DATA('r',"%d\n",*buff);
  }

void
rbcast0_int(buff,type,root,size)
int *buff;
int type;
int root;
int size;
{
  PRINT('r',TYPENOINC(),root,size);
  bcast0(buff,size,TYPE(),root);
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
rbcast0_pointer(buff,type,root,size)
void *buff;
int type;
int root;
int size;
{
  PRINT('r',TYPENOINC(),root,size);
  bcast0(buff,size,TYPE(),root);
#if 0
#if defined(SUN) && defined(NIH)
  CTOHL(buff,size/sizeof(int *));
#endif
#endif
  PRINT_DATA('r',"0x%x\n",*(void**)buff);
  }

int *
rbcast0_test_pointer(type,root,size)
int type;
int root;
int size;
{
  int *buff;
  PRINT('r',TYPENOINC(),root,size);
  bcast0(&buff,size,TYPE(),root);
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
rbcast0_string(buff,type,root,size)
char **buff;
int type;
int root;
int size;
{
  int ilength = sizeof(int);
  int length;
  int i;

  for (i=0; i<size; i+=sizeof(*buff)) {
    /* Get the length of the string. */
    PRINT('r',TYPENOINC(),root,ilength);
    bcast0(&length,ilength,TYPE(),root);
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
      PRINT('r',TYPENOINC(),root,length);
      bcast0(*buff,length,TYPE(),root);
      PRINT_DATA('r',"%s\n",*buff);
      }
    buff++;
    }

  }

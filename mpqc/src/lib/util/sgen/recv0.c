/*
 * recv0.c
 *
 * Copyright (C) 1996 Limit Point Systems, Inc.
 *
 * Author: Curtis Janssen <cljanss@ca.sandia.gov>
 * Maintainer: LPS
 *
 * This file is part of the SC Toolkit.
 *
 * The SC Toolkit is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Library General Public License as published by
 * the Free Software Foundation; either version 2, or (at your option)
 * any later version.
 *
 * The SC Toolkit is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Library General Public License for more details.
 *
 * You should have received a copy of the GNU Library General Public License
 * along with the SC Toolkit; see the file COPYING.LIB.  If not, write to
 * the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 * The U.S. Government is granted a limited license as per AL 91-7.
 */

/* really just a copy of clj's rbcast0 routines */

/* $Log$
 * Revision 1.5  1996/10/25 19:38:39  etseidl
 * add copyleft notice and emacs local variables
 *
 * Revision 1.4  1995/03/18 00:11:28  cljanss
 * Using util/group to provide picl support.  Deleted the comm directory.
 *
 * Revision 1.3  1995/03/17  01:51:42  cljanss
 * Removed -I. and -I$(SRCDIR) from the default include path in
 * GlobalMakefile to avoid name conflicts with system include files.
 * Modified files under src.lib to include all files relative to src.lib.
 * Makefiles under src.bin need to add the -I. and -I$(SRCDIR) back onto
 * INCLUDE and CXXINCLUDE or make other arrangements.
 *
 * Revision 1.2  1994/08/25  22:48:26  etseidl
 * remove rcsids and fix some warnings
 *
 * Revision 1.1.1.1  1993/12/29  12:53:41  etseidl
 * SC source tree 0.1
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


#define NO_TEMPLATES
#include <stdio.h>
#include <stdlib.h>
#include <tmpl.h>
#include <util/group/picl.h>
#include <util/sgen/sgen.h>

#include <util/sgen/sndrcv0.h>

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
  PRINT_DATA('r',"%f\n",*buff);
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

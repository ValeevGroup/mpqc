/*
 * iseq.c
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

/* $Log$
 * Revision 1.4  1996/10/25 19:38:37  etseidl
 * add copyleft notice and emacs local variables
 *
 * Revision 1.3  1995/03/17 01:51:39  cljanss
 * Removed -I. and -I$(SRCDIR) from the default include path in
 * GlobalMakefile to avoid name conflicts with system include files.
 * Modified files under src.lib to include all files relative to src.lib.
 * Makefiles under src.bin need to add the -I. and -I$(SRCDIR) back onto
 * INCLUDE and CXXINCLUDE or make other arrangements.
 *
 * Revision 1.2  1994/08/25  22:48:21  etseidl
 * remove rcsids and fix some warnings
 *
 * Revision 1.1.1.1  1993/12/29  12:53:41  etseidl
 * SC source tree 0.1
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


#define NO_TEMPLATES
#include <stdio.h>
#include <math.h>
#include <util/sgen/sgen.h>


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

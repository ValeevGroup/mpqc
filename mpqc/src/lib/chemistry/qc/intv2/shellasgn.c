/*
 * shellasgn.c
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
 * Revision 1.8  1996/10/25 23:28:13  etseidl
 * add copyleft notice
 *
 * Revision 1.7  1995/11/16 00:47:39  cljanss
 * Removed normalization for individual basis functions.
 *
 * Revision 1.6  1995/10/25 21:19:55  cljanss
 * Adding support for pure am.  Gradients don't yet work.
 *
 * Revision 1.5  1995/03/17  01:49:36  cljanss
 * Removed -I. and -I$(SRCDIR) from the default include path in
 * GlobalMakefile to avoid name conflicts with system include files.
 * Modified files under src.lib to include all files relative to src.lib.
 * Makefiles under src.bin need to add the -I. and -I$(SRCDIR) back onto
 * INCLUDE and CXXINCLUDE or make other arrangements.
 *
 * Revision 1.4  1994/08/26  22:45:43  etseidl
 * fix a bunch of warnings, get rid of rcs id's, get rid of bread/bwrite and
 * fread/fwrite modules
 *
 * Revision 1.3  1994/05/27  18:13:51  cljanss
 * Fixed a bug that could lead to SEGV, etc, when assigning shells.
 *
 * Revision 1.2  1993/12/30  13:32:59  etseidl
 * mostly rcs id stuff
 *
 * Revision 1.3  1992/06/17  22:04:56  jannsen
 * cleaned up for saber-c
 *
 * Revision 1.2  1992/03/31  01:22:55  jannsen
 * Merged in Sandia non CVS codes.
 *
 * Revision 1.1  1992/01/17  22:12:11  cljanss
 * Initial revision
 *
 * Revision 1.1  1992/01/17  12:47:20  seidl
 * Initial revision
 * */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <tmpl.h>
#include <util/sgen/sgen.h>

#include <chemistry/qc/intv2/int_macros.h>

#include <chemistry/qc/intv2/atoms.h>

#include <chemistry/qc/intv2/atomsasgn.h>


int
assign_shell(_shell_1,_shell_2)
shell_t *_shell_1;
shell_t *_shell_2;
{
  typedef int boolean;
  typedef char * string;
  int i,j;
  int ncart;

  _shell_1->nprim=_shell_2->nprim;
  _shell_1->ncon=_shell_2->ncon;
  _shell_1->nfunc=_shell_2->nfunc;

  if(_shell_2->nprim!=0) {
    if(_shell_2->exp!=NULL) {
      _shell_1->exp = (double *) malloc(sizeof(double )*_shell_2->nprim);
      if(_shell_1->exp==NULL) return(-1);
      for(i=0;i<_shell_1->nprim;i++)
        _shell_1->exp[i]=_shell_2->exp[i];
      }
    else _shell_1->exp = NULL;
    }
  else _shell_1->exp = NULL; /* DT. */

  if(_shell_2->ncon!=0) {
    if(_shell_2->type!=NULL) {
      _shell_1->type = 
        (shell_type_t *) malloc(sizeof(shell_type_t )*_shell_2->ncon);
      if(_shell_1->type==NULL) return(-1);
      for (i=0; i<_shell_2->ncon; i++)  {
        assign_shell_type(&(_shell_1->type[i]),&(_shell_2->type[i]));
        }
      }
    else _shell_1->type = NULL;
    }
  else _shell_1->type = NULL; /* DT. */

  if(_shell_2->ncon!=0) {
    if(_shell_2->coef!=NULL) {
      _shell_1->coef = (double **) malloc(sizeof(double *)*_shell_2->ncon);
      if(_shell_1->coef==NULL) return(-1);
      for (i=0; i<_shell_2->ncon; i++)  {
        if(_shell_2->nprim!=0) {
          if(_shell_2->coef[i]!=NULL) {
            _shell_1->coef[i] = 
              (double *) malloc(sizeof(double )*_shell_2->nprim);
            if(_shell_1->coef[i]==NULL) return(-1);
            for(j=0;j<_shell_1->nprim;j++)
              _shell_1->coef[i][j]=_shell_2->coef[i][j];
            }
          }
        else _shell_1->coef[i] = NULL;
        }
      /* Skipped: else _shell_1->coef[i] = NULL;*/ /* DT. */
      }
    else _shell_1->coef = NULL;
    }
  else _shell_1->coef = NULL; /* DT. */

  return 0;
  }

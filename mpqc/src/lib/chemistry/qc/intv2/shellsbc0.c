/*
 * shellsbc0.c
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
 * Revision 1.7  1996/10/25 23:28:15  etseidl
 * add copyleft notice
 *
 * Revision 1.6  1995/11/16 00:47:42  cljanss
 * Removed normalization for individual basis functions.
 *
 * Revision 1.5  1995/10/25 21:19:59  cljanss
 * Adding support for pure am.  Gradients don't yet work.
 *
 * Revision 1.4  1995/03/17  01:49:39  cljanss
 * Removed -I. and -I$(SRCDIR) from the default include path in
 * GlobalMakefile to avoid name conflicts with system include files.
 * Modified files under src.lib to include all files relative to src.lib.
 * Makefiles under src.bin need to add the -I. and -I$(SRCDIR) back onto
 * INCLUDE and CXXINCLUDE or make other arrangements.
 *
 * Revision 1.3  1994/08/26  22:45:49  etseidl
 * fix a bunch of warnings, get rid of rcs id's, get rid of bread/bwrite and
 * fread/fwrite modules
 *
 * Revision 1.2  1993/12/30  13:33:06  etseidl
 * mostly rcs id stuff
 *
 * Revision 1.3  1992/06/17  22:05:15  jannsen
 * cleaned up for saber-c
 *
 * Revision 1.2  1992/03/31  01:24:06  jannsen
 * Merged in Sandia non CVS codes.
 *
 * Revision 1.1  1992/01/17  22:12:11  cljanss
 * Initial revision
 *
 * Revision 1.1  1992/01/17  12:47:20  seidl
 * Initial revision
 * */

#include <stdio.h>
#include <math.h>
#include <tmpl.h>
#include <util/sgen/sgen.h>

#include <chemistry/qc/intv2/int_macros.h>

#include <chemistry/qc/intv2/atoms.h>

#include <chemistry/qc/intv2/atomssbc0.h>

void
sbcast0_shell(_shell,_type,_root)
shell_t *_shell;
int _type;
int _root;
{
  int ncart;
  typedef int boolean;
  typedef char * string;
  int i;

  sbcast0_int(&(_shell->nprim),_type,_root, sizeof(int)*1);
  sbcast0_int(&(_shell->ncon),_type,_root, sizeof(int)*1);
  sbcast0_int(&(_shell->nfunc),_type,_root, sizeof(int)*1);

  if(_shell->nprim!=0) {
    sbcast0_pointer(&(_shell->exp),_type,_root,sizeof(double *));
    if(_shell->exp!=NULL) {
      sbcast0_double((_shell->exp),_type,_root,
        sizeof(double)*_shell->nprim);
      }
    }

  if(_shell->ncon!=0) {
    sbcast0_pointer(&(_shell->type),_type,_root,sizeof(shell_type_t *));
    if(_shell->type!=NULL) {
      for (i=0; i<_shell->ncon; i++)  {
        sbcast0_shell_type(&(_shell->type[i]),_type,_root);
        }
      }
    }

  if(_shell->ncon!=0) {
    sbcast0_pointer(&(_shell->coef),_type,_root,sizeof(double **));
    if(_shell->coef!=NULL) {
      for (i=0; i<_shell->ncon; i++)  {
        if(_shell->nprim!=0) {
          sbcast0_pointer(&(_shell->coef[i]),_type,_root,sizeof(double *));
          if(_shell->coef[i]!=NULL) {
            sbcast0_double((_shell->coef[i]),_type,_root,
              sizeof(double)*_shell->nprim);
            }
          }
        }
      }
    }

  return;
  }


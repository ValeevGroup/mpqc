/*
 * shellzero.c
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
 * Revision 1.7  1996/10/25 23:28:16  etseidl
 * add copyleft notice
 *
 * Revision 1.6  1995/11/16 00:47:43  cljanss
 * Removed normalization for individual basis functions.
 *
 * Revision 1.5  1995/03/17  01:49:41  cljanss
 * Removed -I. and -I$(SRCDIR) from the default include path in
 * GlobalMakefile to avoid name conflicts with system include files.
 * Modified files under src.lib to include all files relative to src.lib.
 * Makefiles under src.bin need to add the -I. and -I$(SRCDIR) back onto
 * INCLUDE and CXXINCLUDE or make other arrangements.
 *
 * Revision 1.4  1994/10/13  22:26:45  etseidl
 * replace bzero with memset
 *
 * Revision 1.3  1994/08/26  22:45:53  etseidl
 * fix a bunch of warnings, get rid of rcs id's, get rid of bread/bwrite and
 * fread/fwrite modules
 *
 * Revision 1.2  1993/12/30  13:33:08  etseidl
 * mostly rcs id stuff
 *
 * Revision 1.3  1992/06/17  22:05:19  jannsen
 * cleaned up for saber-c
 *
 * Revision 1.2  1992/03/31  01:24:36  jannsen
 * Merged in Sandia non CVS codes.
 *
 * Revision 1.1  1992/01/17  22:12:11  cljanss
 * Initial revision
 *
 * Revision 1.1  1992/01/17  12:47:20  seidl
 * Initial revision
 * */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <tmpl.h>
#include <util/sgen/sgen.h>

#include <chemistry/qc/intv2/int_macros.h>

#include <chemistry/qc/intv2/atoms.h>

#include <chemistry/qc/intv2/atomszero.h>

void
zero_shell(_shell)
shell_t *_shell;
{
  int nfunc;
  int i;

  if(_shell->nprim!=0) {
    if(_shell->exp!=NULL) {
      memset(_shell->exp,'\0',sizeof(double)*_shell->nprim);
      }
    }

  if(_shell->ncon!=0) {
    if(_shell->type!=NULL) {
      for (i=0; i<_shell->ncon; i++)  {
        zero_shell_type(&(_shell->type[i]));
        }
      }
    }

  if(_shell->ncon!=0) {
    if(_shell->coef!=NULL) {
      for (i=0; i<_shell->ncon; i++)  {
        if(_shell->nprim!=0) {
          if(_shell->coef[i]!=NULL) {
            memset(_shell->coef[i],'\0',sizeof(double)*_shell->nprim);
            }
          }
        }
      }
    }
  }


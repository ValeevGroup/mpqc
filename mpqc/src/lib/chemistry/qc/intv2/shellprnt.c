/*
 * shellprnt.c
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
 * Revision 1.7  1996/10/25 23:28:13  etseidl
 * add copyleft notice
 *
 * Revision 1.6  1995/11/16 00:47:40  cljanss
 * Removed normalization for individual basis functions.
 *
 * Revision 1.5  1995/10/25 21:19:56  cljanss
 * Adding support for pure am.  Gradients don't yet work.
 *
 * Revision 1.4  1995/03/17  01:49:36  cljanss
 * Removed -I. and -I$(SRCDIR) from the default include path in
 * GlobalMakefile to avoid name conflicts with system include files.
 * Modified files under src.lib to include all files relative to src.lib.
 * Makefiles under src.bin need to add the -I. and -I$(SRCDIR) back onto
 * INCLUDE and CXXINCLUDE or make other arrangements.
 *
 * Revision 1.3  1994/08/26  22:45:45  etseidl
 * fix a bunch of warnings, get rid of rcs id's, get rid of bread/bwrite and
 * fread/fwrite modules
 *
 * Revision 1.2  1993/12/30  13:33:04  etseidl
 * mostly rcs id stuff
 *
 * Revision 1.3  1992/06/17  22:05:06  jannsen
 * cleaned up for saber-c
 *
 * Revision 1.2  1992/03/31  01:23:39  jannsen
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

#include <chemistry/qc/intv2/atomsprnt.h>

extern int sgen_print_nindent;
#define SPI sgen_print_indent(fp)

void
print_shell(fp,_shell)
FILE *fp;
shell_t *_shell;
{
  int ncart;
  int i,j;
  int orig_indent = sgen_print_nindent;

  fprintf(fp,"nprim = ");
  sgen_print_nindent += 8;
  print_int(fp,&(_shell->nprim));
  fprintf(fp,"\n");
  sgen_print_nindent = orig_indent;
  SPI;

  fprintf(fp,"ncon = ");
  sgen_print_nindent += 7;
  print_int(fp,&(_shell->ncon));
  fprintf(fp,"\n");
  sgen_print_nindent = orig_indent;
  SPI;

  fprintf(fp,"nfunc = ");
  sgen_print_nindent += 8;
  print_int(fp,&(_shell->nfunc));
  fprintf(fp,"\n");
  sgen_print_nindent = orig_indent;
  SPI;

  fprintf(fp,"exp = ");
  sgen_print_nindent += 6;
  fprintf(fp,"[");
  sgen_print_nindent += 1;
  for (i=0; i<_shell->nprim; i++)  {
    print_double(fp,&(_shell->exp[i]));
    if(!((i+1)%8) && i) {
      fprintf(fp,"\n");
      SPI;
      }
    }
  fprintf(fp,"]\n");
  sgen_print_nindent += -1;
  sgen_print_nindent = orig_indent;
  SPI;

  fprintf(fp,"type:");
  sgen_print_nindent += 5;
  fprintf(fp,"[");
  sgen_print_nindent += 1;
  for (i=0; i<_shell->ncon; i++)  {
    if (i!=0) SPI;
    fprintf(fp,"(");
    sgen_print_nindent += 1;
    print_shell_type(fp,&(_shell->type[i]));
    SPI;
    fprintf(fp,")\n");
    sgen_print_nindent += -1;
    if(!((i+1)%8) && i) {
      fprintf(fp,"\n");
      SPI;
      }
    }
  SPI;
  fprintf(fp,"]\n");
  sgen_print_nindent += -1;
  sgen_print_nindent = orig_indent;
  SPI;

  fprintf(fp,"coef = ");
  sgen_print_nindent += 7;
  fprintf(fp,"[");
  sgen_print_nindent += 1;
  for (i=0; i<_shell->ncon; i++)  {
    if (i!=0) SPI;
    fprintf(fp,"[");
    sgen_print_nindent += 1;
    for (j=0; j<_shell->nprim; j++)  {
      print_double(fp,&(_shell->coef[i][j]));
      if(!((j+1)%8) && j) {
        fprintf(fp,"\n");
        SPI;
        }
      }
    fprintf(fp,"]\n");
    sgen_print_nindent += -1;
    }
  SPI;
  fprintf(fp,"]\n");
  sgen_print_nindent += -1;
  sgen_print_nindent = orig_indent;
  SPI;
  }

/*
 * basis.c
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

/* These utility routines assist in reading in the basis functions. */

/* $Log$
 * Revision 1.7  1996/10/25 23:27:50  etseidl
 * add copyleft notice
 *
 * Revision 1.6  1995/03/17 01:49:19  cljanss
 * Removed -I. and -I$(SRCDIR) from the default include path in
 * GlobalMakefile to avoid name conflicts with system include files.
 * Modified files under src.lib to include all files relative to src.lib.
 * Makefiles under src.bin need to add the -I. and -I$(SRCDIR) back onto
 * INCLUDE and CXXINCLUDE or make other arrangements.
 *
 * Revision 1.5  1994/08/26  22:45:09  etseidl
 * fix a bunch of warnings, get rid of rcs id's, get rid of bread/bwrite and
 * fread/fwrite modules
 *
 * Revision 1.4  1994/08/24  16:05:25  etseidl
 * get rid of ip functions, they have been replaced by keyval equivalents
 *
 * Revision 1.3  1994/05/27  23:44:12  cljanss
 * Changed some char* to const char*.  Included ipv2 interface from
 * keyval/ipv2c.h.  Fixed a bug with basis->n not being initialized.
 *
 * Revision 1.2  1993/12/30  13:32:43  etseidl
 * mostly rcs id stuff
 *
 * Revision 1.4  1993/04/28  00:30:13  jannsen
 * added int_read_basis global function
 *
 * Revision 1.3  1992/06/17  22:04:23  jannsen
 * cleaned up for saber-c
 *
 * Revision 1.2  1992/03/31  01:21:20  jannsen
 * Merged in Sandia non CVS codes.
 *
 * Revision 1.1  1991/12/14  00:19:36  cljanss
 * Initial revision
 * */
/*
 * Revision 1.6  91/10/31  14:32:24  cljanss
 * The basis set name is kept in the basis structure.
 * 
 * Revision 1.5  91/09/28  19:26:42  cljanss
 * Switch to new file naming scheme
 * 
 * Revision 1.4  91/09/10  19:32:34  cljanss
 * If the basis name is "nothing", an empty basis set
 * will be returned.
 * 
 * Revision 1.3  1991/08/09  16:41:23  cljanss
 * fixed basis set reader, it returned incorrect return codes
 *
 * Revision 1.2  1991/07/16  17:55:27  cljanss
 * slight change to make a character initialization compile on the SGI
 *
 * Revision 1.1  1991/06/16  16:40:07  janssen
 * Initial revision
 * */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <util/sgen/sgen.h>

#include <chemistry/qc/intv2/atoms.h>

/* This will print out a shell type.  The integer am is converted to
 * a character and if puream is nonnull, then the number of functions
 * in the shell is appended to the name.
 */
void
print_shell_type(FILE *fp, shell_type_t *shell_type)
{
  char *amnames = "spdfghijklmnoqrtuvwxyz";

  fprintf(fp," am = %c",amnames[shell_type->am]);
  if (shell_type->puream) fprintf(fp,"%d",2*shell_type->am + 1);

  sgen_print_suppress_indent();
}

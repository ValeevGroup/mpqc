/*
 * comp_erep23.c
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

/* These routines compute two and three center electron repulsion
 * integrals. */

/* $Log$
 * Revision 1.5  1996/10/25 23:28:03  etseidl
 * add copyleft notice
 *
 * Revision 1.4  1996/03/23 02:37:44  cljanss
 * Everything can now be configured with autoconf.
 *
 * Revision 1.3  1995/11/16 00:47:35  cljanss
 * Removed normalization for individual basis functions.
 *
 * Revision 1.2  1995/03/17  01:49:30  cljanss
 * Removed -I. and -I$(SRCDIR) from the default include path in
 * GlobalMakefile to avoid name conflicts with system include files.
 * Modified files under src.lib to include all files relative to src.lib.
 * Makefiles under src.bin need to add the -I. and -I$(SRCDIR) back onto
 * INCLUDE and CXXINCLUDE or make other arrangements.
 *
 * Revision 1.1  1994/05/27  23:51:25  cljanss
 * Added support for 2 and 3 center 2 electron integrals.  Added a test porgram.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <tmpl.h>
#include <math/array/math_lib.h>

#include <chemistry/qc/intv2/atoms.h>
#include <chemistry/qc/intv2/int_flags.h>
#include <chemistry/qc/intv2/int_macros.h>
#include <chemistry/qc/intv2/int_types.h>

#include <chemistry/qc/intv2/inter.h>

#include <chemistry/qc/intv2/comp_erep.gbl>

#include <chemistry/qc/intv2/comp_erep23.gbl>
#include <chemistry/qc/intv2/comp_erep23.lcl>

LOCAL_FUNCTION void
make_int_unit_shell()
{
  int_unit_shell = (shell_t*) malloc(sizeof(shell_t));
  int_unit_shell->nprim = 1;
  int_unit_shell->ncon = 1;
  int_unit_shell->nfunc = 1;
  int_unit_shell->exp = (double*) malloc(sizeof(double));
  int_unit_shell->exp[0] = 0.0;
  int_unit_shell->type = (shell_type_t*) malloc(sizeof(shell_type_t));
  int_unit_shell->type[0].am = 0;
  int_unit_shell->type[0].puream = 0;
  int_unit_shell->coef = (double**) malloc(sizeof(double*));
  int_unit_shell->coef[0] = (double*) malloc(sizeof(double));
  int_unit_shell->coef[0][0] = 1.0;
}

/* Compute a 2 center electron repulsion integral.  Electron 1 is in
 * shell psh1 and electron 2 is in psh2, that is (1 | 2).  To avoid
 * confusing the user of these routines, the INT_NOPERM is set.
 */
GLOBAL_FUNCTION void
int_erep2(flags,psh1,psh2)
    int flags;
    int *psh1;
    int *psh2;
{
  centers_t *cs2 = int_cs2;
  centers_t *cs4 = int_cs4;
  int shd = 0x11111111; /* a dummy shell that will cause death if used */
  if (!int_unit_shell) make_int_unit_shell();
  int_cs2 = (centers_t*) 0x0;
  int_cs4 = (centers_t*) 0x0;
  int_unit2 = 1;
  int_unit4 = 1;
  int_erep(flags,psh1,&shd,psh2,&shd);
  int_unit2 = 0;
  int_unit4 = 0;
  int_cs2 = cs2;
  int_cs4 = cs4;
}

/* This is an alternate interface to int_erep2.  It takes
 * as arguments the flags, an integer vector of shell numbers
 * and an integer vector which will be filled in with size
 * information, if it is non-NULL. */
GLOBAL_FUNCTION void
int_erep2_v(flags,shells,sizes)
int flags;
int *shells;
int  *sizes;
{
  int_erep2(flags|INT_NOPERM,&(shells[0]),&(shells[1]));
  if (sizes) {
      sizes[0] = INT_SH(int_cs1,shells[0]).nfunc;
      sizes[1] = INT_SH(int_cs3,shells[1]).nfunc;
    }
}


/* Computes a 3 center two electron integral.  Electron 1 is in psh1
 * and electron 2 is in psh2 and psh3, that is (1 | 2 3).  To avoid
 * confusing the user of these routines, the INT_NOPERM is set.
 */
GLOBAL_FUNCTION void
int_erep3(flags,psh1,psh2,psh3)
    int flags;
    int *psh1;
    int *psh2;
    int *psh3;
{
  int shd = 0x11111111; /* a dummy shell that will cause death if used */
  if (!int_unit_shell) make_int_unit_shell();
  int_unit2 = 1;
  int_erep(flags|INT_NOPERM,psh1,&shd,psh2,psh3);
  int_unit2 = 0;
}

/* This is an alternate interface to int_erep3.  It takes
 * as arguments the flags, an integer vector of shell numbers
 * and an integer vector which will be filled in with size
 * information, if it is non-NULL. */
GLOBAL_FUNCTION void
int_erep3_v(flags,shells,sizes)
int flags;
int *shells;
int  *sizes;
{
  int_erep3(flags|INT_NOPERM,&(shells[0]),&(shells[1]),&(shells[2]));
  if (sizes) {
      sizes[0] = INT_SH(int_cs1,shells[0]).nfunc;
      sizes[1] = INT_SH(int_cs3,shells[1]).nfunc;
      sizes[2] = INT_SH(int_cs4,shells[2]).nfunc;
    }
}

/*
 * comp_0e.c
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
 * Revision 1.7  1996/10/25 23:27:57  etseidl
 * add copyleft notice
 *
 * Revision 1.6  1996/03/23 02:37:37  cljanss
 * Everything can now be configured with autoconf.
 *
 * Revision 1.5  1995/03/17 01:49:24  cljanss
 * Removed -I. and -I$(SRCDIR) from the default include path in
 * GlobalMakefile to avoid name conflicts with system include files.
 * Modified files under src.lib to include all files relative to src.lib.
 * Makefiles under src.bin need to add the -I. and -I$(SRCDIR) back onto
 * INCLUDE and CXXINCLUDE or make other arrangements.
 *
 * Revision 1.4  1995/03/05  06:05:29  cljanss
 * Added efield integrals.  Changed the dipole moment integral interface.
 *
 * Revision 1.3  1994/08/26  22:45:17  etseidl
 * fix a bunch of warnings, get rid of rcs id's, get rid of bread/bwrite and
 * fread/fwrite modules
 *
 * Revision 1.2  1993/12/30  13:32:46  etseidl
 * mostly rcs id stuff
 *
 * Revision 1.1.1.1  1992/03/17  16:32:35  seidl
 * DOE-NIH Quantum Chemistry Library 0.0
 *
 * Revision 1.1  1992/03/17  16:32:34  seidl
 * Initial revision
 *
 * Revision 1.1  1992/02/03  15:49:49  seidl
 * Initial revision
 *
 * Revision 1.1  1991/12/14  00:19:36  cljanss
 * Initial revision
 * */
/*
 * Revision 1.3  91/09/28  19:26:45  cljanss
 * Switch to new file naming scheme
 * 
 * Revision 1.2  91/06/21  18:38:33  janssen
 * put in first derivatives of nuclear repulsion energy
 * 
 * Revision 1.1  1991/06/16  16:40:07  janssen
 * Initial revision
 * */

#include <stdio.h>
#include <math.h>
#include <tmpl.h>
#include <math/array/math_lib.h>
#include <chemistry/qc/intv2/atoms.h>
#include <chemistry/qc/intv2/int_macros.h>

#include <chemistry/qc/intv2/comp_0e.gbl>
#include <chemistry/qc/intv2/comp_0e.lcl>

/* Compute the nuclear repulsion energy between two sets of centers. */
GLOBAL_FUNCTION double
int_nuclear_repulsion(cs1,cs2)
centers_t *cs1;
centers_t *cs2;
{
  int i,j,k;
  double r,r2;
  double result;

  result = 0.0;
  for (i=0; i<cs1->n; i++) {
    for (j=0; j<(cs1==cs2?i:cs2->n); j++) {
      r2 = 0.0;
      for (k=0; k<3; k++) {
        r = cs1->center[i].r[k] - cs2->center[j].r[k];
        r2 += r*r;
        }
      result += cs1->center[i].charge * cs2->center[j].charge / sqrt(r2);
      }
    }

  return result;
  }

/* This computes the efield at position due to the nuclei.  The result
 * is written to efield.
 */
GLOBAL_FUNCTION void
int_nuclear_efield(cs1,cs2,position,efield)
centers_t *cs1;
centers_t *cs2;
double *position;
double *efield;
{
  int i,j;
  double tmp;
  double r[3];
  centers_t *centers;

  for (i=0; i<3; i++) efield[i] = 0.0;

  centers = cs1;
  while(centers) {
      for (i=0; i<centers->n; i++) {
          tmp = 0.0;
          for (j=0; j<3; j++) {
              r[j] = position[j] - centers->center[i].r[j];
              tmp += r[j]*r[j];
            }
          tmp = centers->center[i].charge/(tmp*sqrt(tmp));
          for (j=0; j<3; j++) {
              efield[j] +=  r[j] * tmp;
            }
        }
      if (centers == cs2) centers = NULL;
      else centers = cs2;
    }
}

/* Compute the nuclear repulsion energy first derivative with respect
 * to the given center. */
GLOBAL_FUNCTION void
int_nuclear_repulsion_1der(cs1,cs2,repder,csder,cdernum)
centers_t *cs1;
centers_t *cs2;
double_vector_t *repder;
centers_t *csder;
int cdernum;
{
  int i,j,k;
  double r[3],r2;
  double factor;

  repder->d[0] = 0.0;
  repder->d[1] = 0.0;
  repder->d[2] = 0.0;
  for (i=0; i<cs1->n; i++) {
    for (j=0; j<(cs1==cs2?i:cs2->n); j++) {
      if (((cs1==csder)&&(cdernum==i))||((cs2==csder)&&(cdernum==j))) {
        r2 = 0.0;
        for (k=0; k<3; k++) {
          r[k] = cs1->center[i].r[k] - cs2->center[j].r[k];
          r2 += r[k]*r[k];
          }
        factor = - cs1->center[i].charge * cs2->center[j].charge * pow(r2,-1.5);
        if ((cs2==csder)&&(cdernum==j)) factor = -factor;
        for (k=0; k<3; k++) {
          repder->d[k] +=  factor * r[k];
          }
        }
      }
    }
  }

/*
 * skeleton.c
 *
 * Copyright (C) 1996 Limit Point Systems, Inc.
 *
 * Author: Edward Seidl <seidl@janed.com>
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

#include <stdio.h>
#include <math.h>
#include <tmpl.h>
#include <math/array/math_lib.h>
#include <chemistry/qc/intv2/int_libv2.h>
#include <math/dmt/libdmt.h>

#include <chemistry/qc/dmtsym/symm.h>

#include <chemistry/qc/dmtsym/skeleton.gbl>
#include <chemistry/qc/dmtsym/skeleton.lcl>

/************************************************************************
 * 
 * this is really just a debugging function.  given a symmetric dmt matrix
 * and the sym struct, create a skeleton matrix a la Dupuis and King.
 * mat is overwritten with the skeleton matrix
 *
 */

GLOBAL_FUNCTION void
sym_skel_matrix(sym_info,mat)
sym_struct_t *sym_info;
dmt_matrix mat;
{
  int i,j;
  int si,sj,sij;
  int isiz,jsiz;
  int nlocal;
  double lij;
  double *tmp;

  if (dmt_distribution(mat) != SCATTERED) {
    fprintf(stderr,"sym_skel_matrix: want scattered matrices only\n");
    exit(-1);
  }

  nlocal = dmt_nlocal(mat);

  for (i=0; i < nlocal ; i++) {
    dmt_get_block_dsc(mat,i,&si,&isiz,&sj,&jsiz,&tmp);
    sij = (si > sj) ? si*(si+1)/2 + sj : sj*(sj+1)/2 +si;
    lij = (double) sym_info->lamij[sij];
    for (j=0; j < isiz*jsiz ; j++) tmp[j] *= lij;
  }
}

/* Local Variables:
 * mode: c
 * eval: (c-set-style "ETS")
 * End:
 */

/*
 * expewdens.c
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

#include <stdio.h>
#include <math.h>
#include <tmpl.h>
#include <math/dmt/libdmt.h>

#include <chemistry/qc/dmtqc/expewdens.gbl>
#include <chemistry/qc/dmtqc/expewdens.lcl>

/* Computes occupied and virtual density matrices with a weighting
 * factor exp(-epsilon * t) for virtuals and exp(epsilon*t) for
 * occupied where epsilon is an orbital eigenvalue.
 */
GLOBAL_FUNCTION void
dmt_expewdensity (nocc, C, eigval, D, E, t)
int nocc;
dmt_matrix C;
double *eigval;
dmt_matrix D;
dmt_matrix E;
double t;
{
  double *col;
  int i;
  loop_t *loop;
  int iind,isize,jsize;

  dmt_fill (D, 0.0);
  dmt_fill (E, 0.0);
  loop = dmt_ngl_create("%mr",C);
  while(dmt_ngl_next(loop)) {
    dmt_ngl_create_inner(loop,0);
    while(dmt_ngl_next_inner_m(loop,&iind,&isize,&i,&jsize,&col)) {
      if (i < nocc) {
        /*Sum contribution from |col| into local |D| blocks.*/
        int localp = dmt_nlocal (D);
        int il;
        int mub, nub, mu, nu, mustart, nustart, musize, nusize;
        double *block;
        for (il = 0; il<localp; il++) {
          dmt_get_block (D, il, &mub, &nub, &block);
          dmt_describe_block (D, mub, &mustart, &musize);
          dmt_describe_block (D, nub, &nustart, &nusize);
          for (mu = 0; mu<musize; mu++) {
            for (nu = 0; nu<nusize; nu++) {
              block[mu*nusize+nu] +=   col[nustart+nu]*col[mustart+mu]
                                     * exp(eigval[i]*t);
              }
            }
          }
        }
      else {
        /*Sum contribution from |col| into local |E| blocks.*/
        int localp = dmt_nlocal (E);
        int il;
        int mub, nub, mu, nu, mustart, nustart, musize, nusize;
        double *block;
        for (il = 0; il<localp; il++) {
          dmt_get_block (E, il, &mub, &nub, &block);
          dmt_describe_block (E, mub, &mustart, &musize);
          dmt_describe_block (E, nub, &nustart, &nusize);
          for (mu = 0; mu<musize; mu++) {
            for (nu = 0; nu<nusize; nu++) {
              block[mu*nusize+nu] +=   col[nustart+nu]*col[mustart+mu]
                                     * exp(-eigval[i]*t);
              }
            }
          }
        }
      }
    }
  dmt_ngl_kill(loop);

  }


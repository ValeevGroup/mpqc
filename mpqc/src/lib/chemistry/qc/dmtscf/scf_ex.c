/*
 * scf_ex.c
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

/* This computes the exchange energy of the wavefunction. */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include <tmpl.h>
#include <util/group/picl.h>
#include <math/array/math_lib.h>
#include <util/misc/libmisc.h>
#include <chemistry/qc/intv2/int_libv2.h>
#include <chemistry/qc/dmtsym/sym_dmt.h>
#include <chemistry/qc/dmtqc/libdmtqc.h>

#include <chemistry/qc/dmtscf/scf.h>

#include <chemistry/qc/dmtscf/scf_ex.gbl>
#include <chemistry/qc/dmtscf/scf_ex.lcl>

/*********************************************************************
 *
 * calculate the exchange energy given an scf vector
 *
 * input:
 *   scf_info = pointer to initialized scf struct
 *   centers  = pointer to initialized centers struct
 *   Scf_Vec  = column distributed matrix containing converged wavefunction
 *
 * on return:
 *   scf_info.e_exc contains the exchange energy
 *
 */

GLOBAL_FUNCTION void
scf_ex(scf_info,centers,Scf_Vec)
scf_struct_t *scf_info;
centers_t*centers;
dmt_matrix Scf_Vec;
{
  loop_t *loop;
  int ish,jsh,ksh,lsh, isz,jsz,ksz,lsz;
  int i,j,k,l;
  int nl;
  double *pkl, *pij;
  double *buf;
  double exchange1 = 0.0;
  double exchange2 = 0.0;
  int ij,kl;
  int index;
  int nlocal;
  int flags;

  dmt_matrix Pmat;

  assert(dmt_distribution(Scf_Vec) == COLUMNS);

 /* let's form the density from the scf vector */
  Pmat = dmt_create("scf_ex density matrix",scf_info->nbfao,SCATTERED);

  dmt_density(Scf_Vec,scf_info->nclosed,Pmat);
  dmt_scale(Pmat,4.0);
  dmt_scale_diagonal(Pmat,0.5);

  nlocal=dmt_nlocal(Pmat);

  if (scf_info->iopen) {
    fprintf(stderr,"scf_ex: can't compute exchange energy for open shells\n");
    return;
  }

  tim_enter("scf_ex");

  /* initialize the integrals routines */
  int_initialize_offsets2(centers,centers,centers,centers);

  flags = INT_EREP|INT_NOSTRB;
  if (!scf_info->int_store1) flags |= INT_NOSTR1;
  if (!scf_info->int_store2) flags |= INT_NOSTR2;
  buf = int_initialize_erep(flags,0,centers,centers,centers,centers);

  loop = dmt_ngl_create("%m",Pmat);
  while (dmt_ngl_next(loop)) {
    dmt_ngl_create_inner(loop,0);
    while (dmt_ngl_next_inner_m(loop,&ksh,&ksz,&lsh,&lsz,&pkl)) {
      double factor,factor2;

      for (nl=0; nl < nlocal; nl++) {
        dmt_get_block_dsc(Pmat,nl,&ish,&isz,&jsh,&jsz,&pij);

        if (ish != jsh && ksh != lsh)
          factor = 2.0;
        else
          factor = 1.0;

        if (ish != jsh) factor *= 0.5;
        if (ksh != lsh) factor *= 0.5;

       /* compute the integrals for this shell block */
        int_erep(INT_EREP|INT_NOBCHK|INT_NOPERM|INT_REDUND,&ish,&ksh,&jsh,&lsh);

        index = 0;
        for (i=0; i<isz; i++) {
          for (k=0; k<ksz; k++) {
            for (j=0; j<jsz; j++) {
              ij = i*jsz + j;
              for (l=0; l<lsz; l++) {
                kl = k*lsz + l;

                factor2 = factor;
                if (ish == jsh && i != j) factor2 *= 0.5;
                if (ksh == lsh && k != l) factor2 *= 0.5;

                exchange1 += buf[index]*pij[ij]*pkl[kl]*factor2;
                index++;
              }
            }
          }
        }

        if (ish == jsh && ksh == lsh) continue;

        int_erep(INT_EREP|INT_NOBCHK|INT_NOPERM|INT_REDUND,&ish,&lsh,&ksh,&jsh);

        index = 0;
        for (i=0; i<isz; i++) {
          for (l=0; l<lsz; l++) {
            for (k=0; k<ksz; k++) {
              kl = k*lsz + l;
              for (j=0; j<jsz; j++) {
                ij = i*jsz + j;

                factor2 = factor;
                if (ish == jsh && i != j) factor2 *= 0.5;
                if (ksh == lsh && k != l) factor2 *= 0.5;

                exchange2 += buf[index]*pij[ij]*pkl[kl]*factor2;
                index++;
              }
            }
          }
        }
      }
    }
  }

 /* we are finished with loop */
  dmt_ngl_kill(loop);

  int_done_erep();
  int_done_offsets2(centers,centers,centers,centers);

  gsum0(&exchange1,1,5,mtype_get(),0);
  gsum0(&exchange2,1,5,mtype_get(),0);

  tim_exit("scf_ex");

  scf_info->e_exc = -0.25*(exchange1+exchange2);
}

/*************************************************************************
 * Local Variables:
 * mode: c
 * eval: (c-set-style "ETS")
 * End:
 */

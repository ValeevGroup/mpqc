
/* $Log$
 * Revision 1.2  1994/08/26 22:47:59  etseidl
 * get rid of rcsids
 *
 * Revision 1.1.1.1  1993/12/29  12:53:00  etseidl
 * SC source tree 0.1
 *
 * Revision 1.1  1992/06/17  21:29:23  jannsen
 * moved from libdmt/matrix.web
 *
 * */

#include <stdio.h>
#include <tmpl.h>
#include <math/dmt/libdmt.h>

#include "dens.gbl"
#include "dens.lcl"

/* This computes the density */

GLOBAL_FUNCTION VOID
dmt_density (C, ndoc, P)
dmt_matrix C;
int ndoc;
dmt_matrix P;
{
  double *col;
  int i;
  loop_t *loop;
  int iind,isize,jsize;
  int localp = dmt_nlocal (P);

  dmt_fill (P, 0.0);
  loop = dmt_ngl_create("%mr",C);
  while(dmt_ngl_next(loop)) {
    dmt_ngl_create_inner(loop,0);
    while(dmt_ngl_next_inner_m(loop,&iind,&isize,&i,&jsize,&col)) {
      if (i < ndoc) {
        /*Sum contribution from |col| into local |P| blocks.*/
        int il;
        int mub, nub, mu, nu, mustart, nustart, musize, nusize;
        double *block;
        for (il = 0; il<localp; il++) {
          dmt_get_block (P, il, &mub, &nub, &block);
          dmt_describe_block (P, mub, &mustart, &musize);
          dmt_describe_block (P, nub, &nustart, &nusize);
          for (mu = 0; mu<musize; mu++) {
            for (nu = 0; nu<nusize; nu++) {
              block[mu*nusize+nu] += col[nustart+nu]*col[mustart+mu];
              }
            }
          }
        }
      }
    }
  dmt_ngl_kill(loop);

  }


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


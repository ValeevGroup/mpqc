
#include <stdio.h>
#include <stdlib.h>
#include <tmpl.h>
#include <math/dmt/libdmt.h>

#include <chemistry/qc/dmtqc/ewdensity.gbl>
#include <chemistry/qc/dmtqc/ewdensity.lcl>

/* This computes the energy weight density for clscf
 * wavefunctions. */
GLOBAL_FUNCTION VOID
dmt_ewdensity (C, Fmo, nocc, W)
dmt_matrix C;
dmt_matrix Fmo;
int nocc;
dmt_matrix W;
{
  double *col;
  int i;
  int nbasis = dmt_size(C);
  double *diag;
  loop_t *loop;
  int iind,isize,jsize;

  diag = (double *) malloc(sizeof(double)*nbasis);
  dmt_get_diagonal(Fmo,diag);

  dmt_fill (W, 0.0);
  loop = dmt_ngl_create("%mr",C);
  while(dmt_ngl_next(loop)) {
    dmt_ngl_create_inner(loop,0);
    while(dmt_ngl_next_inner_m(loop,&iind,&isize,&i,&jsize,&col)) {
      if (i < nocc) {
        /*Sum contribution from |col| into local |W| blocks.*/
        int localp = dmt_nlocal (W);
        int il;
        int mub, nub, mu, nu, mustart, nustart, musize, nusize;
        double *block;
        for (il = 0; il<localp; il++) {
          dmt_get_block (W, il, &mub, &nub, &block);
          dmt_describe_block (W, mub, &mustart, &musize);
          dmt_describe_block (W, nub, &nustart, &nusize);
          for (mu = 0; mu<musize; mu++) {
            for (nu = 0; nu<nusize; nu++) {
              block[mu*nusize+nu] += col[nustart+nu]*col[mustart+mu] * diag[i];
              }
            }
          }
        }
      }
    }
  dmt_ngl_kill(loop);

  free(diag);
  }


/* This computes the energy weighted density for grscf
 * wavefunctions.  The open shell ewd matrix is really just the
 * AO Lagrangian, to within a sign.
 */
GLOBAL_FUNCTION VOID
dmt_open_ewdensity (C, F, Fo, ndocc, nsocc, W)
dmt_matrix C;
dmt_matrix F;
dmt_matrix Fo;
int ndocc;
int nsocc;
dmt_matrix W;
{
  double *Lblk,*Fblk,*FOblk;
  int i,j;
  int nlb,ib,jb,isz,jsz,ist,jst;
  int occi,occj;
  int nbasis = dmt_size(F);
  int nloc=dmt_nlocal (F);
  int nocc=ndocc+nsocc;
  dmt_matrix Ct,WC,Lag,SCR;

/* First form the MO lagrangian.  This has the form
        c    o   v
   c  |2*FC|2*FC|0|
      -------------
   o  |2*FC| FO |0|
      -------------
   v  | 0  |  0 |0|
 */

  Lag = dmt_create("MO lagrangian",nbasis,SCATTERED);

  dmt_fill(Lag,0.0);

  for(nlb=0; nlb < nloc ; nlb++) {
    dmt_get_block(Lag,nlb,&ib,&jb,&Lblk);
    dmt_get_block(F,nlb,&ib,&jb,&Fblk);
    dmt_get_block(Fo,nlb,&ib,&jb,&FOblk);
    dmt_describe_block(F,ib,&ist,&isz);
    dmt_describe_block(F,jb,&jst,&jsz);
  
    for(i=0; i < isz ; i++) {
      occi=0;
      if(i+ist < ndocc) occi=2;
      else if(i+ist < nocc) occi=1;
      for(j=0; j < jsz ; j++) {
        occj=0;
        if(j+jst < ndocc) occj=2;
        else if(j+jst < nocc) occj=1;
        if(occi==2 && occj==2)
          Lblk[i*jsz+j] = 2.0*Fblk[i*jsz+j];
        else if(occi==1 && occj==1)
          Lblk[i*jsz+j] = FOblk[i*jsz+j];
        else if(occi && occj)
          Lblk[i*jsz+j] = 2.0*Fblk[i*jsz+j];
        else
          Lblk[i*jsz+j] = 0.0;
        }
      }
    }

/* now form W=C*Lag(MO)*C~ */
  Ct = dmt_create("transposed vector",nbasis,COLUMNS);
  dmt_copy(C,Ct);
  dmt_transpose(Ct);

  WC = dmt_create("ewdens col",nbasis,COLUMNS);
  SCR = dmt_create("scr col",nbasis,COLUMNS);

  dmt_mult(Lag,Ct,SCR);
  dmt_mult(Ct,SCR,WC);

  dmt_copy(WC,W);
  dmt_scale(W,-1.0);

  dmt_free(Ct);
  dmt_free(WC);
  dmt_free(SCR);
  dmt_free(Lag);
  }


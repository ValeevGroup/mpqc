
/* $Log$
 * Revision 1.1  1993/12/29 12:53:15  etseidl
 * Initial revision
 *
 * Revision 1.5  1992/06/23  20:04:23  seidl
 * change dmt matrices to uppercase,
 * get rid of unnecessary matrice
 *
 * Revision 1.4  1992/06/17  21:54:09  jannsen
 * cleaned up for saber-c
 *
 * Revision 1.3  1992/05/26  20:17:41  jannsen
 * use mtype_get to get message types for global operations
 * check results of memory allocations
 *
 * Revision 1.2  1992/05/19  20:56:24  seidl
 * use message types 7000-7999
 *
 * Revision 1.1.1.1  1992/03/17  16:26:04  seidl
 * DOE-NIH Quantum Chemistry Library 0.0
 *
 * Revision 1.1  1992/03/17  16:26:00  seidl
 * Initial revision
 *
 * Revision 1.2  1992/02/10  16:58:11  seidl
 * broadcast neelec and delta so all nodes will know when converged
 *
 * Revision 1.1  1992/02/04  23:48:08  seidl
 * Initial revision
 *
 * Revision 1.4  1992/01/16  19:53:07  seidl
 * use new integral routines
 *
 * Revision 1.3  1992/01/13  19:11:28  seidl
 * move test for open-shell out of loop, add comments
 *
 * Revision 1.2  1992/01/09  11:38:02  seidl
 * add parallel code
 *
 * Revision 1.1  1991/12/20  16:23:05  seidl
 * Initial revision
 *
 * Revision 1.2  1991/06/22  08:58:55  seidl
 * op is no longer a pointer
 *
 * Revision 1.1  1991/06/19  13:01:36  seidl
 * Initial revision
 * */

static char rcsid[]="$Id$";

#include <stdio.h>
#include <tmpl.h>
#include <math.h>
#include <comm/picl/picl.h>
#include <math/dmt/libdmt.h>
#include "scf.h"

#include "scf_en.gbl"
#include "scf_en.lcl"

static double plimit;

GLOBAL_FUNCTION int
scf_calculate_energy(_scf_info,FOCK,PMAT,DPMAT,GMATO,PMATO,SSCR1,iter,_outfile)
scf_struct_t *_scf_info;
dmt_matrix FOCK;
dmt_matrix PMAT;
dmt_matrix DPMAT;
dmt_matrix GMATO;
dmt_matrix PMATO;
dmt_matrix SSCR1;
int iter;
FILE *_outfile;
{
  int i,j,ij;
  int nlocal,nl;
  int ib,jb,isz,jsz;
  double edif,etot;
  double neelec = 0.0;
  double delta;
  double *fblk,*pblk,*dpblk,*gblk,*poblk;
  dmt_matrix HCORE=dmt_old("libscfv3 hcore matrix");

  if(!iter) plimit = pow(10.0,(double) -(_scf_info->convergence));

/* calculate the total electronic energy
 * for closed shell, 
 *    Eelec = sum(i,j) D(i,j)*{H(i,j)+F(i,j)}
 *
 * for open shell, 
 *    Eelec = sum(i,j) D(i,j)*{H(i,j)+F(i,j)} - 1/2 Do(i,j)*Go(i,j)
 */

  dmt_copy(HCORE,SSCR1);
  dmt_sum(FOCK,SSCR1);

  nlocal = dmt_nlocal(SSCR1);

  delta=0.0;
  neelec=0.0;
  for(nl=0; nl < nlocal ; nl++) {
    dmt_get_block_dsc(SSCR1,nl,&ib,&isz,&jb,&jsz,&fblk);
    dmt_get_block_dsc(PMAT,nl,&ib,&isz,&jb,&jsz,&pblk);

    if(!_scf_info->iopen) {
      if(ib!=jb) {
        for(i=0; i < isz*jsz ; i++) {
          neelec += pblk[i]*fblk[i];
          }
        }
      else {
        for(i=0; i < isz ; i++) {
          for(j=0; j <= i ; j++) {
            ij=i*jsz+j;
            neelec += pblk[ij]*fblk[ij];
            }
          }
        }
      }
    else {
      dmt_get_block_dsc(PMATO,nl,&ib,&isz,&jb,&jsz,&poblk);
      dmt_get_block_dsc(GMATO,nl,&ib,&isz,&jb,&jsz,&gblk);
      if(ib!=jb) {
        for(i=0; i < isz*jsz ; i++) neelec += 
                                pblk[i]*fblk[i]-poblk[i]*gblk[i];
        }
      else {
        for(i=0; i < isz ; i++) {
          for(j=0; j <= i ; j++) {
            ij=i*jsz+j;
            neelec += pblk[ij]*fblk[ij]-poblk[ij]*gblk[ij];
            }
          }
        }
      }

    dmt_get_block_dsc(DPMAT,nl,&ib,&isz,&jb,&jsz,&dpblk);

    if(ib!=jb) for(i=0; i < isz*jsz ; i++) delta += dpblk[i]*dpblk[i];
    else {
      for(i=0; i < isz ; i++) {
        for(j=0; j <= i ; j++) {
          ij=i*jsz+j;
          delta += dpblk[ij]*dpblk[ij];
          }
        }
      }
    }

  neelec*=0.5;
  gsum0(&neelec,1,5,mtype_get(),0);
  bcast0(&neelec,sizeof(double),mtype_get(),0);

  gsum0(&delta,1,5,mtype_get(),0);
  bcast0(&delta,sizeof(double),mtype_get(),0);

  delta = sqrt(delta)/_scf_info->mxcoef2;
  etot = _scf_info->nuc_rep + neelec;
  edif = _scf_info->e_elec - neelec;

  if(mynode0()==0) {
    if (!iter) {
      fprintf(_outfile,"\n  iter       total energy       ");
      fprintf(_outfile," delta E         delta P          diiser\n");
      }

    fprintf(_outfile, "%5d %20.10f %15.6e %15.6e %15.6e\n", 
                   iter+1, etot, edif, delta, _scf_info->diis_er);
    fflush(_outfile);
    }

  if ( delta < plimit ) _scf_info->converged=1;

  _scf_info->e_elec = neelec;

#if 0
 /* this is for TCSCF, it determines when new ci coefficients should
  * be calculated */

  cinext = pow(10.0,-twocut);
  if (delta < cinext && delta && !converged) {
     twocut += incr;
     return(1);
     }
  else return(0);
#endif
  return(0);
  }

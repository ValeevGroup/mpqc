
#include <stdio.h>
#include <tmpl.h>
#include <math.h>
#include <util/misc/libmisc.h>
#include <util/group/picl.h>
#include <math/dmt/libdmt.h>
#include <chemistry/qc/dmtscf/scf.h>

#include <chemistry/qc/dmtscf/scf_en.gbl>
#include <chemistry/qc/dmtscf/scf_en.lcl>

/*****************************************************************************
 *
 * calculate the total electronic energy
 * for closed shell, 
 *    Eelec = sum(i,j) D(i,j)*{H(i,j)+F(i,j)}
 *
 * for open shell, 
 *    Eelec = sum(i,j) D(i,j)*{H(i,j)+F(i,j)} - 1/2 Do(i,j)*Go(i,j)
 *
 * input:
 *   scf_info = pointer to scf struct
 *   Hcore    = scattered dmt matrix containing the core hamiltonian
 *   Fock     = scattered dmt matrix containing the AO Fock matrix
 *   Pmat     = scattered dmt matrix containing the AO density matrix
 *   FockO    = scattered dmt matrix containing the open-shell Fock matrix
 *   PmatO    = scattered dmt matrix containing the open-shell density matrix
 *   SScr1    = scratch scattered dmt matrix
 *   SScr2    = scratch scattered dmt matrix
 *
 * return Eelec
 */

GLOBAL_FUNCTION double
scf_electronic_energy(scf_info,Hcore,Fock,Pmat,FockO,PmatO,SScr1,SScr2)
scf_struct_t *scf_info;
dmt_matrix Hcore;
dmt_matrix Fock;
dmt_matrix Pmat;
dmt_matrix FockO;
dmt_matrix PmatO;
dmt_matrix SScr1;
dmt_matrix SScr2;
{
  int i,j,ij;
  int nlocal,nl;
  int ib,jb,isz,jsz;
  double neelec;
  double *fblk,*pblk,*gblk,*poblk;

  assert(dmt_distribution(Hcore) == SCATTERED);
  assert(dmt_distribution(Fock) == SCATTERED);
  assert(dmt_distribution(Pmat) == SCATTERED);
  assert(dmt_distribution(SScr1) == SCATTERED);
  assert(dmt_distribution(SScr2) == SCATTERED);
  if (scf_info->iopen) {
    assert(dmt_distribution(FockO) == SCATTERED);
    assert(dmt_distribution(PmatO) == SCATTERED);
  }

 /* ok, FockO = H + G - GO, so to get GO back, add -(Fock = H + G) to FockO */
  if (scf_info->iopen) {
    dmt_copy(Fock,SScr1);
    dmt_copy(FockO,SScr2);
    dmt_scale(SScr1,-1.0);
    dmt_sum(SScr1,SScr2);
    dmt_scale(SScr2,-1.0);
  }

 /* now form SScr1 = F + H */
  dmt_copy(Hcore,SScr1);
  dmt_sum(Fock,SScr1);
  
  nlocal = dmt_nlocal(SScr1);

  neelec=0.0;
  for (nl=0; nl < nlocal ; nl++) {
    dmt_get_block_dsc(SScr1,nl,&ib,&isz,&jb,&jsz,&fblk);
    dmt_get_block_dsc(Pmat,nl,&ib,&isz,&jb,&jsz,&pblk);

    if (!scf_info->iopen) {
      if (ib!=jb) {
        for(i=0; i < isz*jsz ; i++) {
          neelec += pblk[i]*fblk[i];
        }
      } else {
        for (i=0; i < isz ; i++) {
          for (j=0; j <= i ; j++) {
            ij=i*jsz+j;
            neelec += pblk[ij]*fblk[ij];
          }
        }
      }
    } else {
      dmt_get_block_dsc(PmatO,nl,&ib,&isz,&jb,&jsz,&poblk);
      dmt_get_block_dsc(SScr2,nl,&ib,&isz,&jb,&jsz,&gblk);
      if (ib!=jb) {
        for (i=0; i < isz*jsz ; i++) neelec += 
                                pblk[i]*fblk[i]-poblk[i]*gblk[i];
      } else {
        for(i=0; i < isz ; i++) {
          for(j=0; j <= i ; j++) {
            ij=i*jsz+j;
            neelec += pblk[ij]*fblk[ij]-poblk[ij]*gblk[ij];
          }
        }
      }
    }
  }

  neelec *= 0.5;
  gsum0(&neelec,1,5,mtype_get(),0);
  bcast0(&neelec,sizeof(double),mtype_get(),0);

  return neelec;
}

/*************************************************************************
 *
 * input:
 *   scf_info = pointer to scf struct
 *   DPmat    = scattered dmt matrix containing the AO density diff. matrix
 *
 * return rms delta(P)
 */

GLOBAL_FUNCTION double
scf_rms_delta_density(scf_info,DPmat)
scf_struct_t *scf_info;
dmt_matrix DPmat;
{
  int i,j,ij;
  int nlocal,nl;
  int ib,jb,isz,jsz;
  double delta;
  double *dpblk;

  assert(dmt_distribution(DPmat) == SCATTERED);

  nlocal = dmt_nlocal(DPmat);

  delta=0.0;
  for (nl=0; nl < nlocal ; nl++) {
    dmt_get_block_dsc(DPmat,nl,&ib,&isz,&jb,&jsz,&dpblk);

    if (ib!=jb) {
      for(i=0; i < isz*jsz ; i++) delta += dpblk[i]*dpblk[i];
    } else {
      for(i=0; i < isz ; i++) {
        for(j=0; j <= i ; j++) {
          ij=i*jsz+j;
          delta += dpblk[ij]*dpblk[ij];
        }
      }
    }
  }

  gsum0(&delta,1,5,mtype_get(),0);
  bcast0(&delta,sizeof(double),mtype_get(),0);

  delta = sqrt(delta/scf_info->mxcoef2);

  return delta;
}

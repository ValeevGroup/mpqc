
/* $Log$
 * Revision 1.1  1993/12/29 12:53:15  etseidl
 * Initial revision
 *
 * Revision 1.3  1992/06/23  20:04:28  seidl
 * change dmt matrices to uppercase,
 * get rid of unnecessary matrice
 *
 * Revision 1.2  1992/06/17  21:54:12  jannsen
 * cleaned up for saber-c
 *
 * Revision 1.1.1.1  1992/03/17  16:26:12  seidl
 * DOE-NIH Quantum Chemistry Library 0.0
 *
 * Revision 1.1  1992/03/17  16:26:11  seidl
 * Initial revision
 *
 * Revision 1.3  1992/02/13  00:42:02  seidl
 * fix small bug for converged wavefunction
 *
 * Revision 1.2  1992/02/07  12:58:11  seidl
 * added pretty pictures, fix bug for converged wavefunction
 *
 * Revision 1.1  1992/02/04  23:48:08  seidl
 * Initial revision
 *
 * Revision 1.3  1992/01/16  19:53:07  seidl
 * use new integral routines
 *
 * Revision 1.2  1992/01/09  11:44:36  seidl
 * no longer include util/ipv2/ip_libv2.h
 *
 * Revision 1.1  1992/01/02  16:19:07  seidl
 * Initial revision
 * */

static char rcsid[]="$Id$";

#include <stdio.h>
#include <tmpl.h>
#include <math/dmt/libdmt.h>
#include <math/array/math_lib.h>
#include <math.h>
#include "scf.h"

#include "scf_fock.gbl"
#include "scf_fock.lcl"


GLOBAL_FUNCTION int
scf_make_fock(_scf_info,FOCK,FOCKO,SCF_VEC,_occ_num,SCR1,SCR2,
              SCR3,iter,_outfile)
scf_struct_t *_scf_info;
dmt_matrix FOCK;
dmt_matrix FOCKO;
dmt_matrix SCF_VEC;
double_vector_t *_occ_num;
dmt_matrix SCR1;
dmt_matrix SCR2;
dmt_matrix SCR3;
int iter;
FILE *_outfile;
{
  int i,j;
  int ib,isz,ist,jb,jsz,jst;
  int nlocalb=dmt_nlocal(FOCK);
  int nn;
  double occi, occj, occ0;
  double *Fblk,*FOblk;

 /* transform fock to mo basis */
  dmt_mult(FOCK,SCF_VEC,SCR3);
  dmt_mult(SCF_VEC,SCR3,SCR1);
  dmt_copy(SCR1,FOCK);

 /* transform fock_open to mo basis */
  dmt_mult(FOCKO,SCF_VEC,SCR3);
  dmt_mult(SCF_VEC,SCR3,SCR2);
  dmt_copy(SCR2,FOCKO);

 /* form effective fock matrix in mo basis */

  occ0 = _occ_num->d[0];
  for (nn=0; nn < nlocalb; nn++ ) {
    dmt_get_block(FOCK,nn,&ib,&jb,&Fblk);
    dmt_get_block(FOCKO,nn,&ib,&jb,&FOblk);
    dmt_describe_block(FOCK,ib,&ist,&isz);
    dmt_describe_block(FOCK,jb,&jst,&jsz);

    for(i=0; i < isz ; i++) {
      for(j=0; j < jsz ; j++) {
        occi = _occ_num->d[ist+i];
        occj = _occ_num->d[jst+j];

     /* default: Guest & Saunders general form 
            C        O         V
        ----------
        |        |
     C  |   fc   |
        |        |
        -------------------
        |        |        |
     O  | 2fc-fo |   fc   |
        |        |        |
        ----------------------------
        |        |        |        |
     V  |   fc   |   fo   |   fc   |
        |        |        |        |
        ----------------------------
      */
        if(iter < _scf_info->maxiter-1 && 
           !_scf_info->converged && !_scf_info->fock_typ) {
          if(occi == occj) 
            Fblk[i*jsz+j] = Fblk[i*jsz+j];
          else if(occi && occj)
            Fblk[i*jsz+j] = 2.0*Fblk[i*jsz+j]-FOblk[i*jsz+j];
          else if(occi==2.0 || occj==2.0)
            Fblk[i*jsz+j] = Fblk[i*jsz+j];
          else
            Fblk[i*jsz+j] = FOblk[i*jsz+j];
          }

     /* Guest & Saunders' form for high spin
            C        O         V
        ----------
        |        |
     C  | 2fc-fo |
        |        |
        -------------------
        |        |        |
     O  | 2fc-fo | 2fc-fo |
        |        |        |
        ----------------------------
        |        |        |        |
     V  |   fc   |   fo   | 2fc-fo |
        |        |        |        |
        ----------------------------
      */
        else if(iter < _scf_info->maxiter-1 && !_scf_info->converged 
                                          && _scf_info->fock_typ == 1) {
          if (occi == occj || occi && occj)
            Fblk[i*jsz+j] = 2.0*Fblk[i*jsz+j]-FOblk[i*jsz+j];
          else if(occi==2.0 || occj==2.0)
            Fblk[i*jsz+j] = Fblk[i*jsz+j];
          else
            Fblk[i*jsz+j] = FOblk[i*jsz+j];
          }

     /* test form
            C        O         V
        ----------
        |        |
     C  |   fo   |
        |        |
        -------------------
        |        |        |
     O  | 2fc-fo |   fo   |
        |        |        |
        ----------------------------
        |        |        |        |
     V  |   fc   |   fo   |   fo   |
        |        |        |        |
        ----------------------------
      */
        else if(iter < _scf_info->maxiter-1 && !_scf_info->converged &&
                                              _scf_info->fock_typ == 2) {
          if(occi == occj) 
            FOblk[i*jsz+j] = FOblk[i*jsz+j];
          else if(occi && occj)
            Fblk[i*jsz+j] = 2.0*Fblk[i*jsz+j]-FOblk[i*jsz+j];
          else if(occi==2.0 || occj==2.0)
            Fblk[i*jsz+j] = Fblk[i*jsz+j];
          else
            Fblk[i*jsz+j] = FOblk[i*jsz+j];
          }

     /* form for converged wavefunction
            C        O         V
        ----------
        |        |
     C  |   fc   |
        |        |
        -------------------
        |        |        |
     O  | 2fc-fo |   fo   |
        |        |        |
        ----------------------------
        |        |        |        |
     V  |   fc   |   fo   |   fo   |
        |        |        |        |
        ----------------------------
      */
        else {
          if((occi+occj)==4.0)
            Fblk[i*jsz+j] = Fblk[i*jsz+j];
          else if(occi == occj) 
            Fblk[i*jsz+j] = FOblk[i*jsz+j];
          else if(occi && occj)
            Fblk[i*jsz+j] = 2.0*Fblk[i*jsz+j]-FOblk[i*jsz+j];
          else if(occi==2.0 || occj==2.0)
            Fblk[i*jsz+j] = Fblk[i*jsz+j];
          else
            Fblk[i*jsz+j] = FOblk[i*jsz+j];
          }

        if(ib==jb && j==i) {
          if (occi == occ0 && occi)
            Fblk[i*jsz+j] -= _scf_info->lvl_shift;
          else 
             if (occi) Fblk[i*jsz+j] -= 0.5*_scf_info->lvl_shift;
          }
        }
      }
    }

  return(0);
  }

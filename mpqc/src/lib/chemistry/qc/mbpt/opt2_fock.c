
/* $Log$
 * Revision 1.2  1995/10/11 21:13:50  cljanss
 * ParsedKeyVal now uses iostream instead of stdio.  Various other cleanups.
 *
 * Revision 1.1  1995/10/09 04:48:53  cljanss
 * The MBPT code from the mpqcic directory has been moved here.  Ida's
 * new gradient code is here as well.
 *
 * Revision 1.2  1995/03/18  00:09:50  cljanss
 * Using util/group to provide picl support.  Deleted the comm directory.
 *
 * Revision 1.1  1995/01/12  02:50:33  ibniels
 * Added opt2 and zapt2(2) capability
 *
 * Revision 1.1.1.1  1993/12/29  12:53:16  etseidl
 * SC source tree 0.1
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

#include <stdio.h>
#include <tmpl.h>
#include <math.h>
#include <util/sgen/sgen.h>
#include <math/dmt/libdmt.h>
#include <util/misc/libmisc.h>
#include <math/array/math_lib.h>
#include <chemistry/qc/intv2/int_libv2.h>
#include <chemistry/qc/dmtsym/sym_dmt.h>
#include <chemistry/qc/dmtscf/scf_dmt.h>
#include <chemistry/qc/force/libforce.h>
#include <util/group/picl.h>

#include <chemistry/qc/mbpt/opt2_fock.h>

int
mbpt_make_opt2_fock(scf_struct_t *scf_info,
                    dmt_matrix FOCK, dmt_matrix FOCKO,
                    dmt_matrix SCF_VEC,
                    dmt_matrix SCR1, dmt_matrix SCR2, dmt_matrix SCR3)
{
  int i,j;
  int ib,isz,ist,jb,jsz,jst;
  int nlocalb=dmt_nlocal(FOCK);
  int nn;
  int errcod;
  double occi, occj;
  double *Fblk,*FOblk;
  double_vector_t occ_num;

/* set up the occupation number array */
  errcod = allocbn_double_vector(&occ_num,"n",scf_info->nbfao);
  if(errcod != 0) {
    fprintf(stderr,"could not allocate memory for occ_num vector\n");
    return(-1);
    }

  for(j=0; j < scf_info->nclosed ; j++)    occ_num.d[j]=2.0;
  for(; j < scf_info->nclosed+scf_info->nopen ; j++)
                                           occ_num.d[j]=1.0;
  for(; j < scf_info->nbfso ; j++)         occ_num.d[j]=0.0;

/* transform fock to mo basis */
  dmt_mult(FOCK,SCF_VEC,SCR3);
  dmt_mult(SCF_VEC,SCR3,SCR1);
  dmt_copy(SCR1,FOCK);

 /* transform fock_open to mo basis */
  dmt_mult(FOCKO,SCF_VEC,SCR3);
  dmt_mult(SCF_VEC,SCR3,SCR2);
  dmt_copy(SCR2,FOCKO);

 /* form effective fock matrix in mo basis */

  for (nn=0; nn < nlocalb; nn++ ) {
    dmt_get_block(FOCK,nn,&ib,&jb,&Fblk);
    dmt_get_block(FOCKO,nn,&ib,&jb,&FOblk);
    dmt_describe_block(FOCK,ib,&ist,&isz);
    dmt_describe_block(FOCK,jb,&jst,&jsz);

    for(i=0; i < isz ; i++) {
      for(j=0; j < jsz ; j++) {
        occi = occ_num.d[ist+i];
        occj = occ_num.d[jst+j];

     /* Guest & Saunders general form.
        This is the form used for an OPT2 calculation.

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

        if(occi == occj) 
          Fblk[i*jsz+j] = Fblk[i*jsz+j];
        else if(occi && occj)
          Fblk[i*jsz+j] = 2.0*Fblk[i*jsz+j]-FOblk[i*jsz+j];
        else if(occi==2.0 || occj==2.0)
          Fblk[i*jsz+j] = Fblk[i*jsz+j];
        else
          Fblk[i*jsz+j] = FOblk[i*jsz+j];
        }
      }
    }

  return(0);
  }

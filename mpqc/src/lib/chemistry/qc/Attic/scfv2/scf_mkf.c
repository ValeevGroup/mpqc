
/* $Log$
 * Revision 1.1  1993/12/29 12:53:04  etseidl
 * Initial revision
 *
 * Revision 1.1.1.1  1992/03/17  17:09:00  seidl
 * DOE-NIH Quantum Chemistry Library 0.0
 *
 * Revision 1.1  1992/03/17  17:08:59  seidl
 * Initial revision
 *
 * Revision 1.2  1992/01/16  19:53:07  seidl
 * use new integral routines
 *
 * Revision 1.1  1991/12/20  16:23:11  seidl
 * Initial revision
 *
 * Revision 1.1  1991/09/23  00:43:55  seidl
 * Initial revision
 * */

static char rcsid[]="$Id$";

#include <stdio.h>
#include <tmpl.h>
#include <math/array/math_lib.h>
#include <math.h>
#include <chemistry/qc/intv2/int_libv2.h>
#include <util/ipv2/ip_libv2.h>
#include "scf_lib.h"
#include "alloc_scf_lib.h"
#include "free_scf_lib.h"
#include "init_scf_lib.h"
#include "ip_scf_lib.h"
#include "print_scf_lib.h"

#include "scf_make_fock.global"
#include "scf_make_fock.local"


GLOBAL_FUNCTION VOID
scf_make_fock(infile,outfile,_irreps,_scf_info,fock_c,fock_ct,fock_o,iter)
FILE *infile;
FILE *outfile;
scf_irreps_t *_irreps;
scf_struct_t *_scf_info;
double **fock_c;
double **fock_ct;
double **fock_o;
int iter;

{
  int i,j,ij,m;
  int nn;
  double occi, occj, occ0;
  double **scr;
  double dampd=1.0;
  double dampo=1.0;
  scf_irrep_t *s;

  scr = (double **) init_matrix(_scf_info->nsomax,_scf_info->nsomax);

  for (m=0; m < _irreps->num_irrep ; m++) {
    s = &_irreps->irrep[m];
    if (nn=s->num_so) {

   /* transform fock to mo basis */
      tri_to_sq(s->fock,fock_ct,nn);
      mxmb(s->scf_vector,nn,1,fock_ct,1,nn,scr,1,nn,nn,nn,nn);
      mxmb(scr,1,nn,s->scf_vector,1,nn,fock_c,1,nn,nn,nn,nn);

   /* transform fock_open to mo basis */
      tri_to_sq(s->op.fock_open,fock_ct,nn);
      mxmb(s->scf_vector,nn,1,fock_ct,1,nn,scr,1,nn,nn,nn,nn);
      mxmb(scr,1,nn,s->scf_vector,1,nn,fock_o,1,nn,nn,nn,nn);

   /* form effective fock matrix in mo basis */

      ij=0;
      occ0 = s->occ_num[0];
      for (i=0; i < nn; i++ ) {
        for (j=0; j <= i; j++) {
          occi = s->occ_num[i];
          occj = s->occ_num[j];

     /* default: Guest & Saunders general form */
          if(iter < _scf_info->maxiter-1 && 
             !_scf_info->converged && 
             !_scf_info->fock_typ) {
            if(occi == occj) 
              s->op.fock_eff[ij] = fock_c[i][j];
            else if(occi)
              s->op.fock_eff[ij] = 2.0*fock_c[i][j]-fock_o[i][j];
            else {
              if(occj==2.0)
                s->op.fock_eff[ij] = fock_c[i][j];
              else
                s->op.fock_eff[ij] = fock_o[i][j];
              }
            }

     /* Guest & Saunders' form for high spin */
          else if(iter < _scf_info->maxiter-1 && !_scf_info->converged 
                                            && _scf_info->fock_typ == 1) {
            if (occi == occj || occi)
              s->op.fock_eff[ij] = 2.0*fock_c[i][j]-fock_o[i][j];
            else if (occj == 2.0) s->op.fock_eff[ij] = fock_c[i][j];
            else s->op.fock_eff[ij] = fock_o[i][j];
            }

     /* test form (fo fo fo) */
          else if(iter < _scf_info->maxiter-1 && !_scf_info->converged &&
                                                _scf_info->fock_typ == 2) {
            if (occi == occj) s->op.fock_eff[ij] = fock_o[i][j];
            else if(occi)
              s->op.fock_eff[ij] = 2.0*fock_c[i][j]-fock_o[i][j];
            else if(occj == 2.0)
              s->op.fock_eff[ij] = fock_c[i][j];
            else s->op.fock_eff[ij] = fock_o[i][j];
            }

     /* test form a*(fc fc fc) */
          else if(iter < _scf_info->maxiter-1 && !_scf_info->converged && 
                                                 _scf_info->fock_typ == 3) {
            if (occi == occj) s->op.fock_eff[ij] = dampd*fock_c[i][j];
            else if(occi)
              s->op.fock_eff[ij] = 
                dampo*(2.0*fock_c[i][j]-fock_o[i][j]);
            else if(occj == 2.0)
              s->op.fock_eff[ij] = dampo*fock_c[i][j];
            else s->op.fock_eff[ij] = dampo*fock_o[i][j];
            }

      /* test form a*(fo fo fo) */
          else if(iter < _scf_info->maxiter-1 && !_scf_info->converged && 
                                                _scf_info->fock_typ == 4) {
            if (occi == occj) s->op.fock_eff[ij] = dampd*fock_o[i][j];
            else if(occi)
              s->op.fock_eff[ij] = 2.0*fock_c[i][j]-fock_o[i][j];
            else if(occj == 2.0)
              s->op.fock_eff[ij] = fock_c[i][j];
            else s->op.fock_eff[ij] = fock_o[i][j];
            }

     /* test form a*(2fc-fo 2fc-fo 2fc-fo) */
          else if(iter < _scf_info->maxiter-1 && !_scf_info->converged && 
                                                  _scf_info->fock_typ == 5) {
            if (occi == occj)
              s->op.fock_eff[ij] = 
                dampd*(2.0*fock_c[i][j]-fock_o[i][j]);
            else if(occi)
              s->op.fock_eff[ij] = 2.0*fock_c[i][j]-fock_o[i][j];
            else if(occj == 2.0)
              s->op.fock_eff[ij] = fock_c[i][j];
              else s->op.fock_eff[ij] = fock_o[i][j];
              }

     /* form for converged wavefunction */
          else {
            if (occi == 2.0) s->op.fock_eff[ij]=fock_c[i][j];
            else if(occj != 2.0)
              s->op.fock_eff[ij]=fock_o[i][j];
            else {
              if(occi)
                s->op.fock_eff[ij] = 2.0*fock_c[i][j]-fock_o[i][j];
              else s->op.fock_eff[ij]=fock_c[i][j];
              }
            }

          if(j==i) {
            if (occi == occ0 && occi)
               s->op.fock_eff[ij] -= _scf_info->lvl_shift;
            else 
               if (occi) s->op.fock_eff[ij] -= 0.5*_scf_info->lvl_shift;
            }
          ij++;
          }
        }
      }
    }
  free_matrix(scr,_scf_info->nsomax);
  }

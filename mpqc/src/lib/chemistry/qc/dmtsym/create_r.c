
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <tmpl.h>
#include <util/misc/libmisc.h>
#include <math/array/math_lib.h>
#include <math/dmt/libdmt.h>
#include <chemistry/qc/intv2/int_libv2.h>

#include "symm.h"
#include "symm_mac.h"

#include "create_r.gbl"
#include "create_r.lcl"

/************************************************************************
 *	
 * this creates the full basis set transformation matrix R for a given
 * symmetry operation g
 *
 * input:
 *   centers = pointer to centers struct
 *   sym_info = pointer to symmetry struct
 *   r = dmt column distributed matrix
 *   g = the index of the symmetry operation
 *
 * on return:
 *   r contains R
 *
 * return 0 on success, -1 on failure
 */

GLOBAL_FUNCTION int
sym_create_r(centers,sym_info,r,g)
centers_t *centers;
sym_struct_t *sym_info;
dmt_matrix r;
int g;
{
  int i,j;
  int atom,shell,bfunc,ffunc,gc;
  int gffunc,glfunc,msize,nlocal,lffunc;
  center_t *c;
  int nat = centers->n;
  double **rp = sym_info->Rp[g];
  double **rd = sym_info->Rd[g];
  double **rf=0;
  double **rg=0;
  double *rcol,**lr;

  assert(dmt_distribution(r) == COLUMNS);

  if (sym_info->Rf != NULL) rf = sym_info->Rf[g];
  if (sym_info->Rg != NULL) rg = sym_info->Rg[g];

  dmt_fill(r,0.0);
  dmt_get_col(r,0,&gffunc,&rcol);

  nlocal = dmt_nlocal(r);
  msize = dmt_size(r);

  glfunc = gffunc+nlocal;

  lr = (double **) malloc(sizeof(double *)*nlocal);
  if (lr==NULL) return(-1);

  for (i=0; i < nlocal ; i++) lr[i] = &rcol[i*msize];

  lffunc=ffunc=0;
  for (atom=0; atom < nat ; atom++) {
    c = &centers->center[atom];
    bfunc = sym_info->first[sym_info->atom_map[atom][g]];

    for (shell=0; shell < c->basis.n ; shell++) {
      for (gc=0; gc < c->basis.shell[shell].ncon ; gc++) {
        switch(c->basis.shell[shell].type[gc].am) {
        case _AM_S:
          if (ffunc >= gffunc && ffunc < glfunc) {
            lr[lffunc][bfunc]=1.0;
            lffunc++;
          }
          bfunc++; ffunc++;
          break;

        case _AM_P:
          for (i=0; i < 3 ; i++) {
            for (j=0; j < 3 ; j++) {
              if ((ffunc+i) >= gffunc && (ffunc+i) < glfunc) {
                lr[lffunc][bfunc+j] = rp[i][j];
              }
            }
            if ((ffunc+i) >= gffunc && (ffunc+i) < glfunc) lffunc++;
          }
          bfunc+=3;
          ffunc+=3;
          break;

        case _AM_D:
          for (i=0; i < 6 ; i++) {
            for (j=0; j < 6 ; j++) {
              if ((ffunc+i) >= gffunc && (ffunc+i) < glfunc) {
                lr[lffunc][bfunc+j] = rd[i][j];
              }
            }
            if((ffunc+i) >= gffunc && (ffunc+i) < glfunc) lffunc++;
          }
          bfunc+=6;
          ffunc+=6;
          break;

        case _AM_F:
          for (i=0; i < 10 ; i++) {
            for (j=0; j < 10 ; j++) {
              if ((ffunc+i) >= gffunc && (ffunc+i) < glfunc) {
                lr[lffunc][bfunc+j] = rf[i][j];
              }
            }
            if ((ffunc+i) >= gffunc && (ffunc+i) < glfunc) lffunc++;
          }
          bfunc+=10;
          ffunc+=10;
          break;

        case _AM_G:
          for (i=0; i < 15 ; i++) {
            for (j=0; j < 15 ; j++) {
              if ((ffunc+i) >= gffunc && (ffunc+i) < glfunc) {
                lr[lffunc][bfunc+j] = rg[i][j];
              }
            }
            if ((ffunc+i) >= gffunc && (ffunc+i) < glfunc) lffunc++;
          }
          bfunc+=15;
          ffunc+=15;
          break;

        default:
          fprintf(stderr,"sym_create_r: cannot yet handle shells with am = %d",
                          c->basis.shell[shell].type[gc].am);
          return -1;
        }
      }
    }
  }

  free(lr);
  return 0;
}

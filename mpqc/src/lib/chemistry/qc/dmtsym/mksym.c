
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <tmpl.h>
#include <math/array/math_lib.h>
#include <chemistry/qc/intv2/int_libv2.h>
#include <math/dmt/libdmt.h>

#include "symm.h"
#include "symm_mac.h"

#include "mksym.gbl"
#include "mksym.lcl"

/* ok, here's the plan.  given a skeleton matrix M, we want to form
 * the symmetrized matrix M'. According to D&K, M' = 1/g* [sum(g) R~ * M * R].
 * We can rewrite each of the g transformations as
 *  M' = (R~*M)*R => M'~ = R~*(R~*M)~, and since M and M' are symmetric,
 * we have M' = R~*(R~*M)~.  Thus our strategy is to form scr=R~*M, 
 * form the transpose of scr, and then form M'=R~*scr.  We could avoid doing
 * the transposition (the most time consuming step) by putting (R~*M)[i][j]
 * in scr[j][i], but since we cannot guarantee that the column boundaries
 * will fall between shells, this can't be done.
 */

GLOBAL_FUNCTION VOID
sym_sym_matrix(centers,sym_info,skel,sym)
centers_t *centers;
sym_struct_t *sym_info;
dmt_matrix skel;
dmt_matrix sym;
{
  int i;
  int ng;
  int nlocal,msize;
  dmt_matrix scr;
  dmt_matrix old_skel=-1,old_sym=-1;
  double *symcol,*skelcol,*scrcol;
  double **lsym,**lskel,**lscr;

 /* this only works for matrices distributed by columns.  if skel or sym
  * are distributed by blocks, then use temporary columns matrices */

  if (dmt_distribution(skel) != COLUMNS) {
    old_skel=skel;
    skel = dmt_columns("sym_sym_matrix: scratch column skel matrix",old_skel);
  }

  if (dmt_distribution(sym) != COLUMNS) {
    old_sym=sym;
    sym = dmt_columns("sym_sym_matrix: scratch column sym matrix",old_sym);
  }

  nlocal = dmt_nlocal(sym);
  msize = dmt_size(sym);

  scr=dmt_create("sym_sym_matrix: scr 1",msize,COLUMNS);

  dmt_fill(sym,0.0);
  dmt_copy(skel,sym);

  dmt_get_col(sym,0,&i,&symcol);
  dmt_get_col(skel,0,&i,&skelcol);
  dmt_get_col(scr,0,&i,&scrcol);

  lsym = (double **) malloc(sizeof(double *)*nlocal);
  if (lsym==NULL) return;
  lskel = (double **) malloc(sizeof(double *)*nlocal);
  if (lskel==NULL) return;
  lscr = (double **) malloc(sizeof(double *)*nlocal);
  if (lscr==NULL) return;

  for (i=0; i < nlocal ; i++) {
    lskel[i] = &skelcol[i*msize];
    lsym[i] = &symcol[i*msize];
    lscr[i] = &scrcol[i*msize];
  }

  for (ng=1; ng < sym_info->g ; ng++) {
    dmt_fill(scr,0.0);
    rtrans_mult(centers,sym_info,lskel,lscr,nlocal,ng);
    dmt_transpose(scr);
    rtrans_mult(centers,sym_info,lscr,lsym,nlocal,ng);
  }

  dmt_scale(sym,(double) 1.0/sym_info->g);

  dmt_free(scr);
  free(lsym);
  free(lskel);
  free(lscr);

  if (old_skel > -1) dmt_free(skel);
  if (old_sym > -1) {
    dmt_copy(sym,old_sym);
    dmt_free(sym);
  }
}


LOCAL_FUNCTION int
rtrans_mult(centers,sym_info,lskel,lsym,nlocal,g)
centers_t *centers;
sym_struct_t *sym_info;
double **lskel;
double **lsym;
int nlocal;
int g;
{
  int i,j,k;
  int atom,shell,bfunc,ffunc,gc;
  center_t *c;
  int nat = centers->n;
  double **rp = sym_info->Rp[g];
  double **rd = sym_info->Rd[g];
  double **rf;
  double **rg;
  char errmsg[81];

  if (sym_info->Rf != NULL) rf = sym_info->Rf[g];
  if (sym_info->Rg != NULL) rg = sym_info->Rg[g];

  ffunc=0;
  for (atom=0; atom < nat ; atom++) {
    c = &centers->center[atom];
    bfunc = sym_info->first[sym_info->atom_map[atom][g]];

    for (shell=0; shell < c->basis.n ; shell++) {
      for (gc=0; gc < c->basis.shell[shell].ncon ; gc++) {
        switch(c->basis.shell[shell].type[gc].am) {
        case _AM_S:
          for (j=0; j < nlocal ; j++) lsym[j][ffunc]+=lskel[j][bfunc];
          bfunc++; ffunc++;
          break;

        case _AM_P:
          for (j=0; j < nlocal ; j++) {
            for (i=0; i < 3 ; i++) {
              for (k=0; k < 3 ; k++) 
                lsym[j][ffunc+i] += rp[i][k]*lskel[j][bfunc+k];
            }
          }
          bfunc+=3;
          ffunc+=3;
          break;

        case _AM_D:
          for (j=0; j < nlocal ; j++) {
            for (i=0; i < 6 ; i++) {
              for (k=0; k < 6 ; k++) 
                lsym[j][ffunc+i] += rd[i][k]*lskel[j][bfunc+k];
            }
          }
          bfunc+=6;
          ffunc+=6;
          break;

        case _AM_F:
          for (j=0; j < nlocal ; j++) {
            for (i=0; i < 10 ; i++) {
              for (k=0; k < 10 ; k++) 
                lsym[j][ffunc+i] += rf[i][k]*lskel[j][bfunc+k];
            }
          }
          bfunc+=10;
          ffunc+=10;
          break;

        case _AM_G:
          for (j=0; j < nlocal ; j++) {
            for (i=0; i < 15 ; i++) {
              for (k=0; k < 15 ; k++) 
                lsym[j][ffunc+i] += rg[i][k]*lskel[j][bfunc+k];
            }
          }
          bfunc+=15;
          ffunc+=15;
          break;

        default:
          fprintf(stderr,"rtrans_mult: cannot yet handle shells with am = %d",
            c->basis.shell[shell].type[gc].am);
          return -1;
        }
      }
    }
  }

  return 0;
}

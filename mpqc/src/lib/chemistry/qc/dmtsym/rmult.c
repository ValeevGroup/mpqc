
/* $Log$
 * Revision 1.1  1993/12/29 12:52:58  etseidl
 * Initial revision
 *
 * Revision 1.1.1.1  1992/03/17  16:27:31  seidl
 * DOE-NIH Quantum Chemistry Library 0.0
 *
 * Revision 1.1  1992/03/17  16:27:29  seidl
 * Initial revision
 *
 * Revision 1.2  1992/01/27  16:51:28  seidl
 * libint needs libmath
 *
 * Revision 1.1  1992/01/27  16:36:57  seidl
 * Initial revision
 * */

static char rcsid[] = "$Id$";

#include <stdio.h>
#include <math.h>
#include <tmpl.h>
#include <math/array/math_lib.h>
#include <chemistry/qc/intv2/int_libv2.h>
#include <math/dmt/libdmt.h>

#include "symm.h"
#include "symm_mac.h"
#include "symerr.gbl"

#include "rmult.gbl"
#include "rmult.lcl"

GLOBAL_FUNCTION int
sym_rtrans_mult(_centers,_sym_info,lskel,lsym,nlocal,g,_outfile)
centers_t *_centers;
sym_struct_t *_sym_info;
double **lskel;
double **lsym;
int nlocal;
int g;
FILE *_outfile;
{
  int i,j,k;
  int atom,shell,count,bfunc,ffunc,gc;
  center_t *c;
  int **s_m = _sym_info->shell_map;
  int **a_m = _sym_info->atom_map;
  int nat = _centers->n;
  double **rp = _sym_info->Rp[g];
  double **rd = _sym_info->Rd[g];
  double **rf;
  double **rg;
  char errmsg[81];

  if(_sym_info->Rf != NULL) rf = _sym_info->Rf[g];
  if(_sym_info->Rg != NULL) rg = _sym_info->Rg[g];

  ffunc=0;
  for(atom=count=0; atom < nat ; atom++) {
    c = &_centers->center[atom];
    bfunc = _sym_info->first[_sym_info->atom_map[atom][g]];

    for(shell=0; shell < c->basis.n ; shell++) {
      for(gc=0; gc < c->basis.shell[shell].ncon ; gc++) {
        switch(c->basis.shell[shell].type[gc].am) {
        case _AM_S:
          for(j=0; j < nlocal ; j++) lsym[j][ffunc]+=lskel[j][bfunc];
          bfunc++; ffunc++;
          break;

        case _AM_P:
          for(j=0; j < nlocal ; j++) {
            for(i=0; i < 3 ; i++) {
              for(k=0; k < 3 ; k++) 
                lsym[j][ffunc+i] += rp[i][k]*lskel[j][bfunc+k];
              }
            }
          bfunc+=3;
          ffunc+=3;
          break;

        case _AM_D:
          for(j=0; j < nlocal ; j++) {
            for(i=0; i < 6 ; i++) {
              for(k=0; k < 6 ; k++) 
                lsym[j][ffunc+i] += rd[i][k]*lskel[j][bfunc+k];
              }
            }
          bfunc+=6;
          ffunc+=6;
          break;

        case _AM_F:
          for(j=0; j < nlocal ; j++) {
            for(i=0; i < 10 ; i++) {
              for(k=0; k < 10 ; k++) 
                lsym[j][ffunc+i] += rf[i][k]*lskel[j][bfunc+k];
              }
            }
          bfunc+=10;
          ffunc+=10;
          break;

        case _AM_G:
          for(j=0; j < nlocal ; j++) {
            for(i=0; i < 15 ; i++) {
              for(k=0; k < 15 ; k++) 
                lsym[j][ffunc+i] += rg[i][k]*lskel[j][bfunc+k];
              }
            }
          bfunc+=15;
          ffunc+=15;
          break;

        default:
          sprintf(errmsg,"cannot yet handle shells with am = %d",
            c->basis.shell[shell].type[gc].am);
          serror(_outfile,__FILE__,errmsg,__LINE__);
          return(-1);
          }
        }
      }
    }

  return(0);
  }

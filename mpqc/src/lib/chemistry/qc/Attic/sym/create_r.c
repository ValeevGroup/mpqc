
/* $Log$
 * Revision 1.1  1993/12/29 12:53:13  etseidl
 * Initial revision
 *
 * Revision 1.1.1.1  1992/03/17  17:10:48  seidl
 * DOE-NIH Quantum Chemistry Library 0.0
 *
 * Revision 1.1  1992/03/17  17:10:46  seidl
 * Initial revision
 *
 * Revision 1.2  1992/01/21  11:34:46  seidl
 * use libintv2
 *
 * Revision 1.1  1991/12/20  15:45:23  seidl
 * Initial revision
 *
 * Revision 1.1  1991/12/03  16:41:54  etseidl
 * Initial revision
 *
 * Revision 1.2  1991/12/02  19:54:05  seidl
 * *** empty log message ***
 *
 * Revision 1.1  1991/11/22  18:28:40  seidl
 * Initial revision
 * */

static char rcsid[] = "$Id$";

#include <stdio.h>
#include <math.h>
#include <tmpl.h>
#include <math/array/math_lib.h>
#include <util/ipv2/ip_libv2.h>
#include <chemistry/qc/intv2/int_libv2.h>

#include "symm.h"
#include "symm_mac.h"
#include "symerr.gbl"

#include "create_r.gbl"
#include "create_r.lcl"

GLOBAL_FUNCTION int
sym_create_r(_centers,_sym_info,r,g,_outfile)
centers_t *_centers;
sym_struct_t *_sym_info;
double_matrix_t *r;
int g;
FILE *_outfile;
{
  int i,j,ij;
  int i1,il,j1,jl;
  int atom,shell,count,bfunc,ffunc,gc;
  int ng;
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

  for(i=0; i < r->n1 ; i++) bzero(r->d[i],sizeof(double)*r->n2);
  
  ffunc=0;
  for(atom=count=0; atom < nat ; atom++) {
    c = &_centers->center[atom];
    bfunc = _sym_info->first[_sym_info->atom_map[atom][g]];
    for(shell=0; shell < c->basis.n ; shell++) {
      for(gc=0; gc < c->basis.shell[shell].ncon ; gc++) {
        switch(c->basis.shell[shell].type[gc].am) {
        case _AM_S:
          r->d[ffunc][bfunc]=1.0;
          bfunc++; ffunc++;
          break;
        case _AM_P:
          for(i=0; i < 3 ; i++)
            for(j=0; j < 3 ; j++) {
              r->d[ffunc+i][bfunc+j] = rp[i][j];
              }
          bfunc+=3;
          ffunc+=3;
          break;
        case _AM_D:
          for(i=0; i < 6 ; i++)
            for(j=0; j < 6 ; j++) {
              r->d[ffunc+i][bfunc+j] = rd[i][j];
              }
          bfunc+=6;
          ffunc+=6;
          break;
        case _AM_F:
          for(i=0; i < 10 ; i++)
            for(j=0; j < 10 ; j++) {
              r->d[ffunc+i][bfunc+j] = rf[i][j];
              }
          bfunc+=10;
          ffunc+=10;
          break;
        case _AM_G:
          for(i=0; i < 15 ; i++)
            for(j=0; j < 15 ; j++) {
              r->d[ffunc+i][bfunc+j] = rg[i][j];
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

GLOBAL_FUNCTION int
sym_create_rp(_centers,_sym_info,r,g,_outfile)
centers_t *_centers;
sym_struct_t *_sym_info;
double_matrix_t *r;
int g;
FILE *_outfile;
{
  int i,j,k,l,f,fp;
  int atom,ffunc,bfunc;
  int shell;
  int nfunc1,nfunc2,neqat;
  int errcod;
  int nbf = _centers->nfunc;
  int nat = _centers->n;
  int gc;

  char errmsg[81];

  center_t *c;

  for(i=0; i < r->n1 ; i++) bzero(r->d[i],sizeof(double)*r->n2);

  ffunc=0;
  for(atom=0; atom < nat ; atom++) {
    c = &_centers->center[atom];
    bfunc = _sym_info->first[_sym_info->atom_map[atom][g]];
    for(shell=0; shell < c->basis.n ; shell++) {
      for(gc=0; gc < c->basis.shell[shell].ncon ; gc++) {
        switch(c->basis.shell[shell].type[gc].am) {
        case _AM_S:
          r->d[bfunc][ffunc]=1.0;
          bfunc++;ffunc++;
          break;
        case _AM_P:
          for(i=0; i < 3 ; i++)
            for(j=0; j < 3 ; j++) {
              r->d[bfunc+i][ffunc+j] = _sym_info->Rp[g][j][i];
              }
          bfunc+=3;
          ffunc+=3;
          break;
        case _AM_D:
          for(i=0; i < 6 ; i++)
            for(j=0; j < 5 ; j++) {
              r->d[bfunc+i][ffunc+j] = _sym_info->Rdp[g][i][j];
              }
          bfunc+=6;
          ffunc+=5;
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
  }

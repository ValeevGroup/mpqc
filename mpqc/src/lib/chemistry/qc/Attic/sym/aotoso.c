
/* $Log$
 * Revision 1.1  1993/12/29 12:53:13  etseidl
 * Initial revision
 *
 * Revision 1.1.1.1  1992/03/17  17:10:44  seidl
 * DOE-NIH Quantum Chemistry Library 0.0
 *
 * Revision 1.1  1992/03/17  17:10:43  seidl
 * Initial revision
 *
 * Revision 1.3  1992/01/21  11:35:02  seidl
 * use libintv2
 *
 * Revision 1.2  1992/01/02  12:36:55  seidl
 * *** empty log message ***
 *
 * Revision 1.1  1991/12/20  15:45:23  seidl
 * Initial revision
 *
 * Revision 1.1  1991/12/03  16:41:01  etseidl
 * Initial revision
 *
 * Revision 1.3  1991/12/02  19:53:49  seidl
 * *** empty log message ***
 *
 * Revision 1.2  1991/11/25  21:42:56  seidl
 * attempt to use pure am d functions was a disaster
 *
 * Revision 1.1  1991/11/22  18:29:34  seidl
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

#include "aotoso.gbl"
#include "aotoso.lcl"

static double t1[1][1] = {1.0};
static double t2[2][2] = {0.70710678118654752440,  0.70710678118654752440,
                          0.70710678118654752440, -0.70710678118654752440};
static double t3[3][3] =
 {0.57735026918962576451,  0.81649658092772603273,  0.0,
  0.57735026918962576451, -0.40824829046386301636, -0.70710678118654752440,
  0.57735026918962576451, -0.40824829046386301636,  0.70710678118654752440};

static double t4[4][4] =
 { 0.5,  0.5,  0.5,  0.5,
   0.5,  0.5, -0.5, -0.5,
   0.5, -0.5, -0.5,  0.5,
   0.5, -0.5,  0.5, -0.5};

GLOBAL_FUNCTION int
sym_aotoso(_cent,_s_info,_irr,_ats,_sta,_outfile)
centers_t *_cent;
sym_struct_t *_s_info;
scf_irreps_t *_irr;
double_matrix_t *_ats;
double_matrix_t *_sta;
FILE *_outfile;
{
  int errcod;
  int nba = _cent->nfunc;
  double_matrix_t chars;
  double_matrix_t scr1;
  double_matrix_t scr2;

  errcod = allocbn_double_matrix(&chars,"n1 n2",nba,_s_info->g);
  if(errcod != 0) {
    fprintf(_outfile,"sym_aotoso: could not allocate memory for chars\n");
    return(-1);
    }
  errcod = allocbn_double_matrix(&scr1,"n1 n2",nba,nba);
  if(errcod != 0) {
    fprintf(_outfile,"sym_aotoso: could not allocate memory for scr2\n");
    return(-1);
    }
  errcod = allocbn_double_matrix(&scr2,"n1 n2",nba,nba);
  if(errcod != 0) {
    fprintf(_outfile,"sym_aotoso: could not allocate memory for scr2\n");
    return(-1);
    }

  errcod = _aotoso(_cent,_s_info,_ats,0,_outfile);
  if(errcod != 0) return(errcod);
    
  errcod = _aotoso(_cent,_s_info,_sta,1,_outfile);
  if(errcod != 0) return(errcod);

#if 0
  fprintf(_outfile,"ats\n");
  math_print_dm(_outfile,_ats);
  fprintf(_outfile,"sta\n");
  math_print_dm(_outfile,_sta);
#endif

#if 1
  errcod = sort_aoso(_cent,_s_info,_irr,_ats,_sta,&chars,&scr1,&scr2,_outfile);
  if(errcod != 0) return(errcod);
#endif

  free_double_matrix(&scr1);
  free_double_matrix(&scr2);
  free_double_matrix(&chars);

  return 0;
  }

LOCAL_FUNCTION int
_aotoso(_cent,_s_info,ats,_inv,_outfile)
centers_t *_cent;
sym_struct_t *_s_info;
double_matrix_t *ats;
int _inv;
FILE *_outfile;
{
  int i,j,k,l,f,fp;
  int atom,ffunc1,ffunc2;
  int nfunc1,nfunc2,neqat;
  int errcod;
  int nba = _cent->nfunc;
  int nat = _cent->n;
  double *tr[15];
  char errmsg[81];
  center_t *c;

  for(i=0; i < ats->n1 ; i++) bzero(ats->d[i],sizeof(double)*ats->n2);

  ffunc1=ffunc2=0;
  for(atom=0; atom < nat ; ) {
    c = &_cent->center[atom];
    nfunc1 = n_func_atom(c,0);
    nfunc2 = n_func_atom(c,0);
    neqat = n_eq_atoms(_s_info,atom);

    switch(neqat) {
    case 1:
      tr[0]=t1[0];
      errcod =
        make_block(_cent,c,_s_info,1,nfunc1,nfunc2,ffunc1,ffunc2,ats,tr,_inv);
      if(errcod != 0) {
        sprintf(errmsg,"cannot yet handle f functions");
        serror(_outfile,__FILE__,errmsg,__LINE__);
        return(-1);
        }
      break;

    case 2:
      tr[0]=t2[0];
      tr[1]=t2[1];
      errcod =
        make_block(_cent,c,_s_info,2,nfunc1,nfunc2,ffunc1,ffunc2,ats,tr,_inv);
      if(errcod != 0) {
        sprintf(errmsg,"cannot yet handle f functions");
        serror(_outfile,__FILE__,errmsg,__LINE__);
        return(-1);
        }
      break;
    case 3:
      tr[0]=t3[0];
      tr[1]=t3[1];
      tr[2]=t3[2];
      errcod =
        make_block(_cent,c,_s_info,3,nfunc1,nfunc2,ffunc1,ffunc2,ats,tr,_inv);
      if(errcod != 0) {
        sprintf(errmsg,"cannot yet handle f functions");
        serror(_outfile,__FILE__,errmsg,__LINE__);
        return(-1);
        }
      break;
    case 4:
      tr[0]=t4[0];
      tr[1]=t4[1];
      tr[2]=t4[2];
      tr[3]=t4[3];
      errcod =
        make_block(_cent,c,_s_info,4,nfunc1,nfunc2,ffunc1,ffunc2,ats,tr,_inv);
      if(errcod != 0) {
        sprintf(errmsg,"cannot yet handle f functions");
        serror(_outfile,__FILE__,errmsg,__LINE__);
        return(-1);
        }
      break;
    default:
      sprintf(errmsg,"cannot yet handle %d equivalent atoms", neqat);
      serror(_outfile,__FILE__,errmsg,__LINE__);
      return(-1);
      }
    ffunc1 += neqat*nfunc1;
    ffunc2 += neqat*nfunc2;
    atom += neqat;
    }
  return(0);
  }

LOCAL_FUNCTION int
n_func_atom(c,puream)
center_t *c;
int puream;
{
  int i,gc,nf;

  nf=0;
  for(i=0; i < c->basis.n ; i++) {
    if(!puream) 
      nf += c->basis.shell[i].nfunc;
    else {
      int nfc = c->basis.shell[i].nfunc;
      for(gc=0; gc < c->basis.shell[i].ncon ; gc++) {
        switch(c->basis.shell[i].type[gc].am) {
        case _AM_D:
          nf += nfc-1;
          break;
        case _AM_F:
          nf += nfc-3;
          break;
        case _AM_G:
          nf += nfc-6;
          break;
        default:
          nf += nfc;
          }
        }
      }
    }

  return nf;
  }

LOCAL_FUNCTION int
n_eq_atoms(_s_info, atom)
sym_struct_t *_s_info;
int atom;
{
  int i,neq;

  neq=0;
  for(i=0; i < _s_info->g ; i++)
    if(_s_info->atom_map[atom][i] == atom) neq++;
  
  neq = _s_info->g/neq;

  return neq;
  }

LOCAL_FUNCTION int
make_block(cs,c,si,neq,nfunc1,nfunc2,ffunc1,ffunc2,ats,tr,inv)
centers_t *cs;
center_t *c;
sym_struct_t *si;
int neq;
int nfunc1;
int nfunc2;
int ffunc1;
int ffunc2;
double_matrix_t *ats;
double **tr;
int inv;
{
  int i,j,k,l,f,fp;
  int kp,lp;
  int errcod;
  int nam;
  int gc;
  double tmp;

  for(i=0; i < neq ; i++) {
    for(j=0; j < neq ; j++) {
      for(f=fp=0; f < c->basis.n ; f++) {
        for(gc=0; gc < c->basis.shell[f].ncon ; gc++) {
          nam = c->basis.shell[f].type[gc].am;
          switch(nam) {
          case _AM_S:
            ats->d[ffunc1+i*nfunc1+fp][ffunc2+j*nfunc2+fp] = tr[i][j];
            break;
          case _AM_P:
            for(k=0; k < 3 ; k++)
              for(l=0; l < 3 ; l++)
                ats->d[ffunc1+i*nfunc1+fp+k][ffunc2+j*nfunc2+fp+l]=
                  tr[i][j]*si->Rp[i][l][k];
            break;
          case _AM_D:
            kp=ffunc1+i*nfunc1+fp;
            lp=ffunc2+j*nfunc2+fp;
            for(k=0; k < 6 ; k++) {
              for(l=0; l < 6 ; l++) {
                ats->d[kp+k][lp+l] = tr[i][j]*si->Rd[i][l][k];
                }
              }
            if(neq==1) {
              ats->d[kp+yy][lp+yy] *= sqrt(0.5);
              ats->d[kp+xx][lp+xx] *= sqrt(0.5);
              ats->d[kp+yy][lp+xx] = -ats->d[kp+yy][lp+yy];
              ats->d[kp+xx][lp+yy] = ats->d[kp+yy][lp+yy];
              }
            if(inv) {
              tmp = ats->d[kp+yy][lp+xy];
              ats->d[kp+yy][lp+xy] = ats->d[kp+xy][lp+xx];
              ats->d[kp+xy][lp+xx] = tmp;
  
              tmp = ats->d[kp+xx][lp+xy];
              ats->d[kp+xx][lp+xy] = ats->d[kp+xy][lp+yy];
              ats->d[kp+xy][lp+yy] = tmp;
              }
            break;
          default:
            return(-1);
            }
          fp+=INT_NCART(nam);
          }
        }
      }
    }
  return(0);
  }

LOCAL_FUNCTION int
sort_aoso(_cent,_s_info,_irr,ats,sta,chars,scr1,scr2,_outfile)
centers_t *_cent;
sym_struct_t *_s_info;
scf_irreps_t *_irr;
double_matrix_t *ats;
double_matrix_t *sta;
double_matrix_t *chars;
double_matrix_t *scr1;
double_matrix_t *scr2;
FILE *_outfile;
{
  int i,j,ij,k;
  int nba = _cent->nfunc;
  int errcod;
  double deg;

  for(i=0; i < _s_info->g ; i++) {
    errcod = sym_create_r(_cent,_s_info,scr1,i,_outfile);
    if(errcod != 0) {
      fprintf(_outfile,"sym_aotoso: sort_aoso:\n");
      fprintf(_outfile,"trouble in sym_create_r\n");
      return(-1);
      }

    math_dmxdm_dm(sta,1,scr1,1,scr2,0,nba,nba,nba,0);
    math_dmxdm_dm(scr2,0,ats,0,scr1,0,nba,nba,nba,0);

    if(i) for(j=0; j < nba ; j++) chars->d[j][i] = scr1->d[j][j]/chars->d[j][0];
    else for(j=0; j < nba ; j++) chars->d[j][i] = scr1->d[j][j];
    }

  for(i=0; i < nba ; i++) chars->d[i][0] = 1.0;

  ij=0;
  for(j=0; j < _s_info->nirrep ; j++) {
    deg = (double) _irr->ir[j].degeneracy;
    for(i=0; i < nba ; i++) {
      if(comp_char(&_s_info->ct,chars->d[i],j,deg)) {
        for(k=0; k < nba ; k++) {
          scr1->d[k][ij]=ats->d[k][i];
          scr2->d[k][ij]=sta->d[k][i];
          }
        ij++;
        }
      }
    }
  
  for(i=0; i < nba ; i++) {
    for(j=0; j < nba ; j++) {
      ats->d[i][j]=scr1->d[i][j];
      sta->d[i][j]=scr2->d[i][j];
      }
    }
  return(0);
  }

LOCAL_FUNCTION int
comp_char(ct, c, irrep, deg)
char_tab_t *ct;
double *c;
int irrep;
double deg;
{
  int ng;
  double tol=1.0e-10;


  for(ng=0; ng < ct->gamma[irrep].g ; ng++)
    if(fabs(deg*c[ng]-(ct->gamma[irrep].rep[ng])) > tol) return(0);
  return(1);
  }

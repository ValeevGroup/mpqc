

/* $Log$
 * Revision 1.2  1994/06/08 01:15:18  cljanss
 * Many changes.  These include: newmat7 and nihmatrix -> scmat
 * and mpqcic -> MPSCF and updated optimize stuff.
 *
 * Revision 1.1.1.1  1993/12/29  12:52:59  etseidl
 * SC source tree 0.1
 *
 * Revision 1.5  1993/04/27  23:58:12  jannsen
 * fixed compiler type conversion complaint
 *
 * Revision 1.4  1992/06/17  21:56:04  jannsen
 * cleaned up for saber-c
 *
 * Revision 1.3  1992/05/26  20:19:00  jannsen
 * small change to support a buggy C compiler
 *
 * Revision 1.2  1992/04/06  12:40:14  seidl
 * merge in sandia changes
 *
 * Revision 1.1.1.1  1992/03/17  16:27:45  seidl
 * DOE-NIH Quantum Chemistry Library 0.0
 *
 * Revision 1.1  1992/03/17  16:27:43  seidl
 * Initial revision
 *
 * Revision 1.2  1992/02/26  14:19:09  seidl
 * free uniq_atoms when done
 *
 * Revision 1.1  1992/01/27  12:53:04  seidl
 * Initial revision
 *
 * Revision 1.4  1992/01/21  11:33:32  seidl
 * use libintv2
 *
 * Revision 1.3  1992/01/02  12:37:17  seidl
 * zero out trans array
 *
 * Revision 1.2  1991/12/20  15:44:37  seidl
 * took out function make_rp_d and put it in the file mkrpd.c
 *
 * Revision 1.1  1991/12/17  21:55:18  seidl
 * Initial revision
 *
 * Revision 1.1  1991/12/03  16:50:10  etseidl
 * Initial revision
 *
 * Revision 1.2  1991/12/02  19:53:33  seidl
 * *** empty log message ***
 *
 * Revision 1.1  1991/11/22  18:28:13  seidl
 * Initial revision
 * */

static char rcsid[] = "$Id$";

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <tmpl.h>
#include <math/array/math_lib.h>
#include <util/keyval/ipv2c.h>
#include <chemistry/qc/intv2/int_libv2.h>

#ifdef INT_SH_AM
#undef INT_SH_AM
#endif

#define INT_SH_AM(c,s) ((c)->center[(c)->center_num[s]].basis.shell[(c)->shell_num[s]].type[0].am)

#ifndef IOFF
#define IOFF(a,b) (a>b)?a*(a+1)/2+b:b*(b+1)/2+a
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include "symm.h"
#include "symm_mac.h"
#include "symmallc.h"
#include "symmfree.h"
#include "symmzero.h"

#include "symerr.gbl"

#include "syminit.gbl"
#include "syminit.lcl"

/* This reads centers info from the input and initializes the centers and
 * sym_info structures. */
GLOBAL_FUNCTION int
sym_init_centers(_centers,_sym_info,point_group,_outfile)
centers_t *_centers;
sym_struct_t *_sym_info;
char *point_group;
FILE *_outfile;
{
  centers_t old_centers;
  int errcod;
  int nat;
  int i,j;

/* let's take a stab at generating the non-unique atoms, etc */

  errcod = ip_count("atoms",&nat,0);
  if(errcod != 0) {
    serror(_outfile,__FILE__,"there don\'t appear to be any atoms\n",
     (__LINE__)-3);
    return(-1);
    }

  i=0;
  errcod = ip_count("geometry",&i,0);
  if (i!=nat) {
    serror(_outfile,__FILE__,"size of atoms != size of geometry\n",
     (__LINE__)-3);
    return(-1);
    }

  allocbn_centers(&old_centers,"n",nat);
  if(errcod != 0) {
    serror(_outfile,__FILE__,"could not allocate memory for centers",
     (__LINE__)-3);
    return(-1);
    }

 /* loop over nat, and initialize the centers structure */
  for(i=0; i < nat ; i++) {
    center_t *_center = &old_centers.center[i];
    char* at_lab;
    char* bset;

    init_center(_center);

    errcod = ip_string("atoms",&at_lab,1,i);
    if (errcod!=IPE_OK) return -1;
    for(j=0; j < strlen(at_lab); j++) at_lab[j] = (char) toupper(at_lab[j]);
    _center->atom = at_lab;

    errcod = ip_data("charges","%lf",&_center->charge,1,i);
    if(errcod != IPE_OK) _center->charge = atom_to_an(_center->atom);
    
    _center->r = (double *) malloc(sizeof(double)*3);
    if (!_center->r) return -1;
    
    errcod = ip_count("basis",&j,0);
    if(errcod == IPE_NOT_AN_ARRAY) {
      errcod = ip_string("basis",&bset,0);
      if (errcod!=IPE_OK) return -1;
      }
    else {
      errcod = ip_string("basis",&bset,1,i);
      if (errcod!=IPE_OK) return -1;
      }
    
    errcod = int_read_basis(sym_to_atom(_center->atom),bset,&(_center->basis));
    if (errcod!=IPE_OK) return -1;


    for(j=0; j < 3 ; j++) {
      errcod = ip_data("geometry","%lf",&_center->r[j],2,i,j);
        if(errcod != 0) {
          serror(_outfile,__FILE__,"problem reading the geometry",
           (__LINE__)-3);
          return(-1);
          }
        }
    }

  sym_init_given_centers(&old_centers,_centers,_sym_info,point_group,_outfile);

  free_centers(&old_centers);
}

/* This inits centers based on the info found in _old_centers.
 * It does not read any input. */
GLOBAL_FUNCTION int
sym_init_given_centers(_old_centers,_centers,_sym_info,point_group,_outfile)
centers_t *_old_centers;
centers_t *_centers;
sym_struct_t *_sym_info;
char *point_group;
FILE *_outfile;
{
  int i,j,ij;
  int g, ng, nirr,gc;
  int gi,gj,gij,nij;
  int leave;
  int nsh, ns;
  int atom, nb;
  int errcod;
  int n;
  int nat,totnat;
  int f_exist=0;
  int g_exist=0;
  int h_exist=0;
  double theta;
  double_array3_t trans;
  int_vector_t n_eq_atoms;
  int_vector_t f_n_eq_atoms;
  int_vector_t first;
  int_vector_t last;
  double_matrix_t uniq_atoms;
  double tr[3];
  double tol = 1.0e-10;
  enum pgroups pg;
  char errmsg[81];

/* the first order of business is to find out what point group we're
 * dealing with */

  if(!point_group || point_group[0] == 0) {
    fprintf(_outfile,"sym_init_centers: no symmetry specified\n");
    return(-1);
    }

/* now parse the point group symbol, this will give us the order of the
 * point group(g), the type of point group (pg), the order of the principle
 * rotation axis (n), and the number of irreps (nirr)
 */
  errcod = parse_symbol(point_group,&g,&pg,&n,&nirr,_outfile);
  if(errcod != 0) return(errcod);

  theta = (double) 2.0*M_PI/n;

/* set up atomic transformation matrices */

  errcod = allocbn_double_array3(&trans,"n1 n2 n3",g,3,3);
  if(errcod != 0) {
   serror(_outfile,__FILE__,"could not allocate memory for trans",(__LINE__)-2);
   return(-1);
   }
  zero_double_array3(&trans);

  errcod = make_trans(&trans,n,theta,pg,_outfile);
  if(errcod != 0) {
    sprintf(errmsg,"cannot handle %s yet",point_group);
    serror(_outfile,__FILE__,errmsg,(__LINE__)-3);
    return(-1);
    }

/* let's take a stab at generating the non-unique atoms, etc */

  nat = _old_centers->n;
  if(!nat) {
    serror(_outfile,__FILE__,"there don\'t appear to be any atoms\n",
     (__LINE__)-3);
    return(-1);
    }

  errcod = allocbn_double_matrix(&uniq_atoms,"n1 n2",nat,3);
  if(errcod != 0) {
    serror(_outfile,__FILE__,"could not allocate memory for uniq_atoms",
     (__LINE__)-3);
    return(-1);
    }
    
  errcod = allocbn_int_vector(&n_eq_atoms,"n",nat);
  if(errcod != 0) {
    serror(_outfile,__FILE__,"could not allocate memory for n_eq_atoms",
     (__LINE__)-3);
    return(-1);
    }

 /* loop over nat, and apply each symop to r, test to see if new */
  totnat=0;
  for(i=0; i < nat ; i++) {
    for(j=0; j < 3 ; j++) {
        uniq_atoms.d[i][j] = _old_centers->center[i].r[j];
        }
    n_eq_atoms.i[i] = test_r(&trans,uniq_atoms.d[i]);
    totnat += n_eq_atoms.i[i];
    }

 /* now that we know how many atoms we have, we can generate their coordinates
  * and get the full centers struct.
  */

  errcod = allocbn_centers(_centers,"n",totnat);
  if(errcod != 0) {
    serror(_outfile,__FILE__,"could not allocate memory for _centers",
      (__LINE__)-3);
    return(-1);
    }

  atom = 0;
  for(i=0; i < nat ; i++) {
    for(ng=0; ng < n_eq_atoms.i[i]; ng++) {
      init_center(&_centers->center[atom]);
      errcod = assign_center(&_centers->center[atom],&_old_centers->center[i]);
      if(errcod != 0) {
        sprintf(errmsg,"trouble making _centers->center[%d]",atom);
        serror(_outfile,__FILE__,errmsg,(__LINE__)-5);
        return(-1);
        }
      for (j=0; j<3; j++) {
        _centers->center[atom].r[j] = uniq_atoms.d[i][j];
        }
      atom++;
      }
    }

  atom = 0;
  for(i=0; i < nat ; i++) {
    for(j=0; j < 3; j++) tr[j]=uniq_atoms.d[i][j];
    if(fabs(tr[0]) < tol && fabs(tr[1]) < tol && fabs(tr[2]) < tol) {
      for(j=0; j < 3; j++) _centers->center[atom].r[j]=tr[j];
      atom++;
      }
    else {
      for(j=0; j < 3; j++) _centers->center[atom].r[j]=tr[j];
      atom++;
      for(ng=1; ng < g; ng++) {
        new_r(uniq_atoms.d[i],tr,trans.d[ng]);
        if(comp_r(_centers,tr) < 0 ) {
          for(j=0; j < 3; j++) _centers->center[atom].r[j]=tr[j];
          atom++;
          }
        }
      }
    }

/* now initialize the full centers struct */

  int_normalize_centers(_centers);
  int_initialize_1e(0,0,_centers,_centers);
  int_initialize_offsets1(_centers,_centers);
  nsh = _centers->nshell;

/* let's see if there are f or g functions. we'll puke if h or higher are
 * used. */

  for(i=0; i < nsh ; i++) {
    j = INT_SH_AM(_centers,i);
    if(j==3) f_exist=1;
    else if(j==4) g_exist=1;
    else if(j >= 5) h_exist=1;
    }  
  if(h_exist) {
    sprintf(errmsg,"shell %d am %d: this program cannot handle am > 4",i,j);
    serror(_outfile,__FILE__,errmsg,(__LINE__)-2);
    return(-1);
    }

 /* and allocate memory for sym_info */

  if(g==1) {
    errcod = allocbn_sym_struct(_sym_info,
       "point_group pg",point_group,pg);
    if(errcod != 0) {
      serror(_outfile,__FILE__,"could not allocate memory for _sym_info",
      (__LINE__)-5);
      return(-1);
      }
    }
  else if(g_exist) {
    errcod = allocbn_sym_struct(_sym_info,
     "natom g nshell nshtri point_group pg nirrep n1 n2 n2p n3 n3p n4 n4p",
      totnat,g,nsh,nsh*(nsh+1)/2,point_group,pg,nirr,3,6,5,10,7,15,9);
    if(errcod != 0) {
      serror(_outfile,__FILE__,"could not allocate memory for _sym_info",
      (__LINE__)-5);
      return(-1);
      }
    }
  else if(f_exist) {
    errcod = allocbn_sym_struct(_sym_info,
     "natom g nshell nshtri point_group pg nirrep n1 n2 n2p n3 n3p",
     totnat,g, nsh,nsh*(nsh+1)/2,point_group,pg,nirr,3,6,5,10,7);
    if(errcod != 0) {
      serror(_outfile,__FILE__,"could not allocate memory for _sym_info",
      (__LINE__)-5);
      return(-1);
      }
    }
  else {
    errcod = allocbn_sym_struct(_sym_info,
     "natom g nshell nshtri point_group pg nirrep n1 n2 n2p",
     totnat,g, nsh,nsh*(nsh+1)/2,point_group,pg,nirr,3,6,5);
    if(errcod != 0) {
      serror(_outfile,__FILE__,"could not allocate memory for _sym_info",
      (__LINE__)-5);
      return(-1);
      }
    }

  zero_sym_struct(_sym_info);

  if(g>1) {
    errcod = allocbn_char_tab(&_sym_info->ct,"nirrep",nirr);
    if(errcod != 0) {
      serror(_outfile,__FILE__,"could not allocate memory for char_tab",
       (__LINE__)-3);
      return(-1);
      }

    errcod = allocbn_int_vector(&f_n_eq_atoms,"n",totnat);
    if(errcod != 0) {
      serror(_outfile,__FILE__,"could not allocate memory for f_n_eq_atoms",
       (__LINE__)-3);
      return(-1);
      }
    errcod = allocbn_int_vector(&first,"n",totnat);
    if(errcod != 0) {
      serror(_outfile,__FILE__,"could not allocate memory for first",
       (__LINE__)-3);
      return(-1);
      }
    errcod = allocbn_int_vector(&last,"n",totnat);
    if(errcod != 0) {
      serror(_outfile,__FILE__,"could not allocate memory for last",
       (__LINE__)-3);
      return(-1);
      }

    for(i=ij=0; i < nat ; i++)
      for(j=0; j < n_eq_atoms.i[i] ; j++,ij++)
        f_n_eq_atoms.i[ij] = n_eq_atoms.i[i];

    ns=nb=0;
    for(i=0; i < totnat ; i++) {
      first.i[i]=ns;
      last.i[i]=ns+_centers->center[i].basis.n;
      ns+=_centers->center[i].basis.n;

      for(ng=0; ng < g ; ng++) {
        new_r(_centers->center[i].r,tr,trans.d[ng]);
        _sym_info->atom_map[i][ng] = comp_r(_centers,tr);
        }

      _sym_info->first[i] = nb;
      for(j=0; j < _centers->center[i].basis.n ; j++) {
        nb += _centers->center[i].basis.shell[j].nfunc;
        }
      }

    nb=0;
    for(i=0; i < totnat ; i++) {
      _sym_info->firstp[i] = nb;
      for(j=0; j < _centers->center[i].basis.n ; j++) {
        for(gc=0; gc < _centers->center[i].basis.shell[j].ncon ; gc++) {
          ns = _centers->center[i].basis.shell[j].type[gc].am;
          nb += 2*ns+1;
          }
        }
      }

    for(i=0; i < totnat ; i++) {
      for(j=0,ns=first.i[i]; ns < last.i[i]; j++,ns++) {
        for(ng=0; ng < g ; ng++) {
          atom=_sym_info->atom_map[i][ng];
          _sym_info->shell_map[ns][ng] = first.i[atom]+j;
          }
        }
      }

   /* now let's set up lambda from Dupuis and King's paper */
  
    for(i=0; i < nsh ; i++) {
      for(ng=leave=0; ng<g ; ng++) if(_sym_info->shell_map[i][ng] > i) leave=1;
      if(leave) continue;
      _sym_info->p1[i]=1;
      for(j=0; j <=i ; j++) {
        ij=IOFF(i,j);
        nij=0;
        for(ng=leave=0; ng < g ; ng++) {
          gi = _sym_info->shell_map[i][ng];
          gj = _sym_info->shell_map[j][ng];
          gij = IOFF(gi,gj);
          if(gij > ij) leave=1;
          if(gij == ij) nij++;
          }
        if (leave) continue;
        _sym_info->lamij[ij] = (char) (g/nij);
        }
      }

   /* finally, let's form Rp and Rd in sym_info */

    errcod = make_rp_d(_sym_info,n,theta,pg,f_exist,g_exist,_outfile);
    if(errcod != 0) {
      sprintf(errmsg,"cannot handle %s yet",point_group);
      serror(_outfile,__FILE__,errmsg,(__LINE__)-3);
      return(-1);
      }

    free_int_vector(&f_n_eq_atoms);
    free_int_vector(&first);
    free_int_vector(&last);
    }

  int_done_1e();
  int_done_offsets1(_centers,_centers);

  free_double_array3(&trans);
  free_double_matrix(&uniq_atoms);
  free_int_vector(&n_eq_atoms);

  return(0);
  }

LOCAL_FUNCTION int
comp_r(c,tr)
centers_t *c;
double *tr;
{
  int i;
  double tol=1.0e-5;

  for(i=0; i < c->n ; i++) {
    if(fabs(c->center[i].r[0]-tr[0]) < tol &&
       fabs(c->center[i].r[1]-tr[1]) < tol &&
       fabs(c->center[i].r[2]-tr[2]) < tol) return(i);
    }
  return(-1);
  }

LOCAL_FUNCTION VOID
new_r(r,tr,tm)
double *r;
double *tr;
double **tm;
{
  int i,j;

  for(i=0; i < 3 ; i++) {
    tr[i]=0.0;
    for(j=0; j < 3 ; j++) tr[i] += tm[i][j]*r[j];
    }
  }

LOCAL_FUNCTION int
test_r(tm,r)
double_array3_t *tm;
double *r;
{
  int i,j,g;
  int nmap;
  double *tm_g_i;
  double tr[3];
  double tol = 1.0e-5;

  tr[0] = tr[1] = tr[2] = 0.0;

  nmap=0;
  for(g=0; g < tm->n1 ; g++) {
    for(i=0; i < 3 ; i++) {
      tm_g_i = tm->d[g][i];
      tr[i]=0.0;
      for(j=0; j < 3 ; j++) tr[i] += tm_g_i[j]*r[j];
      }
    if(fabs(r[0]-tr[0]) < tol &&
       fabs(r[1]-tr[1]) < tol &&
       fabs(r[2]-tr[2]) < tol) nmap++;
    }
  return(tm->n1/nmap);
  }


LOCAL_FUNCTION int
make_trans(trans, n, theta, pg, _outfile)
double_array3_t *trans;
int n;
double theta;
enum pgroups pg;
FILE *_outfile;
{
  int i,j;

  switch(pg) {
  case _PG_C1:
    trans->d[0][0][0] = trans->d[0][1][1] = trans->d[0][2][2] = 1.0;
    break;
  case _PG_CI:
    trans->d[0][0][0] = trans->d[0][1][1] = trans->d[0][2][2] = 1.0;
    trans->d[1][0][0] = trans->d[1][1][1] = trans->d[1][2][2] = -1.0;
    break;
  case _PG_CS:
    trans->d[0][0][0] = trans->d[0][1][1] = trans->d[0][2][2] = 1.0;
    trans->d[1][0][0] = trans->d[1][1][1] = 1.0;
    trans->d[1][2][2] = -1.0;
    break;
  case _PG_CN:
    for(i=0; i < n ; i++) {
      trans->d[i][0][0] = trans->d[i][1][1] = cos((double)theta*i);
      trans->d[i][0][1] = sin((double)theta*i);
      trans->d[i][1][0] = -sin((double)theta*i);
      trans->d[i][2][2] = 1.0;
      }
    break;
  case _PG_CNV:
    for(i=0; i < n ; i++) {
      trans->d[i][0][0] = trans->d[i][1][1] = cos((double)theta*i);
      trans->d[i][0][1] = sin((double)theta*i);
      trans->d[i][1][0] = -sin((double)theta*i);
      trans->d[i][2][2] = 1.0;
      }

    for(i=0,j=n; i < n ; i++,j++) {
      trans->d[j][0][0] = cos((double)theta*i);
      trans->d[j][1][1] = -cos((double)theta*i);
      trans->d[j][1][0] = trans->d[j][0][1] = sin((double)theta*i);
      trans->d[j][2][2] = 1.0;
      }
    break;
  case _PG_CNH:
    for(i=0; i < n ; i++) {
      trans->d[i][0][0] = trans->d[i][1][1] =
      trans->d[i+n][0][0] = trans->d[i+n][1][1] = cos((double)theta*i);
      trans->d[i][0][1] = trans->d[i+n][0][1] = sin((double)theta*i);
      trans->d[i][1][0] = trans->d[i+n][1][0] = -sin((double)theta*i);
      trans->d[i][2][2] = 1.0;
      trans->d[i+n][2][2] = -1.0;
      }
    break;
  case _PG_SN:
    for(i=0; i < n ; i++) {
      trans->d[i][0][0] = trans->d[i][1][1] = cos((double)theta*i);
      trans->d[i][0][1] = sin((double)theta*i);
      trans->d[i][1][0] = -sin((double)theta*i);
      trans->d[i][2][2] = pow(-1.0,(double) i);
      }
    break;
  case _PG_DN:
    for(i=0; i < n ; i++) {
      trans->d[i][0][0] = trans->d[i][1][1] = cos((double) theta*i);
      trans->d[i][0][1] = sin((double) theta*i);
      trans->d[i][1][0] = -sin((double) theta*i);
      trans->d[i][2][2] = 1.0;
      }

    for(i=0,j=n; i < n ; i++,j++) {
      trans->d[j][0][0] = cos((double)theta*i);
      trans->d[j][1][1] = -cos((double)theta*i);
      trans->d[j][1][0] = trans->d[j][0][1] = -sin((double)theta*i);
      trans->d[j][2][2] = -1.0;
      }
    break;
  case _PG_DND:
    for(i=0; i < 2*n ; i++) {
      trans->d[i][0][0] = trans->d[i][1][1] = cos((double)theta*i*0.5);
      trans->d[i][0][1] = sin((double)theta*i*0.5);
      trans->d[i][1][0] = -sin((double)theta*i*0.5);
      trans->d[i][2][2] = pow(-1.0,(double) i);
      }

    for(i=0,j=2*n; i < n ; i++,j++) {
      trans->d[j][0][0] = cos((double)theta*i);
      trans->d[j][1][1] = -cos((double)theta*i);
      trans->d[j][1][0] = trans->d[j][0][1] = -sin((double)theta*i);
      trans->d[j][2][2] = -1.0;
      }

    for(i=0,j=3*n; i < n ; i++,j++) {
      trans->d[j][0][0] = cos((double)theta*i+theta*0.5);
      trans->d[j][1][1] = -cos((double)theta*i+theta*0.5);
      trans->d[j][1][0] = trans->d[j][0][1] = -sin((double)theta*i+theta*0.5);
      trans->d[j][2][2] = 1.0;
      }
    break;
  case _PG_DNH:
    for(i=0; i < n ; i++) {
      trans->d[i][0][0] = trans->d[i][1][1] =
      trans->d[i+n][0][0] = trans->d[i+n][1][1] = cos((double)theta*i);
      trans->d[i][0][1] = trans->d[i+n][0][1] = sin((double)theta*i);
      trans->d[i][1][0] = trans->d[i+n][1][0] = -sin((double)theta*i);
      trans->d[i][2][2] = 1.0;
      trans->d[i+n][2][2] = -1.0;
      }

    for(i=0,j=2*n; i < n ; i++,j++) {
      trans->d[j][0][0] = cos((double)theta*i);
      trans->d[j][1][1] = -cos((double)theta*i);
      trans->d[j][1][0] = trans->d[j][0][1] = -sin((double)theta*i);
      trans->d[j][2][2] = -1.0;
      }

    for(i=0,j=3*n; i < n ; i++,j++) {
      trans->d[j][0][0] = cos((double)theta*i);
      trans->d[j][1][1] = -cos((double)theta*i);
      trans->d[j][1][0] = trans->d[j][0][1] = sin((double)theta*i);
      trans->d[j][2][2] = 1.0;
      }
    break;
  default:
    return(-1);
    }

  return(0);
  }

LOCAL_FUNCTION int
parse_symbol(point_group,g,pg,n,nirr,_outfile)
char *point_group;
int *g;
enum pgroups *pg;
int *n;
int *nirr;
FILE *_outfile;
{
  char errmsg[81];
  sprintf(errmsg,"unknown Schoenflies symbol: %s",point_group);

  *n = 1;
  *g = 1;
  *pg = _PG_C1;

  if(!strcmp(point_group,"ci")) {
    *g = 2;
    *pg = _PG_CI;
    *nirr = 2;
    }
  else if(!strcmp(point_group,"c1")) {
    *g = 1;
    *pg = _PG_C1;
    *nirr = 1;
    }
  else if(!strcmp(point_group,"cs")) {
    *g = 2;
    *pg = _PG_CS;
    *nirr = 2;
    }
  else if(point_group[0]=='c') {
    int nab,ne;

    if(point_group[1] == 0) {
      serror(_outfile,__FILE__,errmsg,(__LINE__)-2);
      return(-1);
      }
    else {
      *n = atoi(&point_group[1]);
      ne = (*n%2) ? *n/2 : *n/2-1;
      nab = (*n%2) ? 1 : 2;
      }
    if(point_group[2] != 0) {
      if(point_group[2] == 'v') {
        *g  = 2* *n;
        *pg = _PG_CNV;
        *nirr = 2*nab + ne;
        }
      else if(point_group[2] == 'h') {
        *g  = 2* *n;
        *pg = _PG_CNH;
        *nirr = 2*(nab+ne);
        }
      else {
        serror(_outfile,__FILE__,errmsg,(__LINE__)-2);
        return(-1);
        }
      }
    else {
      *g = *n;
      *pg = _PG_CN;
      *nirr = nab+ne;
      }
    }
  else if(point_group[0]=='d') {
    int nab,ne;

    if(point_group[1] == 0) {
      serror(_outfile,__FILE__,errmsg,(__LINE__)-2);
      return(-1);
      }
    *n = atoi(&point_group[1]);
    ne = (*n%2) ? *n/2 : *n/2-1;
    nab = (*n%2) ? 1 : 2;

    if(point_group[2] != 0) {
      if(point_group[2] == 'd') {
        *g = 4* *n;
        *pg = _PG_DND;
        *nirr = *n+3;
        }
      else if(point_group[2] == 'h') {
        *g = 4* *n;
        *pg = _PG_DNH;
        *nirr = 4*nab + 2*ne;
        }
      else {
        serror(_outfile,__FILE__,errmsg,(__LINE__)-2);
        return(-1);
        }
      }
    else {
      *g = 2* *n;
      *pg = _PG_DN;
      *nirr = 2*nab + ne;
      }
    }
  else if(point_group[0]=='s') {
    if(point_group[1] == 0) {
      serror(_outfile,__FILE__,errmsg,(__LINE__)-2);
      return(-1);
      }
    *n = atoi(&point_group[1]);
 /* only S2n groups make sense */
    if(*n%2) {
      serror(_outfile,__FILE__,errmsg,(__LINE__)-2);
      return(-1);
      }
    *g = *n;
    *pg=_PG_SN;
    *nirr = *n/2+1;
    }
  else if(point_group[0]=='t') {
    if(point_group[1] != 0) {
      if(point_group[1] == 'd') {
        *g = 24;
        *pg = _PG_TD;
        *nirr = 5;
        }
      else if(point_group[1] == 'h') {
        *g = 24;
        *pg = _PG_TH;
        *nirr = 6;
        }
      else {
        serror(_outfile,__FILE__,errmsg,(__LINE__)-2);
        return(-1);
        }
      }
    else {
      *g = 12;
      *pg = _PG_T;
      *nirr = 3;
      }
    }
  else if(point_group[0]=='o') {
    if(point_group[1] != 0) {
      if(point_group[1] == 'h') {
        *pg = _PG_OH;
        *g = 48;
        *nirr = 10;
        }
      else {
        serror(_outfile,__FILE__,errmsg,(__LINE__)-2);
        return(-1);
        }
      }
    else {
      *g = 24;
      *pg = _PG_O;
      *nirr = 5;
      }
    }
  else if(point_group[0]=='i') {
    if(point_group[1] != 0) {
      if(point_group[1] == 'h') {
        *g = 120;
        *pg = _PG_IH;
        *nirr = 10;
        }
      else {
        serror(_outfile,__FILE__,errmsg,(__LINE__)-2);
        return(-1);
        }
      }
    else {
      *g = 60;
      *pg = _PG_I;
      *nirr = 5;
      }
    }
  else {
    serror(_outfile,__FILE__,errmsg,(__LINE__)-2);
    return(-1);
    }
  return(0);
  }

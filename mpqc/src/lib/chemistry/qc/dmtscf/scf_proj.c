
/* $Log$
 * Revision 1.2  1994/06/08 01:15:13  cljanss
 * Many changes.  These include: newmat7 and nihmatrix -> scmat
 * and mpqcic -> MPSCF and updated optimize stuff.
 *
 * Revision 1.1.1.1  1993/12/29  12:53:17  etseidl
 * SC source tree 0.1
 *
 * Revision 1.6  1993/04/27  23:54:56  jannsen
 * rs/6000 port
 *
 * Revision 1.5  1992/06/29  17:50:01  seidl
 * use dmt matrices now
 *
 * Revision 1.4  1992/06/17  21:53:15  jannsen
 * cleaned up for saber-c and changed to ngl loops
 *
 * Revision 1.3  1992/05/26  20:18:11  jannsen
 * use mtype_get to get message types for global operations
 * check results of memory allocations
 *
 * Revision 1.2  1992/05/04  11:08:38  seidl
 * remove sync0 call
 *
 * Revision 1.1  1992/04/22  16:02:33  seidl
 * add to math/dmt/libdmtscf...i think it works
 * */

/* scf_project_vector tries to take an old vector in one basis, and project
 * it onto the present basis.
 */

static char rcsid[]="$Id$";

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <tmpl.h>
#include <fcntl.h>
#ifdef RS6000
#define O_RDONLY 0
#endif

#include <math/dmt/libdmt.h>
#include <util/bio/libbio.h>
#include <util/misc/libmisc.h>
#include <math/array/math_lib.h>
#include <chemistry/qc/intv2/int_libv2.h>
#include <util/keyval/ipv2c.h>
#include <chemistry/qc/dmtsym/sym_dmt.h>

#include "scf.h"
#include "scf_orth.gbl"

#include "scf_proj.gbl"
#include "scf_proj.lcl"

#define MAX0(a,b) ((a)>(b))?(a):(b)
#define MIN0(a,b) ((a)<(b))?(a):(b)

/* reads in an old basis set and calls scf_project_vector_ to do
 * a projection */
GLOBAL_FUNCTION int
scf_project_vector(centers, irreps, scf_info, Scf_Vec, oldvecfile, outfile)
centers_t *centers;
scf_irreps_t *irreps;
scf_struct_t *scf_info;
dmt_matrix Scf_Vec;
char *oldvecfile;
FILE *outfile;
{
  int i,j;
  int m,n;
  int errcod,count;
  int me=mynode0();
  int basis_array=1;
  char *bset;
  centers_t old_cen;
  center_t *center;

/* construct a centers struct using the old basis set */

  if(me==0) {
    ip_cwk_push();
    ip_cwk_add(":default");
    ip_cwk_add(":project");

    count=0;
    errcod = ip_count("basisfiles",&count,0);
    if(count) ip_append_from_input("basis",outfile);

    errcod = ip_count("oldbasis",&i,0);
    if(errcod == IPE_NOT_AN_ARRAY) {
      basis_array=0;
      errcod = ip_string("oldbasis",&bset,0);
      if (errcod!=IPE_OK) {
        fprintf(stderr,"scf_proj: cannot find oldbasis\n");
        errcod=-1;
        bcast0(&errcod,sizeof(int),mtype_get(),0);
        goto the_place_where_errors_go;
        }
      }

    errcod=allocbn_centers(&old_cen,"n",centers->n);
    if(errcod!=0) {
      fprintf(stderr,"scf_proj: could not allocate old_cen\n");
      errcod=-1;
      bcast0(&errcod,sizeof(int),mtype_get(),0);
      goto the_place_where_errors_go;
      }

    for(i=0; i < centers->n; i++) {
      center = &old_cen.center[i];

      if(basis_array) {
        errcod = ip_string("oldbasis",&bset,1,i);
        if (errcod!=IPE_OK) {
          errcod=-1;
          bcast0(&errcod,sizeof(int),mtype_get(),0);
          goto the_place_where_errors_go;
          }
        }

      errcod=allocbn_center(center,"atom charge",
                            centers->center[i].atom,
                            centers->center[i].charge);
      if (errcod!=0) {
        fprintf(stderr,"scf_proj: could not alloc center %d\n",i);
        errcod=-1;
        bcast0(&errcod,sizeof(int),mtype_get(),0);
        goto the_place_where_errors_go;
        }

      center->r[0]=centers->center[i].r[0];
      center->r[1]=centers->center[i].r[1];
      center->r[2]=centers->center[i].r[2];

      errcod = read_basis(sym_to_atom(center->atom),bset,&(center->basis));
      if (errcod!=0) {
        fprintf(stderr,"scf_proj: could not read basis %d\n",i);
        errcod=-1;
        bcast0(&errcod,sizeof(int),mtype_get(),0);
        goto the_place_where_errors_go;
        }
      }

    int_normalize_centers(&old_cen);

    ip_cwk_pop();

    errcod=0;
    bcast0(&errcod,sizeof(int),mtype_get(),0);

the_place_where_errors_go:
    ;
    }
  else {
    bcast0(&errcod,sizeof(int),mtype_get(),0);
    }

  if(errcod!=0) return(-1);

  return scf_project_vector_(&old_cen, centers, irreps, scf_info, Scf_Vec,
                             oldvecfile, outfile);
}

/* projects the basis set from old_cen to centers, the basis functions
 * in old_cen should be normalized and only on node 0. */
GLOBAL_FUNCTION int
scf_project_vector_(old_cen, centers, irreps, scf_info,
                   Scf_Vec, oldvecfile, outfile)
centers_t *old_cen;
centers_t *centers;
scf_irreps_t *irreps;
scf_struct_t *scf_info;
dmt_matrix Scf_Vec;
char *oldvecfile;
FILE *outfile;
{
  int i,j;
  int m,n;
  int errcod;
  int me=mynode0();
  int nlocal=dmt_nlocal(Scf_Vec);
  int icol,nocc;
  int fd;
  double *col;
  loop_t *loop;

  dmt_matrix S=dmt_old("libscfv3 overlap matrix");
  dmt_matrix Shalf = dmt_old("libscfv3 scf_core_guess scr4");

  dmt_matrix SCR1=dmt_create("math/dmt/libdmtscf proj scr1",scf_info->nbfao,COLUMNS);
  dmt_matrix SCR2=dmt_create("math/dmt/libdmtscf proj scr2",scf_info->nbfao,COLUMNS);
  dmt_matrix SCR3=dmt_create("math/dmt/libdmtscf proj scr3",scf_info->nbfao,COLUMNS);
  dmt_matrix Cprime=dmt_create("math/dmt/libdmtscf proj scr4",scf_info->nbfao,COLUMNS);

  dmt_matrix Cnew,X,D,Sinv,C2,S2;
  center_t *center;

  nocc=irreps->ir[0].nclosed+irreps->ir[0].nopen;

/* form Sinv = (S**-1/2)**2 */

  Sinv=SCR1;
  dmt_copy(Shalf,SCR2);
  dmt_mult(SCR2,Shalf,Sinv);

  dmt_free(Shalf);

 /* broadcast old centers struct to nodes */

  bcast0_centers(old_cen,0,0);

  int_initialize_1e(0,0,centers,old_cen);
  int_initialize_offsets1(centers,old_cen);

  m=old_cen->nfunc;
  n=centers->nfunc;

 /* read old vector */

  C2=SCR2;
  dmt_fill(C2,0.0);

  if(me==0) { 
    if((fd = open(oldvecfile,O_RDONLY))<0) {
      fprintf(stderr,"scf_proj: could not open file %s\n",oldvecfile);
      errcod=-1;
      bcast0(&errcod,sizeof(int),mtype_get(),0);
      } 
    else {
      errcod=0;
      bcast0(&errcod,sizeof(int),mtype_get(),0);
      }
    }
  else {
    errcod=0;
    bcast0(&errcod,sizeof(int),mtype_get(),0);
    }
  if(errcod!=0) return(-1);
 
  loop = dmt_ngl_create("%m",C2);
  while(dmt_ngl_next(loop)) {
    int iind,isize,jsize;

    dmt_ngl_create_inner(loop,0);
    while(dmt_ngl_next_inner_m(loop,&iind,&isize,&i,&jsize,&col)) {
      if(i<nocc && me==0) {
        read(fd,col,sizeof(double)*m);
        }
      }
    }
  dmt_ngl_kill(loop);

  if(me==0) {
    close(fd);
    fprintf(outfile,"  read old scf vector from %s\n",oldvecfile);
    }

 /* form overlap between old basis and new basis */

  S2=SCR3;
  form_s2(centers,old_cen,S2);

 /* form projected vector by multiplying the old vector by the overlap */

  dmt_mult(S2,C2,Cprime);

 /* free some memory */
  int_done_offsets1(centers,old_cen);
  int_done_1e();

  free_centers(old_cen);

 /* form D = Sinv x Cprime */

  D=SCR2;

  dmt_mult(Sinv,Cprime,D);

 /* form X = Cprime~ * D */

  X=SCR3;

  dmt_mult(Cprime,D,X);
  dmt_free(Cprime);

 /* form X**-1/2 */

  errcod = form_x_half(X,nocc);
  if(errcod!=0) {
    fprintf(stderr,"scf_proj: trouble in form_x_half\n");
    return(-1);
    }

 /* finally, form Cnew = D*X**-1/2 */

  Cnew=SCR1;

  dmt_transpose(D);
  dmt_mult(D,X,Cnew);

 /* and put Cnew into Scf_Vec */

  for(i=0; i < nlocal ; i++) {
    double *ncol;

    dmt_get_col(Scf_Vec,i,&icol,&col);
    dmt_get_col(Cnew,i,&icol,&ncol);
    if(icol<nocc) for(j=0; j < n; j++) col[j]=ncol[j];
    }

 /* now orthogonalize Scf_Vec */

  errcod = scf_schmidt(scf_info,irreps,Scf_Vec,S,1,outfile);
  if(errcod !=0 ) {
    fprintf(outfile,"trouble orthogonalizing new vector\n");
    return(-1);
    }

  dmt_free(SCR1);
  dmt_free(SCR2);
  dmt_free(SCR3);

  return(0);
  }

#define COUNT 1
#define READ 2

static char *prb = "encountered a problem while reading the basis set:\n ";

LOCAL_FUNCTION int
read_basis(atom,basisname,_basis)
char *atom;
char *basisname;
basis_t *_basis;
{
  int errcod;

  /* Copy the basisname to the basis structure. */
  _basis->name = (char *)malloc(strlen(basisname)+1);
  check_alloc(_basis->name,"_basis->name");
  strcpy(_basis->name,basisname);

  /* If the basis name is "nothing" then return with an empty basis set. */
  if (!strcmp(basisname,"nothing")) {
    _basis->n = 0;
    _basis->shell = NULL;
    return IPE_OK;
    }

  if ((errcod = read_basis_(COUNT,atom,basisname,_basis))!=IPE_OK) return errcod;

  _basis->shell = (shell_t *) malloc(sizeof(shell_t)*_basis->n);
  check_alloc(_basis->shell,"_basis->shell");
  _basis->n = 0;
  return read_basis_(READ,atom,basisname,_basis);
  }

LOCAL_FUNCTION int
read_basis_(what,atom,basisname,basis)
int what;
char *atom;
const char *basisname;
basis_t *basis;
{
  int i;
  int errcod;
  char key[KEYWORD_LENGTH];
  const char *val;


  i = 0;

  while(1) {
    sprintf(key,":basis:%s:%s:%d",atom,basisname,i);
    errcod = ip_value_v(key,&val,0,NULL);
    if (errcod == IPE_KEY_NOT_FOUND) break;
    sprintf(key,":basis:%s:%s:%d:type",atom,basisname,i);
    errcod = ip_value_v(key,&val,0,NULL);
    if (errcod == IPE_OK) {
      ip_warn("%k has been illegally assigned a value");
      return IPE_VAL_NOT_EXPD;
      }
    else if (errcod == IPE_HAS_NO_VALUE) {
      /* We have found a shell read it in (or count). */
      if (what == READ) {
        sprintf(key,":basis:%s:%s:%d",atom,basisname,i);
        errcod = ip_read_shell_v(key,&basis->shell[basis->n],0,NULL);
        if (errcod != IPE_OK) {
          ip_warn("%scouldn't read a shell: looking for \"%k\"",prb);
          return errcod;
          }
        }
      basis->n++;
      }
    else if (errcod == IPE_KEY_NOT_FOUND) {
      /* Not an atom, see if we have a get. */
      sprintf(key,":basis:%s:%s:%d:get",atom,basisname,i);
      errcod = ip_value_v(key,&val,0,NULL);
      if (errcod != IPE_OK) {
        ip_warn("%scouldn't find a shell or a \"get\" while looking for \"%k\"",
                prb);
        return errcod;
        }
      errcod = read_basis_(what,atom,val,basis);
      if (errcod != IPE_OK) return errcod;
      }
    else {
      return errcod;
      }
    i++;
    }
  /* If i is still zero, then something was missing in the input. */
  if (!i) {
    ip_warn("%scould not find \"%k\"",prb);
    return(-1);
    }
  return IPE_OK;
  }


LOCAL_FUNCTION int
form_x_half(X,nocc)
dmt_matrix X;
int nocc;
{
  int i,j,k;
  int nlocal=dmt_nlocal(X);
  double_matrix_t eigvec,x;
  double_vector_t evals;
  loop_t *loop;

 /* allocate memory */
  if (  (allocbn_double_matrix(&eigvec,"n1 n2",nocc,nocc)!=IPE_OK)
      ||(allocbn_double_vector(&evals,"n",nocc)!=IPE_OK)) {
    fprintf(stderr,"form_x_half: alloc of eigvec and evals failed\n");
    return(-1);
    }

  if(allocbn_double_matrix(&x,"n1 n2",nocc,nocc)!=IPE_OK) {
    fprintf(stderr,"form_x_half: alloc of x failed\n");
    return(-1);
    }

 /* put X in local matrix */

  loop = dmt_ngl_create("%mr",X);
  while(dmt_ngl_next(loop)) {
    int iind,isize,jsize;
    double *col;

    dmt_ngl_create_inner(loop,0);
    while(dmt_ngl_next_inner_m(loop,&iind,&isize,&i,&jsize,&col)) {
      if(i<nocc) for(j=0; j < nocc; j++) x.d[i][j]=col[j];
      }
    }
  dmt_ngl_kill(loop);

 /* now form X**-1/2 */

  math_diag_dm(&x,&evals,&eigvec,1,1.0e-15,1);

  for(i=0; i < nocc; i++) evals.d[i] = 1.0/sqrt(evals.d[i]);

  for(i=0; i < nocc ; i++) {
    for(j=0; j < nocc; j++) {
      x.d[i][j]=0.0;
      for(k=0; k < nocc; k++) 
        x.d[i][j]+=eigvec.d[i][k]*eigvec.d[j][k]*evals.d[k];
      }
    }

 /* put x back into X */
  for(i=0; i < nlocal; i++) {
    int icol;
    double *col;
    dmt_get_col(X,i,&icol,&col);
    if(icol<nocc) for(j=0; j < nocc; j++) col[j]=x.d[icol][j];
    }

  free_double_matrix(&x);
  free_double_matrix(&eigvec);
  free_double_vector(&evals);
  return(0);
  }

LOCAL_FUNCTION VOID
form_s2(new_cen,old_cen,S2)
centers_t *new_cen;
centers_t *old_cen;
dmt_matrix S2;
{
  int i,l,k,ll,kk;
  int ni,oj,lfirst,llast;
  int ofirst,olast;
  int nfirst,nlast;
  int onfunc;
  int nfuncmax=int_find_nfuncmax(new_cen);
  int nlocal=dmt_nlocal(S2);
  int msize=dmt_size(S2);
  double *s2col,**ls2,**ls2col;
  double *block;

  dmt_fill(S2,0.0);

 /* get pointer to locally held columns */
  dmt_get_col(S2,0,&lfirst,&s2col);
  llast=lfirst+nlocal;

 /* make 2d array which points to local columns */
  ls2 = (double **) malloc(sizeof(double *)*nlocal);
  check_alloc(ls2,"scf_proj.c: form_s2: ls2");

 /* offset ls2 by lfirst, so ls2col[lfirst]==ls2[0] */
  ls2col = &ls2[-lfirst];

 /* block will hold calulated overlap integrals */
  block = (double *) malloc(sizeof(double)*nfuncmax*nfuncmax);
  check_alloc(block,"scf_proj.c: form_s2: block");

  for(i=0; i < nlocal; i++) ls2[i] = &s2col[i*msize];

 /* loop over old shells */
  for(oj=0; oj < old_cen->nshell; oj++) {
    onfunc=INT_SH_NFUNC(old_cen,oj);
    ofirst=old_cen->func_num[oj];
    olast=ofirst+onfunc;

   /* is this shell held locally? */
    if((lfirst >= ofirst && lfirst < olast) ||
       (ofirst >= lfirst && ofirst < llast) ||
       (olast > lfirst && olast <= llast) ||
       (llast > ofirst && llast <= olast)) {

     /* loop over new shells */
      for(ni=0; ni < new_cen->nshell; ni++) {
        nfirst=new_cen->func_num[ni];
        nlast=nfirst+INT_SH_NFUNC(new_cen,ni);

        int_shell_overlap(new_cen,old_cen,block,ni,oj);

        ll=MAX0(0,lfirst-ofirst);
        for(l=MAX0(ofirst,lfirst); l < (MIN0(llast,olast)); ll++,l++) {
          for(kk=0,k=nfirst; k < nlast; kk++,k++) {
            ls2col[l][k]=block[kk*onfunc+ll];
            }
          }
        }
      }
    }

 /* we'll be wanting the transpose */
  dmt_transpose(S2);

  free(ls2);
  free(block);
  }

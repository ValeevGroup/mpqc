
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fcntl.h>
#include <unistd.h>

#ifndef O_RDONLY
#define O_RDONLY 0
#endif

#include <tmpl.h>
#include <comm/picl/picl.h>
#include <comm/picl/ext/piclext.h>
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

/*****************************************************************************
 *
 * construct a centers struct using the old basis set.  
 *
 * input:
 *   centers = pointer to initialized centers struct
 *   oldcenters = pointer to uninitialized centers struct
 *
 * on return:
 *   oldcenters contains copy of centers but with old basis set
 *
 * return 0 on success, -1 on failure
 */

GLOBAL_FUNCTION int
scf_make_old_centers(centers, oldcenters)
centers_t *centers;
centers_t *oldcenters;
{
  int i;
  int count;
  int errcod;
  int basis_array;
  char *bset;

  if (mynode0() == 0) {
    ip_cwk_push();
    ip_cwk_add(":default");
    ip_cwk_add(":project");

    count=0;
    errcod = ip_count("basisfiles",&count,0);
    if (count) ip_append_from_input("basis",stderr);

   /* read the value of oldbasis.  it can be an array or a single string */
    errcod = ip_count("oldbasis",&i,0);
    if (errcod == IPE_NOT_AN_ARRAY) {
      basis_array=0;

      errcod = ip_string("oldbasis",&bset,0);
      if (errcod != IPE_OK) {
        fprintf(stderr,"scf_make_old_centers: cannot find oldbasis\n");
        goto the_place_where_errors_go;
      }
    }

   /* no allocate memory for oldcenters */
    errcod=allocbn_centers(oldcenters,"n",centers->n);
    if (errcod != 0) {
      fprintf(stderr,"scf_proj: could not allocate oldcenters\n");
      goto the_place_where_errors_go;
    }

   
  /* and then copy much of centers into oldcenters, but hold off on the
   * basis set info
   */
    for (i=0; i < centers->n; i++) {
      center_t *center = &oldcenters->center[i];

      if (basis_array) {
        errcod = ip_string("oldbasis",&bset,1,i);
        if (errcod!=IPE_OK) {
          fprintf(stderr,"scf_make_old_centers:  "
                         "could not read basis[%d]\n",i);
          goto the_place_where_errors_go;
        }
      }

      errcod = allocbn_center(center,"atom charge",
                              centers->center[i].atom,
                              centers->center[i].charge);
      if (errcod!=0) {
        fprintf(stderr,"scf_make_old_centers: could not alloc center %d\n",i);
        goto the_place_where_errors_go;
      }

      center->r[0]=centers->center[i].r[0];
      center->r[1]=centers->center[i].r[1];
      center->r[2]=centers->center[i].r[2];

      errcod = read_basis(sym_to_atom(center->atom),bset,&(center->basis));
      if (errcod != 0) {
        fprintf(stderr,"scf_make_old_centers: could not read basis %d\n",i);
        goto the_place_where_errors_go;
      }
    }

    int_normalize_centers(oldcenters);

    ip_cwk_pop();

    errcod=0;
    bcast0(&errcod,sizeof(int),mtype_get(),0);

  /* broadcast old centers struct to nodes */

    bcast0_centers(oldcenters,0,0);

   /* if we get here all is done */

    return 0;

   /* if we get here something went wrong */
the_place_where_errors_go:
    errcod=-1;
    bcast0(&errcod,sizeof(int),mtype_get(),0);

 /* the other nodes just sit around and wait to see what happened */
  } else {
    bcast0(&errcod,sizeof(int),mtype_get(),0);
  }

  if (errcod != 0) return -1;

 /* get oldcenters from node 0 */
  bcast0_centers(oldcenters,0,0);

  return 0;
}

/***************************************************************************
 * given the present basis and an old basis, project the vector in oldcenters
 * into the space of the new basis set and place the result in Scf_Vec
 *
 * input:
 *   centers = pointer to initialized centers struct
 *   scf_info = pointer to scf struct
 *   Scf_Vec  = column distributed matrix which contains core hamiltonion
 *              guess scf vector (for virtual orbitals)
 *   S        = scattered dmt matrix containing overlap integrals
 *   Sahalf   = column distributed matrix containing S-1/2
 *   oldvecfile = string containing path to old vector file
 *   oldcenters = pointer to initialized centers struct with old basis set
 *   outfile    = FILE pointer to output
 *
 * on return:
 *   Scf_Vec contains projected scf vector, schmidt orthogonalized
 *   
 * return 0 on success, -1 on failure
 */

GLOBAL_FUNCTION int
scf_project_vector(centers, scf_info, Scf_Vec, S, Sahalf, 
                   oldvecfile, oldcenters, outfile)
centers_t *centers;
scf_struct_t *scf_info;
dmt_matrix Scf_Vec;
dmt_matrix S;
dmt_matrix Sahalf;
char *oldvecfile;
centers_t *oldcenters;
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

  dmt_matrix Scr1=dmt_create("libdmtscf proj scr1",scf_info->nbfao,COLUMNS);
  dmt_matrix Scr2=dmt_create("libdmtscf proj scr2",scf_info->nbfao,COLUMNS);
  dmt_matrix Scr3=dmt_create("libdmtscf proj scr3",scf_info->nbfao,COLUMNS);
  dmt_matrix Cprime=dmt_create("libdmtscf proj scr4",scf_info->nbfao,COLUMNS);

  dmt_matrix Cnew,X,D,Sinv,C2,S2;

  assert(dmt_distribution(Scf_Vec) == COLUMNS);
  assert(dmt_distribution(Sahalf) == COLUMNS);
  assert(dmt_distribution(S) == SCATTERED);

  nocc = scf_info->nclosed + scf_info->nopen;

 /* form Sinv = (S**-1/2)**2 */

  Sinv=Scr1;
  dmt_copy(Sahalf,Scr2);
  dmt_mult(Scr2,Sahalf,Sinv);

 /* initialize the integral routines */
  int_initialize_1e(0,0,centers,oldcenters);
  int_initialize_offsets1(centers,oldcenters);

  m=oldcenters->nfunc;
  n=centers->nfunc;

 /* read old vector */

  C2=Scr2;
  dmt_fill(C2,0.0);

 /* node 0 opens the old vector file so that we can read it in */
  if (me==0) {
    if ((fd = open(oldvecfile,O_RDONLY)) < 0) {
      fprintf(stderr,"scf_project_vector: " 
                     "could not open file %s\n",oldvecfile);
      errcod=-1;
      bcast0(&errcod,sizeof(int),mtype_get(),0);
    } else {
      errcod=0;
      bcast0(&errcod,sizeof(int),mtype_get(),0);
    }
  } else {
    errcod=0;
    bcast0(&errcod,sizeof(int),mtype_get(),0);
  }

 /* if we couldn't open the vector file, return */
  if (errcod < 0) return -1;
 
 /* now pass the matrix for the old vector around on the loop.
  * node 0 will read in the columns of occupied MOs from oldvecfile.
  */
  loop = dmt_ngl_create("%m",C2);
  while (dmt_ngl_next(loop)) {
    int iind,isize,jsize;

    dmt_ngl_create_inner(loop,0);
    while(dmt_ngl_next_inner_m(loop,&iind,&isize,&i,&jsize,&col)) {
      if (i < nocc && me==0) {
        read(fd,col,sizeof(double)*m);
      }
    }
  }

  dmt_ngl_kill(loop);

  if (me==0) {
    close(fd);
    fprintf(outfile,"  read old scf vector from %s\n",oldvecfile);
  }

 /* form overlap between old basis and new basis */

  S2=Scr3;
  form_s2(centers,oldcenters,S2);

 /* form projected vector by multiplying the old vector by the overlap */

  dmt_mult(S2,C2,Cprime);

 /* free some memory */
  int_done_offsets1(centers,oldcenters);
  int_done_1e();

  free_centers(oldcenters);

 /* form D = Sinv x Cprime */

  D=Scr2;

  dmt_mult(Sinv,Cprime,D);

 /* form X = Cprime~ * D */

  X=Scr3;

  dmt_mult(Cprime,D,X);
  dmt_free(Cprime);

 /* form X**-1/2 */

  errcod = form_x_half(X,nocc);
  if (errcod!=0) {
    fprintf(stderr,"scf_project_vector: trouble in form_x_half\n");
    return -1;
  }

 /* finally, form Cnew = D*X**-1/2 */

  Cnew=Scr1;

  dmt_transpose(D);
  dmt_mult(D,X,Cnew);

 /* and put Cnew into Scf_Vec */

  for (i=0; i < nlocal ; i++) {
    double *ncol;

    dmt_get_col(Scf_Vec,i,&icol,&col);
    dmt_get_col(Cnew,i,&icol,&ncol);
    if (icol < nocc) for(j=0; j < n; j++) col[j]=ncol[j];
  }

 /* now orthogonalize Scf_Vec */

  errcod = scf_schmidt(scf_info,Scf_Vec,S,1);
  if (errcod !=0 ) {
    fprintf(stderr,"scf_project_vector: trouble orthogonalizing new vector\n");
    return -1;
  }

  dmt_free(Scr1);
  dmt_free(Scr2);
  dmt_free(Scr3);

  return 0;
}

/***************************************************************************
 * these are utility functions stolen out of libintv2
 */

#define COUNT 1
#define READ 2

static char *prb = "encountered a problem while reading the basis set:\n ";

LOCAL_FUNCTION int
read_basis(atom,basisname,basis)
char *atom;
char *basisname;
basis_t *basis;
{
  int errcod;

  /* Copy the basisname to the basis structure. */
  basis->name = (char *)malloc(strlen(basisname)+1);
  check_alloc(basis->name,"basis->name");
  strcpy(basis->name,basisname);

  /* If the basis name is "nothing" then return with an empty basis set. */
  if (!strcmp(basisname,"nothing")) {
    basis->n = 0;
    basis->shell = NULL;
    return IPE_OK;
  }

  if ((errcod = read_basis_(COUNT,atom,basisname,basis))!=IPE_OK) return errcod;

  basis->shell = (shell_t *) malloc(sizeof(shell_t)*basis->n);
  check_alloc(basis->shell,"basis->shell");
  basis->n = 0;
  return read_basis_(READ,atom,basisname,basis);
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
    } else if (errcod == IPE_HAS_NO_VALUE) {
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
    } else if (errcod == IPE_KEY_NOT_FOUND) {
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
    } else {
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


/*************************************************************************
 * these are utility functions used by scf_project_vector
 */

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
  if ((allocbn_double_matrix(&eigvec,"n1 n2",nocc,nocc)!=0) ||
      (allocbn_double_vector(&evals,"n",nocc) != 0)) {
    fprintf(stderr,"form_x_half: alloc of eigvec and evals failed\n");
    return -1;
  }

  if (allocbn_double_matrix(&x,"n1 n2",nocc,nocc) != 0) {
    fprintf(stderr,"form_x_half: alloc of x failed\n");
    return -1;
  }

 /* put X in local matrix */

  loop = dmt_ngl_create("%mr",X);
  while (dmt_ngl_next(loop)) {
    int iind,isize,jsize;
    double *col;

    dmt_ngl_create_inner(loop,0);
    while(dmt_ngl_next_inner_m(loop,&iind,&isize,&i,&jsize,&col)) {
      if (i < nocc) for(j=0; j < nocc; j++) x.d[i][j]=col[j];
    }
  }

  dmt_ngl_kill(loop);

 /* now form X**-1/2 */

  math_diag_dm(&x,&evals,&eigvec,1,1.0e-15,1);

  for (i=0; i < nocc; i++) evals.d[i] = 1.0/sqrt(evals.d[i]);

  for (i=0; i < nocc ; i++) {
    for (j=0; j < nocc; j++) {
      x.d[i][j]=0.0;
      for (k=0; k < nocc; k++) 
        x.d[i][j] += eigvec.d[i][k] * eigvec.d[j][k] * evals.d[k];
    }
  }

 /* put x back into X */
  for (i=0; i < nlocal; i++) {
    int icol;
    double *col;
    dmt_get_col(X,i,&icol,&col);
    if (icol < nocc) for (j=0; j < nocc; j++) col[j]=x.d[icol][j];
  }

  free_double_matrix(&x);
  free_double_matrix(&eigvec);
  free_double_vector(&evals);
  return(0);
}

LOCAL_FUNCTION VOID
form_s2(new_cen,oldcenters,S2)
centers_t *new_cen;
centers_t *oldcenters;
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

  for (i=0; i < nlocal; i++) ls2[i] = &s2col[i*msize];

 /* loop over old shells */
  for (oj=0; oj < oldcenters->nshell; oj++) {
    onfunc=INT_SH_NFUNC(oldcenters,oj);
    ofirst=oldcenters->func_num[oj];
    olast=ofirst+onfunc;

   /* is this shell held locally? */
    if ((lfirst >= ofirst && lfirst < olast) ||
        (ofirst >= lfirst && ofirst < llast) ||
        (olast > lfirst && olast <= llast) ||
        (llast > ofirst && llast <= olast)) {

     /* loop over new shells */
      for (ni=0; ni < new_cen->nshell; ni++) {
        nfirst=new_cen->func_num[ni];
        nlast=nfirst+INT_SH_NFUNC(new_cen,ni);

        int_shell_overlap(new_cen,oldcenters,block,ni,oj);

        ll=MAX0(0,lfirst-ofirst);
        for (l=MAX0(ofirst,lfirst); l < (MIN0(llast,olast)); ll++,l++) {
          for (kk=0,k=nfirst; k < nlast; kk++,k++) {
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

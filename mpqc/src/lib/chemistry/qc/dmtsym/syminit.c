
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <tmpl.h>
#include <math/array/math_lib.h>
#include <chemistry/qc/intv2/int_libv2.h>

#ifndef IOFF
#define IOFF(a,b) (a>b)?a*(a+1)/2+b:b*(b+1)/2+a
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include <chemistry/qc/dmtsym/symm.h>
#include <chemistry/qc/dmtsym/symm_mac.h>
#include <chemistry/qc/dmtsym/symmallc.h>
#include <chemistry/qc/dmtsym/symmfree.h>
#include <chemistry/qc/dmtsym/symmzero.h>

#include <chemistry/qc/dmtsym/mkrpd.gbl>

#include <chemistry/qc/dmtsym/syminit.gbl>
#include <chemistry/qc/dmtsym/syminit.lcl>

/***************************************************************************
 *
 * This initializes centers based on the info found in unique_centers.
 * It does not read any input.
 *
 * input:
 *   unique_centers = pointer to centers struct containing unique atoms
 *   centers        = pointer to uninitialized centers struct
 *   sym_info       = pointer to uninitialized sym_struct
 *   point_group    = string containing schoenflies symbol
 *
 * on return:
 *   centers contains all redundant atoms
 *   sym_info contains symmetry information
 *
 * returns 0 on success, -1 on failure
 */

GLOBAL_FUNCTION int
sym_init_given_centers(unique_centers,centers,sym_info,point_group)
centers_t *unique_centers;
centers_t *centers;
sym_struct_t *sym_info;
char *point_group;
{
  int g,n,nirr;
  int errcod;
  enum pgroups pg;
  double_array3_t trans;

 /* the first order of business is to find out what point group we're
  * dealing with
  */

  if (!point_group || point_group[0] == 0) {
    fprintf(stderr,"sym_init_centers: no symmetry specified\n");
    return(-1);
  }

/* now parse the point group symbol, this will give us the order of the
 * point group(g), the type of point group (pg), the order of the principle
 * rotation axis (n), and the number of irreps (nirr)
 */

  errcod = sym_parse_symbol(point_group,&g,&pg,&n,&nirr);
  if (errcod != 0) {
    fprintf(stderr,"sym_init_given_centers: "
                   "could not parse point group symbol %s\n",point_group);
    return -1;
  }

 /* set up atomic transformation matrices */

  errcod = allocbn_double_array3(&trans,"n1 n2 n3",g,3,3);
  if (errcod != 0) {
    fprintf(stderr,"sym_init_given_centers: "
                  "could not allocate memory for trans\n");
    return -1;
  }
  zero_double_array3(&trans);

  errcod = make_trans(&trans,n,pg);
  if(errcod != 0) {
    fprintf(stderr,"sym_init_given_centers: cannot handle %s yet",point_group);
    return -1;
  }

 /* let's take a stab at generating the non-unique atoms, etc */

  sym_redund_centers_from_unique(unique_centers,centers,&trans);

 /* now initialize the full centers struct */

  sym_make_sym_struct(centers,sym_info,point_group,&trans);

  free_double_array3(&trans);

  return(0);
}

/**********************************************************************
 *
 * given two double vectors containing x, y, and z coordinates, return
 * sqrt(delta_x^2 + delta_y^2 + delta_z^2)
 *
 */

LOCAL_FUNCTION double
dist(a,b)
double *a;
double *b;
{
  double x,y,z;
  return sqrt((x=(a[0]-b[0]))*x +
              (y=(a[1]-b[1]))*y +
              (z=(a[2]-b[2]))*z);
}

/**********************************************************************
 *
 * this takes a pointer to a centers struct and a pointer to a double vector
 * containing the x, y, and z coordinates of an atom.  it returns the
 * index of the center in c which has the same coordinates as tr, or
 * returns -1 if no atom in c matchs tr
 */

LOCAL_FUNCTION int
comp_r(c,tr)
centers_t *c;
double *tr;
{
  int i;

  for (i=0; i < c->n ; i++)
    if (dist(c->center[i].r,tr) < 0.1) return (i);

  return -1;
}

/**************************************************************************
 *
 * given an old set of coordinates r, and a transformation matrix tr,
 * return the transformed coordinates in tr
 */

LOCAL_FUNCTION void
new_r(r,tr,tm)
double *r;
double *tr;
double **tm;
{
  int i,j;

  for (i=0; i < 3 ; i++) {
    tr[i]=0.0;
    for(j=0; j < 3 ; j++) tr[i] += tm[i][j]*r[j];
  }
}

/**************************************************************************
 * 
 * given a set of transformation matrices, and a set of cartesian coordinates,
 * find the number of times r maps into itself, and then return the
 * number of atoms equivalent to the one at r.
 */

LOCAL_FUNCTION int
n_eq(tm,r)
double_array3_t *tm;
double *r;
{
  int g;
  int nmap;
  double tr[3];

  for (nmap=g=0; g < tm->n1 ; g++) {
    new_r(r,tr,tm->d[g]);
    if (dist(r,tr) < 0.1) nmap++;
  }

  return(tm->n1/nmap);
}

/*************************************************************************
 *
 * this function will produce a set of transformation matrices describing
 * all of the symmetry operations in a given point group
 *
 * input:
 *   trans = pointer to an initialized double_array_3
 *   n     = the order of the principal rotation axis
 *   pg    = the generic type of the point group (eg Cnh, Dnd, etc...)
 *   
 * on return:
 *   trans contains the transformation matrices
 *
 * return 0 on success, -1 on failure
 */

LOCAL_FUNCTION int
make_trans(trans, n, pg)
double_array3_t *trans;
int n;
enum pgroups pg;
{
  int i,j;
  double theta;

  theta = (double) 2.0*M_PI/n;

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
    for (i=0; i < n ; i++) {
      trans->d[i][0][0] = trans->d[i][1][1] = cos((double)theta*i);
      trans->d[i][0][1] = sin((double)theta*i);
      trans->d[i][1][0] = -sin((double)theta*i);
      trans->d[i][2][2] = 1.0;
    }
    break;
  case _PG_CNV:
    for (i=0; i < n ; i++) {
      trans->d[i][0][0] = trans->d[i][1][1] = cos((double)theta*i);
      trans->d[i][0][1] = sin((double)theta*i);
      trans->d[i][1][0] = -sin((double)theta*i);
      trans->d[i][2][2] = 1.0;
    }

    for (i=0,j=n; i < n ; i++,j++) {
      trans->d[j][0][0] = cos((double)theta*i);
      trans->d[j][1][1] = -cos((double)theta*i);
      trans->d[j][1][0] = trans->d[j][0][1] = sin((double)theta*i);
      trans->d[j][2][2] = 1.0;
    }
    break;
  case _PG_CNH:
    for (i=0; i < n ; i++) {
      trans->d[i][0][0] = trans->d[i][1][1] =
      trans->d[i+n][0][0] = trans->d[i+n][1][1] = cos((double)theta*i);
      trans->d[i][0][1] = trans->d[i+n][0][1] = sin((double)theta*i);
      trans->d[i][1][0] = trans->d[i+n][1][0] = -sin((double)theta*i);
      trans->d[i][2][2] = 1.0;
      trans->d[i+n][2][2] = -1.0;
    }
    break;
  case _PG_SN:
    for (i=0; i < n ; i++) {
      trans->d[i][0][0] = trans->d[i][1][1] = cos((double)theta*i);
      trans->d[i][0][1] = sin((double)theta*i);
      trans->d[i][1][0] = -sin((double)theta*i);
      trans->d[i][2][2] = pow(-1.0,(double) i);
    }
    break;
  case _PG_DN:
    for (i=0; i < n ; i++) {
      trans->d[i][0][0] = trans->d[i][1][1] = cos((double) theta*i);
      trans->d[i][0][1] = sin((double) theta*i);
      trans->d[i][1][0] = -sin((double) theta*i);
      trans->d[i][2][2] = 1.0;
    }

    for (i=0,j=n; i < n ; i++,j++) {
      trans->d[j][0][0] = cos((double)theta*i);
      trans->d[j][1][1] = -cos((double)theta*i);
      trans->d[j][1][0] = trans->d[j][0][1] = -sin((double)theta*i);
      trans->d[j][2][2] = -1.0;
    }
    break;
  case _PG_DND:
    for (i=0; i < 2*n ; i++) {
      trans->d[i][0][0] = trans->d[i][1][1] = cos((double)theta*i*0.5);
      trans->d[i][0][1] = sin((double)theta*i*0.5);
      trans->d[i][1][0] = -sin((double)theta*i*0.5);
      trans->d[i][2][2] = pow(-1.0,(double) i);
    }

    for (i=0,j=2*n; i < n ; i++,j++) {
      trans->d[j][0][0] = cos((double)theta*i);
      trans->d[j][1][1] = -cos((double)theta*i);
      trans->d[j][1][0] = trans->d[j][0][1] = -sin((double)theta*i);
      trans->d[j][2][2] = -1.0;
    }

    for (i=0,j=3*n; i < n ; i++,j++) {
      trans->d[j][0][0] = cos((double)theta*i+theta*0.5);
      trans->d[j][1][1] = -cos((double)theta*i+theta*0.5);
      trans->d[j][1][0] = trans->d[j][0][1] = -sin((double)theta*i+theta*0.5);
      trans->d[j][2][2] = 1.0;
    }
    break;
  case _PG_DNH:
    for (i=0; i < n ; i++) {
      trans->d[i][0][0] = trans->d[i][1][1] =
      trans->d[i+n][0][0] = trans->d[i+n][1][1] = cos((double)theta*i);
      trans->d[i][0][1] = trans->d[i+n][0][1] = sin((double)theta*i);
      trans->d[i][1][0] = trans->d[i+n][1][0] = -sin((double)theta*i);
      trans->d[i][2][2] = 1.0;
      trans->d[i+n][2][2] = -1.0;
    }

    for (i=0,j=2*n; i < n ; i++,j++) {
      trans->d[j][0][0] = cos((double)theta*i);
      trans->d[j][1][1] = -cos((double)theta*i);
      trans->d[j][1][0] = trans->d[j][0][1] = -sin((double)theta*i);
      trans->d[j][2][2] = -1.0;
    }

    for (i=0,j=3*n; i < n ; i++,j++) {
      trans->d[j][0][0] = cos((double)theta*i);
      trans->d[j][1][1] = -cos((double)theta*i);
      trans->d[j][1][0] = trans->d[j][0][1] = sin((double)theta*i);
      trans->d[j][2][2] = 1.0;
    }
    break;
  default:
    return -1;
  }

  return 0;
}

/************************************************************************
 *
 * given the schoenflies symbol of a point group, return various information
 * about that point group
 *
 * input:
 *   point_group = string containing the schoenflies symbol
 *   g           = the order of the group
 *   pg          = the generic point group type
 *   n           = the order of the principal rotation axis
 *   nirr        = the number of irreps
 *
 * on return:
 *   g, pg, n, and nirr contain what they should
 *
 * return 0 on success, -1 on failure
 */

GLOBAL_FUNCTION int
sym_parse_symbol(point_group,g,pg,n,nirr)
char *point_group;
int *g;
enum pgroups *pg;
int *n;
int *nirr;
{
  *n = 1;
  *g = 1;
  *pg = _PG_C1;

  if (!strcmp(point_group,"ci")) {
    *g = 2;
    *pg = _PG_CI;
    *nirr = 2;
  } else if (!strcmp(point_group,"c1")) {
    *g = 1;
    *pg = _PG_C1;
    *nirr = 1;
  } else if (!strcmp(point_group,"cs")) {
    *g = 2;
    *pg = _PG_CS;
    *nirr = 2;
  }

 /* the cn groups */

  else if (point_group[0]=='c') {
    int nab,ne;

   /* n will be the second element of point_group */
    if (point_group[1] == 0) {
      goto parse_sym_err;
    } else {
      *n = atoi(&point_group[1]);
      ne = (*n%2) ? *n/2 : *n/2-1;
      nab = (*n%2) ? 1 : 2;
    }

   /* check to see if this is cnv or cnh */
    if (point_group[2] != 0) {
      *g  = 2 * (*n);

      if (point_group[2] == 'v') {
        *pg = _PG_CNV;
        *nirr = 2*nab + ne;
      } else if(point_group[2] == 'h') {
        *pg = _PG_CNH;
        *nirr = 2*(nab+ne);
      } else {
        goto parse_sym_err;
      }
    } else {
      *g = *n;
      *pg = _PG_CN;
      *nirr = nab+ne;
    }
  }

 /* the dn groups */
  else if (point_group[0]=='d') {
    int nab,ne;

    if (point_group[1] == 0)
      goto parse_sym_err;

    *n = atoi(&point_group[1]);
    ne = (*n%2) ? *n/2 : *n/2-1;
    nab = (*n%2) ? 1 : 2;

    if (point_group[2] != 0) {
      *g = 4 * (*n);

      if (point_group[2] == 'd') {
        *pg = _PG_DND;
        *nirr = *n+3;
      } else if (point_group[2] == 'h') {
        *pg = _PG_DNH;
        *nirr = 4*nab + 2*ne;
      } else {
        goto parse_sym_err;
      }
    } else {
      *g = 2 * (*n);
      *pg = _PG_DN;
      *nirr = 2*nab + ne;
    }
  }

 /* the sn groups */
  else if (point_group[0]=='s') {
    if (point_group[1] == 0)
      goto parse_sym_err;

  /* only S2n groups make sense */
    *n = atoi(&point_group[1]);
    if((*n)%2)
      goto parse_sym_err;

    *g = *n;
    *pg=_PG_SN;
    *nirr = *n/2+1;
    }

 /* hard wired tetrahedral groups */

  else if (point_group[0]=='t') {
    if (point_group[1] != 0) {
      if (point_group[1] == 'd') {
        *g = 24;
        *pg = _PG_TD;
        *nirr = 5;
      } else if(point_group[1] == 'h') {
        *g = 24;
        *pg = _PG_TH;
        *nirr = 6;
      } else {
        goto parse_sym_err;
      }
    } else {
      *g = 12;
      *pg = _PG_T;
      *nirr = 3;
    }
  }

 /* hard wired octahedral groups */

  else if (point_group[0]=='o') {
    if (point_group[1] != 0) {
      if (point_group[1] == 'h') {
        *pg = _PG_OH;
        *g = 48;
        *nirr = 10;
      } else {
        goto parse_sym_err;
      }
    } else {
      *g = 24;
      *pg = _PG_O;
      *nirr = 5;
    }
  }

 /* hard wired icosahedral groups */

  else if (point_group[0]=='i') {
    if (point_group[1] != 0) {
      if (point_group[1] == 'h') {
        *g = 120;
        *pg = _PG_IH;
        *nirr = 10;
      } else {
        goto parse_sym_err;
      }
    } else {
      *g = 60;
      *pg = _PG_I;
      *nirr = 5;
    }
  }

 /* this is bad */
  else {
    goto parse_sym_err;
  }

 /* this is good */
  return 0;

parse_sym_err:
  fprintf(stderr,"unknown Schoenflies symbol: %s",point_group);
  return -1;

}

/************************************************************************
 *
 * given the unique atoms in unique_centers, form all atoms in centers
 *
 * input:
 *   unique_centers = pointer to centers struct containing unique atoms
 *   centers        = uninitialized centers struct
 *   trans          = double_array_3 containing transformation matrices
 *                      for all symmetry operations of the point group
 *
 * on return:
 *   centers contains redundant atoms
 *
 * return 0 on success, -1 on failure
 */

GLOBAL_FUNCTION int
sym_redund_centers_from_unique(unique_centers, centers, trans)
centers_t *unique_centers;
centers_t *centers;
double_array3_t *trans;
{
  int i,j,ng;
  int errcod;
  int atom, totnat;
  int nat = unique_centers->n;
  int *n_eq_atoms;

  double tr[3];
  double tol = 1.0e-10;

  if (!nat) {
    fprintf(stderr,"sym_redund_centers_from_unique: "
                   "there don\'t appear to be any atoms\n");
    return -1;
  }

  n_eq_atoms = (int *) malloc(sizeof(int)*nat);
  if (!n_eq_atoms) {
    fprintf(stderr,"sym_redund_centers_from_unique: "
                   "could not allocate memory for n_eq_atoms");
    return -1;
  }

 /* loop over nat, and apply each symop to r, test to see if new */
  for (i=totnat=0; i < nat ; i++) {
    n_eq_atoms[i] = n_eq(trans,unique_centers->center[i].r);
    totnat += n_eq_atoms[i];
  }

 /* now that we know how many atoms we have, we can generate their coordinates
  * and get the full centers struct.
  */

  errcod = allocbn_centers(centers,"n",totnat);
  if(errcod != 0) {
    fprintf(stderr,"sym_redund_centers_from_unique"
                   "could not allocate memory for centers");
    return -1;
  }

  for (i=atom=0; i < nat ; i++) {
    for (ng=0; ng < n_eq_atoms[i]; ng++,atom++) {
      init_center(&centers->center[atom]);
      errcod = assign_center(&centers->center[atom],
                             &unique_centers->center[i]);
      if (errcod != 0) {
        fprintf(stderr,"sym_redund_centers_from_unique"
                       "trouble making centers->center[%d]",atom);
        return -1;
      }

      for (j=0; j < 3; j++) centers->center[atom].r[j] =
                                 unique_centers->center[i].r[j];
    }
  }

  for (i=atom=0; i < nat ; i++) {
    for (j=0; j < 3; j++) centers->center[atom].r[j] = tr[j] =
                                              unique_centers->center[i].r[j];
    atom++;

    if (!(fabs(tr[0]) < tol && fabs(tr[1]) < tol && fabs(tr[2]) < tol)) {
      for (ng=1; ng < trans->n1; ng++) {
        new_r(unique_centers->center[i].r,tr,trans->d[ng]);
        if (comp_r(centers,tr) < 0 ) {
          for (j=0; j < 3; j++) centers->center[atom].r[j]=tr[j];
          atom++;
        }
      }
    }
  }

  int_initialize_centers(centers);

  return 0;
}

/************************************************************************
 *
 * given a center struct and a set of transformation matrices, construct
 * the sym_struct
 *
 * input:
 *   centers  = centers struct containing all atoms
 *   sym_info = uninitialized sym_struct
 *   trans    = pointer to double_array_3 containing transformation matrices
 *               for all symmetry operations of the point group
 *
 * on return:
 *   sym_info contains all symmetry information
 *
 * return 0 on success, -1 on failure
 */

GLOBAL_FUNCTION int
sym_make_sym_struct(centers,sym_info,point_group,trans)
centers_t *centers;
sym_struct_t *sym_info;
char *point_group;
double_array3_t *trans;
{
  int i,j,nsh,errcod;
  int g,n,nirr,ns,nb,gc,ng,nij;
  int f_exist=0, g_exist=0, h_exist=0;
  int *first,*last;
  enum pgroups ptgrp_sym;
  double tr[3];

 /* we need some more information about the basis set, etc */
  int_initialize_1e(0,0,centers,centers);
  int_initialize_offsets1(centers,centers);

  nsh = centers->nshell;

 /* let's see if there are f or g functions. we'll puke if h or higher are
  * used.
  */

  for (i=0; i < nsh ; i++) {
    j = centers->center[centers->center_num[i]].basis.
        shell[centers->shell_num[i]].type[0].am;

    if (j==3) f_exist=1;
    else if (j==4) g_exist=1;
    else if (j >= 5) h_exist=1;
  }

  if (h_exist) {
    fprintf(stderr,"sym_make_sym_struct: cannot handle am greater than 4\n");
    return -1;
  }

 /* parse the point group symbol */
  errcod = sym_parse_symbol(point_group,&g,&ptgrp_sym,&n,&nirr);
  if (errcod != 0) {
    fprintf(stderr,"sym_make_sym_struct: "
                   "could not parse point group symbol %s\n",point_group);
    return -1;
  }

 /* and allocate memory for sym_info */

  if (g==1) {
    errcod =
      allocbn_sym_struct(sym_info,"point_group pg",point_group,ptgrp_sym);

  } else if (g_exist) {
    errcod = allocbn_sym_struct(sym_info,
         "natom g nshell nshtri point_group pg nirrep n1 n2 n2p n3 n3p n4 n4p",
         centers->n, g, nsh, nsh*(nsh+1)/2, point_group, ptgrp_sym, nirr,
         3, 6, 5, 10, 7, 15, 9);

  } else if (f_exist) {
    errcod = allocbn_sym_struct(sym_info,
         "natom g nshell nshtri point_group pg nirrep n1 n2 n2p n3 n3p",
         centers->n , g, nsh, nsh*(nsh+1)/2, point_group, ptgrp_sym, nirr,
         3, 6, 5, 10, 7);

  } else {
    errcod = allocbn_sym_struct(sym_info,
                       "natom g nshell nshtri point_group pg nirrep n1 n2 n2p",
                       centers->n, g, nsh, nsh*(nsh+1)/2, point_group,
                       ptgrp_sym, nirr, 3, 6, 5);
  }

  if (errcod != 0) {
    fprintf(stderr,"sym_make_sym_struct:  cannot alloc sym_struct\n");
    return -1;
  }

  zero_sym_struct(sym_info);

  if (g > 1) {
    errcod = allocbn_char_tab(&sym_info->ct,"nirrep",nirr);
    if (errcod != 0) {
      fprintf(stderr,"sym_make_sym_struct: "
                     "could not allocate memory for char_tab");
      return -1;
    }

    first = (int *) malloc(sizeof(int)*centers->n);
    last = (int *) malloc(sizeof(int)*centers->n);

    if (!first || !last) {
      fprintf(stderr,"sym_make_sym_struct: cannot alloc first or last\n");
      return -1;
    }

    ns=nb=0;
    for (i=0; i < centers->n ; i++) {
      first[i] = ns;
      last[i]  = ns + centers->center[i].basis.n;
      ns += centers->center[i].basis.n;

      for (ng=0; ng < g ; ng++) {
        new_r(centers->center[i].r,tr,trans->d[ng]);
        sym_info->atom_map[i][ng] = comp_r(centers,tr);
      }

      sym_info->first[i] = nb;
      for (j=0; j < centers->center[i].basis.n ; j++) {
        nb += centers->center[i].basis.shell[j].nfunc;
      }
    }

    for (i=nb=0; i < centers->n ; i++) {
      sym_info->firstp[i] = nb;
      for (j=0; j < centers->center[i].basis.n ; j++) {
        for (gc=0; gc < centers->center[i].basis.shell[j].ncon ; gc++) {
          ns = centers->center[i].basis.shell[j].type[gc].am;
          nb += 2*ns+1;
        }
      }
    }

    for (i=0; i < centers->n ; i++) {
      for (j=0,ns=first[i]; ns < last[i]; j++,ns++) {
        for (ng=0; ng < g ; ng++) {
          int atom = sym_info->atom_map[i][ng];
          sym_info->shell_map[ns][ng] = first[atom]+j;
        }
      }
    }

   /* now let's set up lambda from Dupuis and King's paper */
  
    for (i=0; i < nsh ; i++) {
      int leave;

      for (ng=leave=0; ng < g ; ng++)
        if(sym_info->shell_map[i][ng] > i) leave=1;

      if(leave) continue;

      sym_info->p1[i]=1;

      for (j=0; j <=i ; j++) {
        int ij=IOFF(i,j);
        nij=0;
        for (ng=leave=0; ng < g ; ng++) {
          int gi = sym_info->shell_map[i][ng];
          int gj = sym_info->shell_map[j][ng];
          int gij = IOFF(gi,gj);
          if (gij > ij) leave=1;
          if (gij == ij) nij++;
        }
        if (leave) continue;
        sym_info->lamij[ij] = (char) (g/nij);
      }
    }

   /* finally, let's form Rp and Rd in sym_info */

    errcod = sym_make_rp_d(sym_info,n,ptgrp_sym,f_exist,g_exist);
    if (errcod != 0) {
      fprintf(stderr,"sym_make_sym_struct: cannot handle %s yet",point_group);
      return -1;
    }

    free(first);
    free(last);
  }

  int_done_1e();
  int_done_offsets1(centers,centers);

  return 0;
}

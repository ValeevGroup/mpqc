
#if defined(__GNUC__)
#pragma implementation
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

extern "C" {
#include <math/array/matrix.h>
#include <math/array/matrixzero.h>
#include <math/array/matrixprnt.h>
#include <math/array/matrixallc.h>
}

#include <chemistry/qc/intv2/int_macros.h>
#include <chemistry/qc/intv2/tformv2.h>
#include <chemistry/qc/intv2/integralv2.h>

#include <chemistry/qc/intv2/utils.gbl>

#define PRINT 0

static double *source = 0;
static int nsourcemax = 0;

static void
do_copy1(double *source, double *target, int chunk,
         int n1, int s1, int offset1,
         int n2, int s2, int offset2)
{
  int i1, i2;

  for (i1=0; i1<n1; i1++) {
      int off = ((offset1 + i1)*s2 + offset2)*chunk;
      for (i2=0; i2<n2*chunk; i2++, off++) {
          target[off] = source[off];
        }
    }
}

static void
do_copy2(double *source, double *target,
         int n1, int s1, int offset1,
         int n2, int s2, int offset2,
         int n3, int s3, int offset3,
         int n4, int s4, int offset4)
{
  int i1, i2, i3, i4;

  for (i1=0; i1<n1; i1++) {
      for (i2=0; i2<n2; i2++) {
          for (i3=0; i3<n3; i3++) {
              int off = (((offset1 + i1)*s2 + offset2 + i2)
                         *s3 + offset3 + i3)*s4 + offset4;
              for (i4=0; i4<n4; i4++, off++) {
                  target[off] = source[off];
                }
            }
        }
    }
}

static void
do_sparse_transform11(double *source, double *target, int chunk,
                      SphericalTransformIter& trans,
                      int offsetcart1,
                      int offsetpure1,
                      int n2, int s2, int offset2)
{
  int i2;

  for (trans.begin(); trans.ready(); trans.next()) {
      double coef = trans.coef();
      int pure = trans.pureindex();
      int cart = trans.cartindex();
      int offtarget = ((offsetpure1 + pure)*s2 + offset2)*chunk;
      int offsource = ((offsetcart1 + cart)*s2 + offset2)*chunk;
      for (i2=0; i2<n2*chunk; i2++) {
          target[offtarget++] += coef * source[offsource++];
        }
    }
}

static void
do_sparse_transform12(double *source, double *target, int chunk,
                      SphericalTransformIter& trans,
                      int n1, int offset1,
                      int s2cart, int offsetcart2,
                      int s2pure, int offsetpure2)
{
  int i1, ichunk;

  for (trans.begin(); trans.ready(); trans.next()) {
      double coef = trans.coef();
      int pure = trans.pureindex();
      int cart = trans.cartindex();
      for (i1=0; i1<n1; i1++) {
          int offtarget = ((offset1 + i1)*s2pure + offsetpure2 + pure)*chunk;
          int offsource = ((offset1 + i1)*s2cart + offsetcart2 + cart)*chunk; 
          for (ichunk=0; ichunk<chunk; ichunk++) {
              target[offtarget++] += coef * source[offsource++];
            }
        }
    }
}

void
do_abort()
{
  abort();
}


void
do_sparse_transform2(double *source, double *target,
                     int index, SphericalTransformIter& trans,
                     int stcart, int stpure,
                     int ogctcart, int ogctpure,
                     int n1, int s1, int ogc1,
                     int n2, int s2, int ogc2,
                     int n3, int s3, int ogc3,
                     int n4, int s4, int ogc4)
{
  int i1, i2, i3, i4;
  int offtarget, offsource;

  switch (index) {
  case 0:
      n1 = 1;
      s1 = 1;
      break;
  case 1:
      n2 = 1;
      s2 = 1;
      break;
  case 2:
      n3 = 1;
      s3 = 1;
      break;
  case 3:
      n4 = 1;
      s4 = 1;
      break;
    }

  for (trans.begin(); trans.ready(); trans.next()) {
      double coef = trans.coef();
      int pure = trans.pureindex();
      int cart = trans.cartindex();
      for (i1=0; i1<n1; i1++) {
          for (i2=0; i2<n2; i2++) {
              for (i3=0; i3<n3; i3++) {
                  for (i4=0; i4<n4; i4++) {
                      if (index == 0) {
                          offtarget = ogctpure + pure;
                          offsource = ogctcart + cart;
                        }
                      else {
                          offtarget = ogc1 + i1;
                          offsource = ogc1 + i1;
                        }
                      if (index == 1) {
                          offtarget = offtarget*stpure + ogctpure + pure;
                          offsource = offsource*stcart + ogctcart + cart;
                        }
                      else {
                          offtarget = offtarget*s2 + ogc2 + i2;
                          offsource = offsource*s2 + ogc2 + i2;
                        }                
                      if (index == 2) {  
                          offtarget = offtarget*stpure + ogctpure + pure;
                          offsource = offsource*stcart + ogctcart + cart;
                        }                
                      else {             
                          offtarget = offtarget*s3 + ogc3 + i3;
                          offsource = offsource*s3 + ogc3 + i3;
                        }                
                      if (index == 3) {  
                          offtarget = offtarget*stpure + ogctpure + pure;
                          offsource = offsource*stcart + ogctcart + cart;
                        }                
                      else {             
                          offtarget = offtarget*s4 + ogc4 + i4;
                          offsource = offsource*s4 + ogc4 + i4;
                        }
                      target[offtarget] += coef * source[offsource];
#if PRINT
                      if (fabs(coef * source[offsource])>1.0e-15) {
                          printf("%3.1f * %15.11f [%d] += [%d] -> %15.11f\n",
                                 coef, source[offsource], offsource,
                                 offtarget, target[offtarget]);
                        }
#endif
                    }
                }
            }
        }
    }
}

/* make sure enough space exists for the source integrals */
static void
source_space(int nsource)
{
  if (nsourcemax < nsource) {
      if (source) free(source);
      source = (double*) malloc(sizeof(double)*nsource);
      nsourcemax = nsource;
    }
}

/* Compute the number of cartesian functions in a shell. */
static int
have_pure(shell_t *sh)
{
  int i;
  for (i=0; i<sh->ncon; i++) {
      if (sh->type[i].puream) return 1;
    }
  return 0;
}

static void
copy_to_source(double *integrals, int nsource)
{
  int i;
  double *tmp, *tmp2;

  /* Allocate more temporary space if needed. */
  source_space(nsource);

  tmp = source;
  tmp2 = integrals;
  for (i=0; i<nsource; i++) *tmp++ = *tmp2++;
}

static void
do_transform_1e(double *integrals, shell_t *sh1, shell_t *sh2, int chunk)
{
  int i, j;
  int ogc1, ogc2;
  int ogc1pure, ogc2pure;
  int am1, am2;
  int pure1 = have_pure(sh1);
  int pure2 = have_pure(sh2);
  int ncart1 = int_ncart(sh1);
  int ncart2 = int_ncart(sh2);
  int nfunc1 = sh1->nfunc;
  int nfunc2 = sh2->nfunc;
  int nfunci, nfuncj;

  if (!pure1 && !pure2) return;

  /* Loop through the generalized general contractions,
   * transforming the first index. */
  if (pure1) {
      copy_to_source(integrals, ncart1*ncart2*chunk);
      memset(integrals, 0, sizeof(double)*sh1->nfunc*ncart2*chunk);

      ogc1 = 0;
      ogc1pure = 0;
      for (i=0; i<sh1->ncon; i++) {
          am1 = sh1->type[i].am;
          nfunci = INT_NFUNC(sh1->type[i].puream, sh1->type[i].am);
          ogc2 = 0;
          for (j=0; j<sh2->ncon; j++) {
              am2 = sh2->type[j].am;
              nfuncj = INT_NFUNC(sh2->type[j].puream, sh2->type[j].am);

              if (sh1->type[i].puream) {
                  SphericalTransformIterV2 trans(sh1->type[i].am);
                  do_sparse_transform11(source, integrals, chunk,
                                        trans,
                                        ogc1,
                                        ogc1pure,
                                        INT_NCART(am2), ncart2, ogc2);
                }
              else {
                  do_copy1(source, integrals, chunk,
                           nfunci, nfunc1, ogc1pure,
                           INT_NCART(am2), ncart2, ogc2);
                }
              ogc2 += INT_NCART(am2);
            }
          ogc1 += INT_NCART(am1);
          ogc1pure += INT_NPURE(am1);
        }
    }

  if (pure2) {
      copy_to_source(integrals, nfunc1*ncart2*chunk);
      memset(integrals, 0, sizeof(double)*sh1->nfunc*sh2->nfunc*chunk);

      ogc1 = 0;
      for (i=0; i<sh1->ncon; i++) {
          am1 = sh1->type[i].am;
          nfunci = INT_NFUNC(sh1->type[i].puream, sh1->type[i].am);
          ogc2 = 0;
          ogc2pure = 0;
          for (j=0; j<sh2->ncon; j++) {
              am2 = sh2->type[j].am;
              nfuncj = INT_NFUNC(sh2->type[j].puream, sh2->type[j].am);

              if (sh2->type[j].puream) {
                  SphericalTransformIterV2 trans(sh2->type[j].am);
                  do_sparse_transform12(source, integrals, chunk,
                                        trans,
                                        INT_NPURE(am1), ogc1,
                                        ncart2, ogc2,
                                        sh2->nfunc, ogc2pure);
                }
              else {
                  do_copy1(source, integrals, chunk,
                           nfunci, nfunc1, ogc1,
                           nfuncj, nfunc2, ogc2pure);
                }
              ogc2 += INT_NCART(am2);
              ogc2pure += INT_NPURE(am2);
            }
          ogc1 += INT_NPURE(am1);
        }
    }
}

/* it is ok for integrals and target to overlap */
void
transform_1e(double *integrals, double *target,
             shell_t *sh1, shell_t *sh2, int chunk)
{
  int ntarget;

  do_transform_1e(integrals, sh1, sh2, chunk);

  /* copy the integrals to the target, if necessary */
  ntarget = sh1->nfunc * sh2->nfunc;
  if (integrals != target) {
      memmove(target, integrals, ntarget*sizeof(double)*chunk);
    }
}

/* it is not ok for integrals and target to overlap */
void
accum_transform_1e(double *integrals, double *target,
                   shell_t *sh1, shell_t *sh2, int chunk)
{
  int i, ntarget;

  do_transform_1e(integrals, sh1, sh2, chunk);

  /* accum the integrals to the target */
  ntarget = sh1->nfunc * sh2->nfunc * chunk;
  for (i=0; i<ntarget; i++) target[i] += integrals[i];
}

void
int_transform_1e(double *integrals, double *target,
             shell_t *sh1, shell_t *sh2)
{
  transform_1e(integrals, target, sh1, sh2, 1);
}

void
int_accum_transform_1e(double *integrals, double *target,
                   shell_t *sh1, shell_t *sh2)
{
  accum_transform_1e(integrals, target, sh1, sh2, 1);
}

void
int_transform_1e_xyz(double *integrals, double *target,
             shell_t *sh1, shell_t *sh2)
{
  transform_1e(integrals, target, sh1, sh2, 3);
}

void
int_accum_transform_1e_xyz(double *integrals, double *target,
                   shell_t *sh1, shell_t *sh2)
{
  accum_transform_1e(integrals, target, sh1, sh2, 3);
}

void
do_gencon_sparse_transform_2e(double *integrals, double *target, int index,
                              shell_t *sh1, shell_t *sh2,
                              shell_t *sh3, shell_t *sh4)
{
  int ncart[4];
  int nfunc[4];
  int nfunci, nfuncj, nfunck, nfuncl;
  int ncarti, ncartj, ncartk, ncartl;
  int i, j, k, l;
  int ogccart[4];
  int ogcfunc[4];
  int am1, am2, am3, am4;

  int ntarget1;
  int ntarget2;
  int ntarget3;
  int ntarget4;
              
  int nsource1;
  int nsource2;
  int nsource3;
  int nsource4;

  int *ni = &ncarti;
  int *nj = &ncartj;
  int *nk = &ncartk;
  int *nl = &ncartl;

  int *ogc1 = &ogccart[0];
  int *ogc2 = &ogccart[1];
  int *ogc3 = &ogccart[2];
  int *ogc4 = &ogccart[3];

  shell_t *shell;

  int *tgencon;

  ncart[0] = int_ncart(sh1);
  ncart[1] = int_ncart(sh2);
  ncart[2] = int_ncart(sh3);
  ncart[3] = int_ncart(sh4);
  nfunc[0] = sh1->nfunc;
  nfunc[1] = sh2->nfunc;
  nfunc[2] = sh3->nfunc;
  nfunc[3] = sh4->nfunc;

  ntarget1 = ncart[0];
  ntarget2 = ncart[1];
  ntarget3 = ncart[2];
  ntarget4 = ncart[3];
  
  nsource1 = ncart[0];
  nsource2 = ncart[1];
  nsource3 = ncart[2];
  nsource4 = ncart[3];

  if (index >= 0) {
      ntarget1 = nfunc[0];
      if (index >= 1) {
          ntarget2 = nfunc[1];
          nsource1 = nfunc[0];
          ni = &nfunci;
          ogc1 = &ogcfunc[0];
          if (index >= 2) {
              ntarget3 = nfunc[2];
              nsource2 = nfunc[1];
              nj = &nfuncj;
              ogc2 = &ogcfunc[1];
              if (index >= 3) {
                  ntarget4 = nfunc[3];
                  nsource3 = nfunc[2];
                  nk = &nfunck;
                  ogc3 = &ogcfunc[2];
                }
            }
        }
    }

  switch (index) {
  case 0:
      shell = sh1;
      tgencon = &i;
      break;
  case 1:
      shell = sh2;
      tgencon = &j;
      break;
  case 2:
      shell = sh3;
      tgencon = &k;
      break;
  case 3:
      shell = sh4;
      tgencon = &l;
      break;
    }

#if PRINT
    {
      double *tmp = integrals;
      printf("Before transform of index %d (%dx%dx%dx%d)\n",
             index, nsource1, nsource2, nsource3, nsource4);
      for (i=0; i<nsource1; i++) {
          for (j=0; j<nsource2; j++) {
              for (k=0; k<nsource3; k++) {
                  for (l=0; l<nsource4; l++) {
                      if (fabs(*tmp)>1.e-15) {
                          printf("(%d %d|%d %d) = %15.11lf\n",i,j,k,l,*tmp);
                        }
                      tmp++;
                    }
                }
            }
        }
    }
#endif

  copy_to_source(integrals, nsource1*nsource2*nsource3*nsource4);
  memset(target, 0, sizeof(double)*ntarget1*ntarget2*ntarget3*ntarget4);

  ogccart[0] = 0;
  ogcfunc[0] = 0;
  for (i=0; i<sh1->ncon; i++) {
      am1 = sh1->type[i].am;
      nfunci = INT_NFUNC(sh1->type[i].puream,sh1->type[i].am);
      ncarti = INT_NCART(am1);
      ogccart[1] = 0;
      ogcfunc[1] = 0;
      for (j=0; j<sh2->ncon; j++) {
          am2 = sh2->type[j].am;
          nfuncj = INT_NFUNC(sh2->type[j].puream,sh2->type[j].am);
          ncartj = INT_NCART(am2);
          ogccart[2] = 0;
          ogcfunc[2] = 0;
          for (k=0; k<sh3->ncon; k++) {
              am3 = sh3->type[k].am;
              nfunck = INT_NFUNC(sh3->type[k].puream,sh3->type[k].am);
              ncartk = INT_NCART(am3);
              ogccart[3] = 0;
              ogcfunc[3] = 0;
              for (l=0; l<sh4->ncon; l++) {
                  am4 = sh4->type[l].am;
                  nfuncl = INT_NFUNC(sh4->type[l].puream,sh4->type[l].am);
                  ncartl = INT_NCART(am4);

                  if (shell->type[*tgencon].puream) {
                      SphericalTransformIterV2
                          trans(shell->type[*tgencon].am);
                      do_sparse_transform2(source, target,
                                           index, trans,
                                           ncart[index], nfunc[index],
                                           ogccart[index], ogcfunc[index],
                                           *ni, nsource1, *ogc1,
                                           *nj, nsource2, *ogc2,
                                           *nk, nsource3, *ogc3,
                                           *nl, nsource4, *ogc4);
                    }
                  else {
                      do_copy2(source, integrals,
                               *ni, nsource1, *ogc1,
                               *nj, nsource2, *ogc2,
                               *nk, nsource3, *ogc3,
                               *nl, nsource4, *ogc4);
                    }
                  ogccart[3] += ncartl;
                  ogcfunc[3] += nfuncl;
                }
              ogccart[2] += ncartk;
              ogcfunc[2] += nfunck;
            }
          ogccart[1] += ncartj;
          ogcfunc[1] += nfuncj;
        }
      ogccart[0] += ncarti;
      ogcfunc[0] += nfunci;
    }

#if PRINT
    {
      double *tmp = integrals;
      printf("After transform of index %d (%dx%dx%dx%d)\n",
             index, ntarget1, ntarget2, ntarget3, ntarget4);
      for (i=0; i<ntarget1; i++) {
          for (j=0; j<ntarget2; j++) {
              for (k=0; k<ntarget3; k++) {
                  for (l=0; l<ntarget4; l++) {
                      if (fabs(*tmp)>1.e-15) {
                          printf("(%d %d|%d %d) = %15.11lf\n",i,j,k,l,*tmp);
                        }
                      tmp++;
                    }
                }
            }
        }
    }
#endif
}

void
int_transform_2e(double *integrals, double *target,
                 shell_t *sh1, shell_t *sh2,
                 shell_t *sh3, shell_t *sh4)
{
  int pure1 = have_pure(sh1);
  int pure2 = have_pure(sh2);
  int pure3 = have_pure(sh3);
  int pure4 = have_pure(sh4);

  if (pure1) {
      do_gencon_sparse_transform_2e(integrals, target, 0, sh1, sh2, sh3, sh4);
      integrals = target;
    }
  if (pure2) {
      do_gencon_sparse_transform_2e(integrals, target, 1, sh1, sh2, sh3, sh4);
      integrals = target;
    }
  if (pure3) {
      do_gencon_sparse_transform_2e(integrals, target, 2, sh1, sh2, sh3, sh4);
      integrals = target;
    }
  if (pure4) {
      do_gencon_sparse_transform_2e(integrals, target, 3, sh1, sh2, sh3, sh4);
      integrals = target;
    }

  if (integrals != target) {
      int nint = sh1->nfunc * sh2->nfunc * sh3->nfunc * sh4->nfunc;
      memmove(target, integrals, sizeof(double)*nint);
    }
}

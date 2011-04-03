//
// tformv3.cc
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@limitpt.com>
// Maintainer: LPS
//
// This file is part of the SC Toolkit.
//
// The SC Toolkit is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The SC Toolkit is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the SC Toolkit; see the file COPYING.LIB.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//
// The U.S. Government is granted a limited license as per AL 91-7.
//

#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <util/misc/formio.h>
#include <chemistry/qc/basis/integral.h>
#include <chemistry/qc/intv3/macros.h>
#include <chemistry/qc/intv3/tformv3.h>
#include <chemistry/qc/intv3/utils.h>

using namespace sc;

////////////////////////////////////////////////////////////////////////////

#define PRINT 0

void
Int2eV3::transform_init()
{
  source = 0;
  nsourcemax = 0;
}

void
Int2eV3::transform_done()
{
  delete[] source;
}

void
Int1eV3::transform_init()
{
  source = 0;
  nsourcemax = 0;
}

void
Int1eV3::transform_done()
{
  delete[] source;
}

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

static void
do_sparse_transform2_1(double *source, double *target,
                       SphericalTransformIter& trans,
                       int stcart, int stpure,
                       int ogctcart, int ogctpure,
                       int n2, int s2, int ogc2,
                       int n3, int s3, int ogc3,
                       int n4, int s4, int ogc4)
{
  int i2, i3, i4;

  for (trans.begin(); trans.ready(); trans.next()) {
      double coef = trans.coef();
      int pure = trans.pureindex();
      int cart = trans.cartindex();
      int offtarget1 = ogctpure + pure;
      int offsource1 = ogctcart + cart;
      int offtarget2 = offtarget1*s2 + ogc2;
      int offsource2 = offsource1*s2 + ogc2;
      for (i2=0; i2<n2; i2++,offtarget2++,offsource2++) {
          int offtarget3 = offtarget2*s3 + ogc3;
          int offsource3 = offsource2*s3 + ogc3;
          for (i3=0; i3<n3; i3++,offtarget3++,offsource3++) {
              int offtarget4 = offtarget3*s4 + ogc4;
              int offsource4 = offsource3*s4 + ogc4;
              for (i4=0; i4<n4; i4++,offtarget4++,offsource4++) {
                  target[offtarget4] += coef * source[offsource4];
                }
            }
        }
    }
}

static void
do_sparse_transform2_2(double *source, double *target,
                       SphericalTransformIter& trans,
                       int stcart, int stpure,
                       int ogctcart, int ogctpure,
                       int n1, int s1, int ogc1,
                       int n3, int s3, int ogc3,
                       int n4, int s4, int ogc4)
{
  int i1, i3, i4;

  for (trans.begin(); trans.ready(); trans.next()) {
      double coef = trans.coef();
      int pure = trans.pureindex();
      int cart = trans.cartindex();
      int offtarget1 = ogc1;
      int offsource1 = ogc1;
      for (i1=0; i1<n1; i1++,offtarget1++,offsource1++) {
          int offtarget2 = offtarget1*stpure + ogctpure + pure;
          int offsource2 = offsource1*stcart + ogctcart + cart;
          int offtarget3 = offtarget2*s3 + ogc3;
          int offsource3 = offsource2*s3 + ogc3;
          for (i3=0; i3<n3; i3++,offtarget3++,offsource3++) {
              int offtarget4 = offtarget3*s4 + ogc4;
              int offsource4 = offsource3*s4 + ogc4;
              for (i4=0; i4<n4; i4++,offtarget4++,offsource4++) {
                  target[offtarget4] += coef * source[offsource4];
                }
            }
        }
    }
}

static void
do_sparse_transform2_3(double *source, double *target,
                       SphericalTransformIter& trans,
                       int stcart, int stpure,
                       int ogctcart, int ogctpure,
                       int n1, int s1, int ogc1,
                       int n2, int s2, int ogc2,
                       int n4, int s4, int ogc4)
{
  int i1, i2, i4;

  for (trans.begin(); trans.ready(); trans.next()) {
      double coef = trans.coef();
      int pure = trans.pureindex();
      int cart = trans.cartindex();
      int offtarget1 = ogc1;
      int offsource1 = ogc1;
      for (i1=0; i1<n1; i1++,offtarget1++,offsource1++) {
          int offtarget2 = offtarget1*s2 + ogc2;
          int offsource2 = offsource1*s2 + ogc2;
          for (i2=0; i2<n2; i2++,offtarget2++,offsource2++) {
              int offtarget3 = offtarget2*stpure + ogctpure + pure;
              int offsource3 = offsource2*stcart + ogctcart + cart;
              int offtarget4 = offtarget3*s4 + ogc4;
              int offsource4 = offsource3*s4 + ogc4;
              for (i4=0; i4<n4; i4++,offtarget4++,offsource4++) {
                  target[offtarget4] += coef * source[offsource4];
                }
            }
        }
    }
}

static void
do_sparse_transform2_4(double *source, double *target,
                       SphericalTransformIter& trans,
                       int stcart, int stpure,
                       int ogctcart, int ogctpure,
                       int n1, int s1, int ogc1,
                       int n2, int s2, int ogc2,
                       int n3, int s3, int ogc3)
{
  int i1, i2, i3;

  for (trans.begin(); trans.ready(); trans.next()) {
      double coef = trans.coef();
      int pure = trans.pureindex();
      int cart = trans.cartindex();
      int offtarget1 = ogc1;
      int offsource1 = ogc1;
      for (i1=0; i1<n1; i1++,offtarget1++,offsource1++) {
          int offtarget2 = offtarget1*s2 + ogc2;
          int offsource2 = offsource1*s2 + ogc2;
          for (i2=0; i2<n2; i2++,offtarget2++,offsource2++) {
              int offtarget3 = offtarget2*s3 + ogc3;
              int offsource3 = offsource2*s3 + ogc3;
              int offtarget4 = offtarget3*stpure + ogctpure + pure;
              int offsource4 = offsource3*stcart + ogctcart + cart;
              for (i3=0; i3<n3; i3++,offtarget4+=stpure,offsource4+=stcart) {
                  //for (i3=0; i3<n3; i3++,offtarget3++,offsource3++) {
                  //int offtarget4 = offtarget3*stpure + ogctpure + pure;
                  //int offsource4 = offsource3*stcart + ogctcart + cart;
                  target[offtarget4] += coef * source[offsource4];
                }
            }
        }
    }
}

static void
do_sparse_transform2(double *source, double *target,
                     int index, SphericalTransformIter& trans,
                     int stcart, int stpure,
                     int ogctcart, int ogctpure,
                     int n1, int s1, int ogc1,
                     int n2, int s2, int ogc2,
                     int n3, int s3, int ogc3,
                     int n4, int s4, int ogc4)
{
  switch (index) {
  case 0:
      do_sparse_transform2_1(source, target, trans,
                             stcart,stpure,
                             ogctcart, ogctpure,
                             n2, s2, ogc2,
                             n3, s3, ogc3,
                             n4, s4, ogc4);
      break;
  case 1:
      do_sparse_transform2_2(source, target, trans,
                             stcart,stpure,
                             ogctcart, ogctpure,
                             n1, s1, ogc1,
                             n3, s3, ogc3,
                             n4, s4, ogc4);
      break;
  case 2:
      do_sparse_transform2_3(source, target, trans,
                             stcart,stpure,
                             ogctcart, ogctpure,
                             n1, s1, ogc1,
                             n2, s2, ogc2,
                             n4, s4, ogc4);
      break;
  case 3:
      do_sparse_transform2_4(source, target, trans,
                             stcart,stpure,
                             ogctcart, ogctpure,
                             n1, s1, ogc1,
                             n2, s2, ogc2,
                             n3, s3, ogc3);
      break;
    }
}

/* make sure enough space exists for the source integrals */
void
Int2eV3::source_space(int nsource)
{
  if (nsourcemax < nsource) {
      delete[] source;
      source = new double[nsource*3];
      target = &source[nsource];
      scratch = &source[nsource*2];
      nsourcemax = nsource;
    }
}

void
Int2eV3::copy_to_source(double *integrals, int nsource)
{
  int i;
  double *tmp, *tmp2;

  /* Allocate more temporary space if needed. */
  source_space(nsource);

  tmp = source;
  tmp2 = integrals;
  for (i=0; i<nsource; i++) *tmp++ = *tmp2++;
}

/* make sure enough space exists for the source integrals */
void
Int1eV3::source_space(int nsource)
{
  if (nsourcemax < nsource) {
      delete[] source;
      source = new double[nsource];
      nsourcemax = nsource;
    }
}

void
Int1eV3::copy_to_source(double *integrals, int nsource)
{
  int i;
  double *tmp, *tmp2;

  /* Allocate more temporary space if needed. */
  source_space(nsource);

  tmp = source;
  tmp2 = integrals;
  for (i=0; i<nsource; i++) *tmp++ = *tmp2++;
}

void
Int1eV3::do_transform_1e(Integral *integ,
                         double *integrals,
                         GaussianShell *sh1, GaussianShell *sh2,
                         int chunk)
{
  int i, j;
  int ogc1, ogc2;
  int ogc1pure, ogc2pure;
  int am1, am2;
  int pure1 = sh1->has_pure();
  int pure2 = sh2->has_pure();
  int ncart1 = sh1->ncartesian();
  int ncart2 = sh2->ncartesian();
  int nfunc1 = sh1->nfunction();
  int nfunc2 = sh2->nfunction();
  int nfunci, nfuncj;

  if (!pure1 && !pure2) return;

  /* Loop through the generalized general contractions,
   * transforming the first index. */
  if (pure1) {
      copy_to_source(integrals, ncart1*ncart2*chunk);
      memset(integrals, 0, sizeof(double)*sh1->nfunction()*ncart2*chunk);

      ogc1 = 0;
      ogc1pure = 0;
      for (i=0; i<sh1->ncontraction(); i++) {
          am1 = sh1->am(i);
          nfunci = sh1->nfunction(i);
          ogc2 = 0;
          for (j=0; j<sh2->ncontraction(); j++) {
              am2 = sh2->am(j);
              nfuncj = sh2->nfunction(j);

              if (sh1->is_pure(i)) {
                  SphericalTransformIter
                      trans(integ->spherical_transform(sh1->am(i)));
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
      memset(integrals, 0,
             sizeof(double)*sh1->nfunction()*sh2->nfunction()*chunk);

      ogc1 = 0;
      for (i=0; i<sh1->ncontraction(); i++) {
          am1 = sh1->am(i);
          nfunci = sh1->nfunction(i);
          ogc2 = 0;
          ogc2pure = 0;
          for (j=0; j<sh2->ncontraction(); j++) {
              am2 = sh2->am(j);
              nfuncj = sh2->nfunction(j);

              if (sh2->is_pure(j)) {
                  SphericalTransformIter
                      trans(integ->spherical_transform(sh2->am(j)));
                  do_sparse_transform12(source, integrals, chunk,
                                        trans,
                                        INT_NPURE(am1), ogc1,
                                        ncart2, ogc2,
                                        sh2->nfunction(), ogc2pure);
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
Int1eV3::transform_1e(Integral *integ,
                      double *integrals, double *target,
                      GaussianShell *sh1, GaussianShell *sh2, int chunk)
{
  int ntarget;

  do_transform_1e(integ, integrals, sh1, sh2, chunk);

  /* copy the integrals to the target, if necessary */
  ntarget = sh1->nfunction() * sh2->nfunction();
  if (integrals != target) {
      memmove(target, integrals, ntarget*sizeof(double)*chunk);
    }
}

/* it is not ok for integrals and target to overlap */
void
Int1eV3::accum_transform_1e(Integral *integ,
                            double *integrals, double *target,
                            GaussianShell *sh1, GaussianShell *sh2, int chunk)
{
  int i, ntarget;

  do_transform_1e(integ, integrals, sh1, sh2, chunk);

  /* accum the integrals to the target */
  ntarget = sh1->nfunction() * sh2->nfunction() * chunk;
  for (i=0; i<ntarget; i++) target[i] += integrals[i];
}

void
Int1eV3::transform_1e(Integral*integ,
                            double *integrals, double *target,
                            GaussianShell *sh1, GaussianShell *sh2)
{
  transform_1e(integ, integrals, target, sh1, sh2, 1);
}

void
Int1eV3::accum_transform_1e(Integral*integ,
                                  double *integrals, double *target,
                                  GaussianShell *sh1, GaussianShell *sh2)
{
  accum_transform_1e(integ, integrals, target, sh1, sh2, 1);
}

void
Int1eV3::transform_1e_xyz(Integral*integ,
                          double *integrals, double *target,
                          GaussianShell *sh1, GaussianShell *sh2)
{
  transform_1e(integ, integrals, target, sh1, sh2, 3);
}

void
Int1eV3::accum_transform_1e_xyz(Integral*integ,
                                double *integrals, double *target,
                                GaussianShell *sh1, GaussianShell *sh2)
{
  accum_transform_1e(integ, integrals, target, sh1, sh2, 3);
}

void
Int2eV3::do_gencon_sparse_transform_2e(Integral*integ,
                                       double *integrals, double *target,
                                       int index,
                                       GaussianShell *sh1, GaussianShell *sh2,
                                       GaussianShell *sh3, GaussianShell *sh4)
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

  GaussianShell *shell;

  int *tgencon;

  ncart[0] = sh1->ncartesian();
  ncart[1] = sh2->ncartesian();
  ncart[2] = sh3->ncartesian();
  ncart[3] = sh4->ncartesian();
  nfunc[0] = sh1->nfunction();
  nfunc[1] = sh2->nfunction();
  nfunc[2] = sh3->nfunction();
  nfunc[3] = sh4->nfunction();

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
  default:
      shell = 0;
      tgencon = 0;
      break;
    }

#if PRINT
    {
      double *tmp = integrals;
      ExEnv::outn() << scprintf("Before transform of index %d (%dx%dx%dx%d)\n",
             index, nsource1, nsource2, nsource3, nsource4);
      for (i=0; i<nsource1; i++) {
          for (j=0; j<nsource2; j++) {
              for (k=0; k<nsource3; k++) {
                  for (l=0; l<nsource4; l++) {
                      if (fabs(*tmp)>1.e-15) {
                          ExEnv::outn() << scprintf("(%d %d|%d %d) = %15.11lf\n",i,j,k,l,*tmp);
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
  for (i=0; i<sh1->ncontraction(); i++) {
      am1 = sh1->am(i);
      nfunci = sh1->nfunction(i);
      ncarti = INT_NCART(am1);
      ogccart[1] = 0;
      ogcfunc[1] = 0;
      for (j=0; j<sh2->ncontraction(); j++) {
          am2 = sh2->am(j);
          nfuncj = sh2->nfunction(j);
          ncartj = INT_NCART(am2);
          ogccart[2] = 0;
          ogcfunc[2] = 0;
          for (k=0; k<sh3->ncontraction(); k++) {
              am3 = sh3->am(k);
              nfunck = sh3->nfunction(k);
              ncartk = INT_NCART(am3);
              ogccart[3] = 0;
              ogcfunc[3] = 0;
              for (l=0; l<sh4->ncontraction(); l++) {
                  am4 = sh4->am(l);
                  nfuncl = sh4->nfunction(l);
                  ncartl = INT_NCART(am4);

                  if (shell->is_pure(*tgencon)) {
                      SphericalTransformIter
                        trans(integ->spherical_transform(shell->am(*tgencon)));
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
      ExEnv::outn() << scprintf("After transform of index %d (%dx%dx%dx%d)\n",
             index, ntarget1, ntarget2, ntarget3, ntarget4);
      for (i=0; i<ntarget1; i++) {
          for (j=0; j<ntarget2; j++) {
              for (k=0; k<ntarget3; k++) {
                  for (l=0; l<ntarget4; l++) {
                      if (fabs(*tmp)>1.e-15) {
                          ExEnv::outn()
                            << scprintf("(%d %d|%d %d) = %15.11lf\n",
                                        i,j,k,l,*tmp);
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
Int2eV3::transform_2e_slow(Integral *integ, double *integrals, double *target,
                           GaussianShell *sh1, GaussianShell *sh2,
                           GaussianShell *sh3, GaussianShell *sh4)
{
  int pure1 = sh1->has_pure();
  int pure2 = sh2->has_pure();
  int pure3 = sh3->has_pure();
  int pure4 = sh4->has_pure();

  if (pure1) {
      do_gencon_sparse_transform_2e(integ,
                                    integrals, target, 0, sh1, sh2, sh3, sh4);
      integrals = target;
    }
  if (pure2) {
      do_gencon_sparse_transform_2e(integ,
                                    integrals, target, 1, sh1, sh2, sh3, sh4);
      integrals = target;
    }
  if (pure3) {
      do_gencon_sparse_transform_2e(integ,
                                    integrals, target, 2, sh1, sh2, sh3, sh4);
      integrals = target;
    }
  if (pure4) {
      do_gencon_sparse_transform_2e(integ,
                                    integrals, target, 3, sh1, sh2, sh3, sh4);
      integrals = target;
    }

  if (integrals != target) {
      int nint = sh1->nfunction() * sh2->nfunction()
               * sh3->nfunction() * sh4->nfunction();
      memmove(target, integrals, sizeof(double)*nint);
    }
}

/////////////////////////////////////////////////////////////////////////////

static void
do_sparse_transform2_1new(double *source, double *target,
                          SphericalTransformIter& trans,
                          int stcart, int stpure,
                          int n2,
                          int n3,
                          int n4)
{
  int i234, n234=n2*n3*n4;

  for (trans.begin(); trans.ready(); trans.next()) {
    double coef = trans.coef();
    int pure = trans.pureindex();
    int cart = trans.cartindex();
    int offtarget4 = pure*n234;
    int offsource4 = cart*n234;
    for (i234=0; i234<n234; i234++,offtarget4++,offsource4++) {
      target[offtarget4] += coef * source[offsource4];
      }
    }
}

static void
do_sparse_transform2_2new(double *source, double *target,
                          SphericalTransformIter& trans,
                          int stcart, int stpure,
                          int n1,
                          int n3,
                          int n4)
{
  int i1, i34, n34=n3*n4;
  int n34stpure=n34*(stpure-1); // -1 because of the increment int the loop
  int n34stcart=n34*(stcart-1); // ditto

  for (trans.begin(); trans.ready(); trans.next()) {
    double coef = trans.coef();
    int pure = trans.pureindex();
    int cart = trans.cartindex();
    int offtarget4 = pure*n34;
    int offsource4 = cart*n34;
    for (i1=0; i1<n1; i1++) {
      for (i34=0; i34<n34; i34++,offtarget4++,offsource4++) {
        target[offtarget4] += coef * source[offsource4];
        }
      offtarget4 += n34stpure;
      offsource4 += n34stcart;
      }
    }
}

static void
do_sparse_transform2_3new(double *source, double *target,
                          SphericalTransformIter& trans,
                          int stcart, int stpure,
                          int n1,
                          int n2,
                          int n4)
{
  int i12, i4, n12=n1*n2, n4stpure=n4*stpure, n4stcart=n4*stcart;

  for (trans.begin(); trans.ready(); trans.next()) {
    double coef = trans.coef();
    int pure = trans.pureindex();
    int cart = trans.cartindex();
    int offtarget3 = pure*n4;
    int offsource3 = cart*n4;
    for (i12=0; i12<n12; i12++) {
      int offtarget4 = offtarget3;
      int offsource4 = offsource3;
      for (i4=0; i4<n4; i4++,offtarget4++,offsource4++) {
        target[offtarget4] += coef * source[offsource4];
        }
      offtarget3 += n4stpure;
      offsource3 += n4stcart;
      }
    }
}

static void
do_sparse_transform2_4new(double *source, double *target,
                          SphericalTransformIter& trans,
                          int stcart, int stpure,
                          int n1,
                          int n2,
                          int n3)
{
  int n123=n1*n2*n3;

  for (trans.begin(); trans.ready(); trans.next()) {
    double coef = trans.coef();
    int pure = trans.pureindex();
    int cart = trans.cartindex();
    int offtarget4 = pure;
    int offsource4 = cart;
    for (int i123=0; i123<n123; i123++) {
      target[offtarget4] += coef * source[offsource4];
      offtarget4 += stpure;
      offsource4 += stcart;
      }
    }
}

// Cartint and pureint may overlap.  The must be enough space
// in pureint to hold all the cartints.  The cartint buffer
// will be overwritten in any case.
void
Int2eV3::transform_2e(Integral *integ,
                      double *cartint, double *pureint,
                      GaussianShell *sh1, GaussianShell *sh2,
                      GaussianShell *sh3, GaussianShell *sh4)
{
  int pure1 = sh1->has_pure();
  int pure2 = sh2->has_pure();
  int pure3 = sh3->has_pure();
  int pure4 = sh4->has_pure();

  int nfunc1=sh1->nfunction();
  int nfunc2=sh2->nfunction();
  int nfunc3=sh3->nfunction();
  int nfunc4=sh4->nfunction();
  int nfunc34 =  nfunc3 * nfunc4;
  int nfunc234 = nfunc2 * nfunc34;
  int nfunc1234 = nfunc1 * nfunc234;

  if (!pure1 && !pure2 && !pure3 && !pure4) {
    if (pureint!=cartint) memmove(pureint, cartint, sizeof(double)*nfunc1234);
    return;
    }

  int ncart1=sh1->ncartesian();
  int ncart2=sh2->ncartesian();
  int ncart3=sh3->ncartesian();
  int ncart4=sh4->ncartesian();
  int ncart34 =  ncart3 * ncart4;
  int ncart234 = ncart2 * ncart34;

  // allocate the scratch arrays, if needed
  source_space(ncart1*ncart234);

  int ncon1 = sh1->ncontraction();
  int ncon2 = sh2->ncontraction();
  int ncon3 = sh3->ncontraction();
  int ncon4 = sh4->ncontraction();

  if (ncon1==1 && ncon2==1 && ncon3==1 && ncon4==1) {
    double *sourcebuf = cartint;
    double *targetbuf = target;
    // transform indices
    if (pure1) {
      SphericalTransformIter transi(integ->spherical_transform(sh1->am(0)));
      memset(targetbuf,0,sizeof(double)*nfunc1*ncart2*ncart3*ncart4);
      do_sparse_transform2_1new(sourcebuf, targetbuf, transi,
                                ncart1, nfunc1,
                                ncart2,
                                ncart3,
                                ncart4);
      double*tmp=sourcebuf; sourcebuf=targetbuf; targetbuf=tmp;
      }
    if (pure2) {
      SphericalTransformIter transj(integ->spherical_transform(sh2->am(0)));
      memset(targetbuf,0,sizeof(double)*nfunc1*nfunc2*ncart3*ncart4);
      do_sparse_transform2_2new(sourcebuf, targetbuf, transj,
                                ncart2, nfunc2,
                                nfunc1,
                                ncart3,
                                ncart4);
      double*tmp=sourcebuf; sourcebuf=targetbuf; targetbuf=tmp;
      }
    if (pure3) {
      SphericalTransformIter transk(integ->spherical_transform(sh3->am(0)));
      memset(targetbuf,0,sizeof(double)*nfunc1*nfunc2*nfunc3*ncart4);
      do_sparse_transform2_3new(sourcebuf, targetbuf, transk,
                                ncart3, nfunc3,
                                nfunc1,
                                nfunc2,
                                ncart4);
      double*tmp=sourcebuf; sourcebuf=targetbuf; targetbuf=tmp;
      }
    if (pure4) {
      SphericalTransformIter transl(integ->spherical_transform(sh4->am(0)));
      memset(targetbuf,0,sizeof(double)*nfunc1234);
      do_sparse_transform2_4new(sourcebuf, targetbuf, transl,
                                ncart4, nfunc4,
                                nfunc1,
                                nfunc2,
                                nfunc3);
      double*tmp=sourcebuf; sourcebuf=targetbuf; targetbuf=tmp;
      }
    if (sourcebuf!=pureint)
      memmove(pureint, sourcebuf, sizeof(double)*nfunc1234);
    }
  else {
    // begin gc loop
    int ogccart1 = 0;
    int ogcfunc1 = 0;
    for (int i=0; i<ncon1; i++) {
      int am1 = sh1->am(i);
      int nfunci = sh1->nfunction(i);
      int ispurei = sh1->is_pure(i);
      int ncarti = INT_NCART_NN(am1);
      int ogccart2 = 0;
      int ogcfunc2 = 0;
      SphericalTransformIter transi(integ->spherical_transform(am1));
      for (int j=0; j<ncon2; j++) {
        int am2 = sh2->am(j);
        int nfuncj = sh2->nfunction(j);
        int ispurej = sh2->is_pure(j);
        int ncartj = INT_NCART_NN(am2);
        int ogccart3 = 0;
        int ogcfunc3 = 0;
        SphericalTransformIter transj(integ->spherical_transform(am2));
        for (int k=0; k<ncon3; k++) {
          int am3 = sh3->am(k);
          int nfunck = sh3->nfunction(k);
          int ispurek = sh3->is_pure(k);
          int ncartk = INT_NCART_NN(am3);
          int ogccart4 = 0;
          int ogcfunc4 = 0;
          SphericalTransformIter transk(integ->spherical_transform(am3));
          for (int l=0; l<ncon4; l++) {
            int am4 = sh4->am(l);
            int nfuncl = sh4->nfunction(l);
            int ispurel = sh4->is_pure(l);
            int ncartl = INT_NCART_NN(am4);
    
    ;
    // copy to source buffer
    int cartindex1 = ogccart1*ncart234
                   + ogccart2*ncart34 + ogccart3*ncart4 + ogccart4;
    double *tmp_source = source;
    int is;
    for (is=0; is<ncarti; is++) {
      int cartindex12 = cartindex1;
      for (int js=0; js<ncartj; js++) {
        int cartindex123 = cartindex12;
        for (int ks=0; ks<ncartk; ks++) {
          double *tmp_cartint = &cartint[cartindex123];
          for (int ls=0; ls<ncartl; ls++) {
            *tmp_source++ = *tmp_cartint++;
            }
          cartindex123 += ncart4;
          }
        cartindex12 += ncart34;
        }
      cartindex1 += ncart234;
      }
    
    
    double *sourcebuf = source;
    double *targetbuf = target;
    
    // transform indices
    if (ispurei) {
      memset(targetbuf,0,sizeof(double)*nfunci*ncartj*ncartk*ncartl);
      do_sparse_transform2_1new(sourcebuf, targetbuf, transi,
                                ncarti, nfunci,
                                ncartj,
                                ncartk,
                                ncartl);
      double*tmp=sourcebuf; sourcebuf=targetbuf; targetbuf=tmp;
      }
    if (ispurej) {
      memset(targetbuf,0,sizeof(double)*nfunci*nfuncj*ncartk*ncartl);
      do_sparse_transform2_2new(sourcebuf, targetbuf, transj,
                                ncartj, nfuncj,
                                nfunci,
                                ncartk,
                                ncartl);
      double*tmp=sourcebuf; sourcebuf=targetbuf; targetbuf=tmp;
      }
    if (ispurek) {
      memset(targetbuf,0,sizeof(double)*nfunci*nfuncj*nfunck*ncartl);
      do_sparse_transform2_3new(sourcebuf, targetbuf, transk,
                                ncartk, nfunck,
                                nfunci,
                                nfuncj,
                                ncartl);
      double*tmp=sourcebuf; sourcebuf=targetbuf; targetbuf=tmp;
      }
    if (ispurel) {
      memset(targetbuf,0,sizeof(double)*nfunci*nfuncj*nfunck*nfuncl);
      SphericalTransformIter transl(integ->spherical_transform(am4));
      do_sparse_transform2_4new(sourcebuf, targetbuf, transl,
                                ncartl, nfuncl,
                                nfunci,
                                nfuncj,
                                nfunck);
      double*tmp=sourcebuf; sourcebuf=targetbuf; targetbuf=tmp;
      }
    
    // copy to scratch buffer
    int funcindex1 = ogcfunc1*nfunc234
                   + ogcfunc2*nfunc34 + ogcfunc3*nfunc4 + ogcfunc4;
    tmp_source = sourcebuf;
    for (is=0; is<nfunci; is++) {
      int funcindex12 = funcindex1;
      for (int js=0; js<nfuncj; js++) {
        int funcindex123 = funcindex12;
        for (int ks=0; ks<nfunck; ks++) {
          double *tmp_scratch = &scratch[funcindex123];
          for (int ls=0; ls<nfuncl; ls++) {
            *tmp_scratch++ = *tmp_source++;
            }
          funcindex123 += nfunc4;
          }
        funcindex12 += nfunc34;
        }
      funcindex1 += nfunc234;
      }
    
    // end gc loop
            ogccart4 += ncartl;
            ogcfunc4 += nfuncl;
            }
          ogccart3 += ncartk;
          ogcfunc3 += nfunck;
          }
        ogccart2 += ncartj;
        ogcfunc2 += nfuncj;
        }
      ogccart1 += ncarti;
      ogcfunc1 += nfunci;
      }
    memcpy(pureint, scratch, sizeof(double)*nfunc1234);
    }
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:

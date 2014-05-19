//
// tformv3.cc
//
// Copyright (C) 2001 Edward Valeev
//
// Author: Edward Valeev <evaleev@vt.edu>
// Maintainer: EV
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
#include <chemistry/qc/libint2/macros.h>
#include <chemistry/qc/libint2/tform.h>
#include <chemistry/qc/libint2/int1e.h>
#include <chemistry/qc/libint2/int2e.h>
#include <chemistry/qc/libint2/gto.h>
#include <chemistry/qc/libint2/tform.timpl.h>

using namespace std;
using namespace sc;

inline int max(int a,int b) { return (a > b) ? a : b;}
inline void fail()
{
  ExEnv::errn() << scprintf("failing module:\n%s",__FILE__) << endl;
  abort();
}

static void transform1e_1(SphericalTransformIter&, double*, double*, int);
static void transform1e_2(SphericalTransformIter&, double*, double*, int, int);
static void transform1e_vec_2(const int, SphericalTransformIter&, double*, double*, int, int);
static void transform2e_1(SphericalTransformIter&, double*, double*, int);
static void transform2e_2(SphericalTransformIter&, double*, double*, int, int, int);
static void transform2e_3(SphericalTransformIter&, double*, double*, int, int, int);
static void transform2e_4(SphericalTransformIter&, double*, double*, int, int);

void Int1eLibint2::transform_contrquartets_(double * source_ints_buf, double *target_ints_buf)
{
  double *source1, *target1;
  double *source2, *target2;
  double *source = source_ints_buf;
  double *target = target_ints_buf;
  double *tmpbuf = tformbuf_;

  for (int gc1=0; gc1<int_shell1_->ncontraction(); gc1++) {
    int am1 = int_shell1_->am(gc1);
    int is_pure1 = int_shell1_->is_pure(gc1) ? 1 : 0;
    int ncart1 = int_shell1_->ncartesian(gc1);
    int nbf1 = int_shell1_->nfunction(gc1);

    int target_bf2_offset = 0;
    for (int gc2=0; gc2<int_shell2_->ncontraction(); gc2++) {
      int am2 = int_shell2_->am(gc2);
      int is_pure2 = int_shell2_->is_pure(gc2) ? 1 : 0;
      int ncart2 = int_shell2_->ncartesian(gc2);
      int nbf2 = int_shell2_->nfunction(gc2);
      
      int transform_index = 2*is_pure1 + is_pure2;
      switch (transform_index) {
      case 0:
	break;

      case 1:
	source2 = source;
	target2 = target;
	break;
	
      case 2:
	source1 = source;
	target1 = target;
	break;

      case 3:
	source2 = source;
	target2 = tmpbuf;
	source1 = tmpbuf;
	target1 = target;
	break;
      }

      if (is_pure2) {
	SphericalTransformIter stiter(integral_->spherical_transform(am2));
	transform1e_2(stiter,source2, target2, ncart1,ncart2);
      }
      if (is_pure1) {
	SphericalTransformIter stiter(integral_->spherical_transform(am1));
	transform1e_1(stiter,source1, target1, nbf2);
      }
      
      source += (ncart1*ncart2);
      target += (nbf1*nbf2);
    }
  }
}

void Int1eLibint2::transform_contrquartets_vec_(const int ntypes, double * source_ints_buf, double *target_ints_buf)
{
  double *source1, *target1;
  double *source2, *target2;
  double *source = source_ints_buf;
  double *target = target_ints_buf;
  double *tmpbuf = tformbuf_;

  for (int gc1=0; gc1<int_shell1_->ncontraction(); gc1++) {
    int am1 = int_shell1_->am(gc1);
    int is_pure1 = int_shell1_->is_pure(gc1) ? 1 : 0;
    int ncart1 = int_shell1_->ncartesian(gc1);
    int nbf1 = int_shell1_->nfunction(gc1);

    int target_bf2_offset = 0;
    for (int gc2=0; gc2<int_shell2_->ncontraction(); gc2++) {
      int am2 = int_shell2_->am(gc2);
      int is_pure2 = int_shell2_->is_pure(gc2) ? 1 : 0;
      int ncart2 = int_shell2_->ncartesian(gc2);
      int nbf2 = int_shell2_->nfunction(gc2);
      
      int transform_index = 2*is_pure1 + is_pure2;
      switch (transform_index) {
      case 0:
	break;

      case 1:
	source2 = source;
	target2 = target;
	break;
	
      case 2:
	source1 = source;
	target1 = target;
	break;

      case 3:
	source2 = source;
	target2 = tmpbuf;
	source1 = tmpbuf;
	target1 = target;
	break;
      }

      if (is_pure2) {
	SphericalTransformIter stiter(integral_->spherical_transform(am2));
	transform1e_vec_2(ntypes, stiter, source2, target2, ncart1,ncart2);
      }
      if (is_pure1) {
	SphericalTransformIter stiter(integral_->spherical_transform(am1));
	transform1e_1(stiter, source1, target1, nbf2*ntypes);
      }
      
      source += (ntypes*ncart1*ncart2);
      target += (ntypes*nbf1*nbf2);
    }
  }
}

void Int2eLibint2::transform_contrquartets_(double * source_ints_buf, double *target_ints_buf)
{
  double *source1, *target1;
  double *source2, *target2;
  double *source3, *target3;
  double *source4, *target4;
  double *source = source_ints_buf;
  double *target = target_ints_buf;
  double *tmpbuf = &tformbuf_[0];

  for (int gc1=0; gc1<int_shell1_->ncontraction(); gc1++) {
    int am1 = int_shell1_->am(gc1);
    int is_pure1 = int_shell1_->is_pure(gc1) ? 1 : 0;
    int ncart1 = int_shell1_->ncartesian(gc1);
    int nbf1 = int_shell1_->nfunction(gc1);

    int target_bf2_offset = 0;
    for (int gc2=0; gc2<int_shell2_->ncontraction(); gc2++) {
      int am2 = int_shell2_->am(gc2);
      int is_pure2 = int_shell2_->is_pure(gc2) ? 1 : 0;
      int ncart2 = int_shell2_->ncartesian(gc2);
      int nbf2 = int_shell2_->nfunction(gc2);
      
      int target_bf3_offset = 0;
      for (int gc3=0; gc3<int_shell3_->ncontraction(); gc3++) {
	int am3 = int_shell3_->am(gc3);
	int is_pure3 = int_shell3_->is_pure(gc3) ? 1 : 0;
	int ncart3 = int_shell3_->ncartesian(gc3);
	int nbf3 = int_shell3_->nfunction(gc3);
	
	int target_bf4_offset = 0;
	for (int gc4=0; gc4<int_shell4_->ncontraction(); gc4++) {
	  int am4 = int_shell4_->am(gc4);
	  int is_pure4 = int_shell4_->is_pure(gc4) ? 1 : 0;
	  int ncart4 = int_shell4_->ncartesian(gc4);
	  int nbf4 = int_shell4_->nfunction(gc4);

	  int transform_index = 8*is_pure1 + 4*is_pure2 + 2*is_pure3 + is_pure4;
	  switch (transform_index) {
	    case 0:
	      break;

	    case 1:
	      source4 = source;
	      target4 = target;
	      break;

	    case 2:
	      source3 = source;
	      target3 = target;
	      break;

	    case 3:
	      source4 = source;
	      target4 = tmpbuf;
	      source3 = tmpbuf;
	      target3 = target;
	      break;
	      
	    case 4:
	      source2 = source;
	      target2 = target;
	      break;

	    case 5:
	      source4 = source;
	      target4 = tmpbuf;
	      source2 = tmpbuf;
	      target2 = target;
	      break;

	    case 6:
	      source3 = source;
	      target3 = tmpbuf;
	      source2 = tmpbuf;
	      target2 = target;
	      break;
	    
	    case 7:
	      source4 = source;
	      target4 = tmpbuf;
	      source3 = tmpbuf;
	      target3 = source;
	      source2 = source;
	      target2 = target;
	      break;

	    case 8:
	      source1 = source;
	      target1 = target;
	      break;
	      
	    case 9:
	      source4 = source;
	      target4 = tmpbuf;
	      source1 = tmpbuf;
	      target1 = target;
	      break;

	    case 10:
	      source3 = source;
	      target3 = tmpbuf;
	      source1 = tmpbuf;
	      target1 = target;
	      break;

	    case 11:
	      source4 = source;
	      target4 = tmpbuf;
	      source3 = tmpbuf;
	      target3 = source;
	      source1 = source;
	      target1 = target;
	      break;
	      
	    case 12:
	      source2 = source;
	      target2 = tmpbuf;
	      source1 = tmpbuf;
	      target1 = target;
	      break;

	    case 13:
	      source4 = source;
	      target4 = tmpbuf;
	      source2 = tmpbuf;
	      target2 = source;
	      source1 = source;
	      target1 = target;
	      break;

	    case 14:
	      source3 = source;
	      target3 = tmpbuf;
	      source2 = tmpbuf;
	      target2 = source;
	      source1 = source;
	      target1 = target;
	      break;
	    
	    case 15:
	      source4 = source;
	      target4 = tmpbuf;
	      source3 = tmpbuf;
	      target3 = source;
	      source2 = source;
	      target2 = tmpbuf;
	      source1 = tmpbuf;
	      target1 = target;
	      break;
	  }

	  if (is_pure4) {
	    SphericalTransformIter stiter(integral_->spherical_transform(am4));
	    transform2e_4(stiter, source4, target4, ncart1*ncart2*ncart3,ncart4);
	  }
	  if (is_pure3) {
	    SphericalTransformIter stiter(integral_->spherical_transform(am3));
	    transform2e_3(stiter,source3, target3, ncart1*ncart2,ncart3,nbf4);
	  }
	  if (is_pure2) {
	    SphericalTransformIter stiter(integral_->spherical_transform(am2));
	    transform2e_2(stiter,source2, target2, ncart1,ncart2,nbf3*nbf4);
	  }
	  if (is_pure1) {
	    SphericalTransformIter stiter(integral_->spherical_transform(am1));
	    transform2e_1(stiter,source1, target1, nbf2*nbf3*nbf4);
	  }
	  
	  source += (ncart1*ncart2*ncart3*ncart4);
	  target += (nbf1*nbf2*nbf3*nbf4);
	}
      }
    }
  }
}


static void transform1e_1(SphericalTransformIter& sti, double *s, double *t, int nl)
{
  memset(t,0,sti.n()*nl*sizeof(double));

  for (sti.begin(); sti.ready(); sti.next()) {
    double *sptr = s + sti.cartindex()*nl;
    double *tptr = t + sti.pureindex()*nl;
    double coef = sti.coef();
    for(int l=0; l<nl; l++)
      *(tptr++) += coef * *(sptr++);
  }
}

static void transform1e_2(SphericalTransformIter& sti, double *s, double *t, int nk, int nl)
{
  const int sl = nl;
  const int tl = sti.n();

  memset(t,0,nk*tl*sizeof(double));

  for (sti.begin(); sti.ready(); sti.next()) {
    double *sptr = s + sti.cartindex();
    double *tptr = t + sti.pureindex();
    double coef = sti.coef();
    for(int k=0; k<nk; k++,sptr+=sl,tptr+=tl) {
	*(tptr) += coef * *(sptr);
    }
  }
}

static void transform1e_vec_2(const int ntypes, SphericalTransformIter& sti, double *s, double *t, int nk, int nl)
{
  const int sl = nl*ntypes;
  const int tl = sti.n()*ntypes;

  memset(t,0,nk*tl*sizeof(double));

  for (sti.begin(); sti.ready(); sti.next()) {
    double *sptr = s + sti.cartindex()*ntypes;
    double *tptr = t + sti.pureindex()*ntypes;
    double coef = sti.coef();
    for(int k=0; k<nk; k++,sptr+=sl,tptr+=tl) {
      for(int type=0; type<ntypes; type++)
	tptr[type] += coef * sptr[type];
    }
  }
}

static void transform2e_1(SphericalTransformIter& sti, double *s, double *t, int njkl)
{
  memset(t,0,sti.n()*njkl*sizeof(double));

  for (sti.begin(); sti.ready(); sti.next()) {
    double *sptr = s + sti.cartindex()*njkl;
    double *tptr = t + sti.pureindex()*njkl;
    double coef = sti.coef();
    for(int jkl=0; jkl<njkl; jkl++)
      *(tptr++) += coef * *(sptr++);
  }
}

static void transform2e_2(SphericalTransformIter& sti, double *s, double *t, int ni, int nj, int nkl)
{
  int sj = sti.n();
  const int sjkl = nj*nkl;
  const int tjkl = sj*nkl;

  memset(t,0,ni*tjkl*sizeof(double));

  for (sti.begin(); sti.ready(); sti.next()) {
    double *sptr = s + sti.cartindex()*nkl;
    double *tptr = t + sti.pureindex()*nkl;
    double coef = sti.coef();
    for(int i=0; i<ni; i++,sptr+=sjkl,tptr+=tjkl) {
      for(int kl=0; kl<nkl; kl++)
	tptr[kl] += coef * sptr[kl];
    }
  }
}

static void transform2e_3(SphericalTransformIter& sti, double *s, double *t, int nij, int nk, int nl)
{
  int sk = sti.n();
  const int skl = nk*nl;
  const int tkl = sk*nl;

  memset(t,0,nij*tkl*sizeof(double));

  for (sti.begin(); sti.ready(); sti.next()) {
    double *sptr = s + sti.cartindex()*nl;
    double *tptr = t + sti.pureindex()*nl;
    double coef = sti.coef();
    for(int ij=0; ij<nij; ij++,sptr+=skl,tptr+=tkl) {
      for(int l=0; l<nl; l++)
	tptr[l] += coef * sptr[l];
    }
  }
}

static void transform2e_4(SphericalTransformIter& sti, double *s, double *t, int nijk, int nl)
{
  const int sl = nl;
  const int tl = sti.n();

  memset(t,0,nijk*tl*sizeof(double));

  for (sti.begin(); sti.ready(); sti.next()) {
    double *sptr = s + sti.cartindex();
    double *tptr = t + sti.pureindex();
    double coef = sti.coef();
    for(int ijk=0; ijk<nijk; ijk++,sptr+=sl,tptr+=tl) {
	*(tptr) += coef * *(sptr);
    }
  }
}


/*!------------------------------------------------------------------------------------
    Normalizes cartesian components according to appropriate convention (see GTOInfo)
 ------------------------------------------------------------------------------------*/
void
Int1eLibint2::norm_contrcart1_(double* data) {
  this->norm_contrcart_<1u>(data);
}

void
Int2eLibint2::norm_contrcart1_(double *data)
{
  this->norm_contrcart_<1u>(data);
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:

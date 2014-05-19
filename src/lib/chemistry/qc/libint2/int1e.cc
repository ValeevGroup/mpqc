//
// int1e.cc
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

#include <chemistry/qc/libint2/int1e.h>
#include <chemistry/qc/libint2/macros.h>

using namespace std;
using namespace sc;

inline int max(int a,int b) { return (a > b) ? a : b;}
inline void fail()
{
  ExEnv::errn() << scprintf("failing module:\n%s",__FILE__) << endl;
  abort();
}


Int1eLibint2::Int1eLibint2(Integral *integral,
		       const Ref<GaussianBasisSet>&b1,
		       const Ref<GaussianBasisSet>&b2,
		       int order, bool need_overlap, bool need_coulomb,
		       int ntypes) :
  integral_(integral), bs1_(b1), bs2_(b2), operset_params_(),
  EdotV_origin_(0), Q_origin_(0), need_overlap_(need_overlap),
  need_coulomb_(need_coulomb), ntypes_(ntypes)
{
 if (order > 0) {
    // Complain here
  }

  max_doublet_size_ = bs1_->max_nfunction_in_shell() * bs2_->max_nfunction_in_shell();
  target_ints_buffer_ = new double[ntypes_*max_doublet_size_];

  max_cart_doublet_size_ = bs1_->max_ncartesian_in_shell() * bs2_->max_ncartesian_in_shell();
  // These are target integrals in Cartesian basis and in by-contraction-doublets order
  cart_ints_ = new double[ntypes_*max_cart_doublet_size_];
  if (bs1_->has_pure() || bs2_->has_pure() || bs1_->max_ncontraction() != 1 || bs2_->max_ncontraction() != 1) {
    // These are target integrals in spherical harmonics basis and in by-contraction-doublets order
    sphharm_ints_ = new double[ntypes_*max_doublet_size_];
    // compute how much space one contraction doublet may need
    int nshell1 = bs1_->nshell();
    int maxncart1 = 0;
    for(int sh1=0; sh1<nshell1;sh1++) {
      int maxncart = bs1_->shell(sh1).max_cartesian();
      if (maxncart > maxncart1) maxncart1 = maxncart;
    }
    int nshell2 = bs2_->nshell();
    int maxncart2 = 0;
    for(int sh2=0; sh2<nshell2;sh2++) {
      int maxncart = bs2_->shell(sh2).max_cartesian();
      if (maxncart > maxncart2) maxncart2 = maxncart;
    }
    tformbuf_ = new double[ntypes_*maxncart1*maxncart2];
  }
  else {
    sphharm_ints_ = 0;
    tformbuf_ = 0;
  }

  int max_am = max(bs1_->max_angular_momentum(),bs2_->max_angular_momentum());
  if (need_overlap_) {
    // Allocate OIXYZ
    // max_am+1 - the range of exponents of x, y, and z
    // 2 - to get kinetic energy or p4 integrals
    // order - to allow for derivatives
    OIX_ = init_block_(max_am+1+2+order,max_am+1+2+order);
    OIY_ = init_block_(max_am+1+2+order,max_am+1+2+order);
    OIZ_ = init_block_(max_am+1+2+order,max_am+1+2+order);
  }

  if (need_coulomb_) {
    int efield_order = 0;
    if (ntypes_ > 1) efield_order = 1;
    if (ntypes_ > 3) efield_order = 2;
    const int mmax = bs1_->max_angular_momentum() + bs2_->max_angular_momentum() + order + efield_order;
    Fm_Eval_ = CoreIntsEngine<_FmEvalType>::instance(mmax);
    Fm_table_ = new double[mmax+1];
    indmax_ = (max_am+order)*(max_am+1+order)*(max_am+1+order)+1;
    // Allocate AI0
    AI0_ = init_box_(indmax_,indmax_,2*(max_am+order)+3);
    if (efield_order) {
      AIX_ = init_box_(indmax_,indmax_,2*(max_am+order)+2);
      AIY_ = init_box_(indmax_,indmax_,2*(max_am+order)+2);
      AIZ_ = init_box_(indmax_,indmax_,2*(max_am+order)+2);
    }
    if (efield_order) {
      AIXX_ = init_box_(indmax_,indmax_,2*(max_am+order)+1);
      AIXY_ = init_box_(indmax_,indmax_,2*(max_am+order)+1);
      AIXZ_ = init_box_(indmax_,indmax_,2*(max_am+order)+1);
      AIYY_ = init_box_(indmax_,indmax_,2*(max_am+order)+1);
      AIYZ_ = init_box_(indmax_,indmax_,2*(max_am+order)+1);
      AIZZ_ = init_box_(indmax_,indmax_,2*(max_am+order)+1);
    }
  }
}

Int1eLibint2::~Int1eLibint2()
{
  delete[] cart_ints_;
  if (sphharm_ints_) {
    delete[] sphharm_ints_;
    sphharm_ints_ = 0;
  }
  if (tformbuf_) {
    delete[] tformbuf_;
    tformbuf_ = 0;
  }
  if (need_coulomb_) {
    free_box_(AI0_);
    if (ntypes_ > 1) {
      free_box_(AIX_);
      free_box_(AIY_);
      free_box_(AIZ_);
    }
    if (ntypes_ > 3) {
      free_box_(AIXX_);
      free_box_(AIXY_);
      free_box_(AIXZ_);
      free_box_(AIYY_);
      free_box_(AIYZ_);
      free_box_(AIZZ_);
    }
    delete[] Fm_table_;
  }
  if (need_overlap_) {
    free_block_(OIX_);
    free_block_(OIY_);
    free_block_(OIZ_);
  }
  delete[] target_ints_buffer_;
}

void Int1eLibint2::set_params(const Ref<IntParams>& params)
{
  operset_params_ = params;
}

void Int1eLibint2::set_EdotV_origin(const Ref<EfieldDotVectorData>& origin)
{
  EdotV_origin_ = origin;
}

void Int1eLibint2::set_Q_origin(const Ref<PointChargeData>& origin)
{
  Q_origin_ = origin;
}

Ref<IntParams>
Int1eLibint2::params()
{
  return operset_params_;
}

Ref<IntParamsOrigin>
Int1eLibint2::origin()
{
  return require_dynamic_cast<IntParamsOrigin*>(operset_params_.pointer(), "need multipole origin, but not in multipole evaluator");
}

Ref<EfieldDotVectorData>
Int1eLibint2::EdotV_origin()
{
  return EdotV_origin_;
}

Ref<PointChargeData>
Int1eLibint2::Q_origin()
{
  return Q_origin_;
}

void Int1eLibint2::zero_buffers_()
{
  double *buf1 = cart_ints_;
  for(int i=0; i<max_cart_doublet_size_; i++,buf1++) {
    *buf1 = 0.0;
  }
  buf1 = target_ints_buffer_;
  for(int i=0; i<max_doublet_size_; i++,buf1++) {
    *buf1 = 0.0;
  }
}

void Int1eLibint2::compute_doublet_info_(int sh1, int sh2)
{
  int_shell1_ = &bs1_->shell(sh1);
  int_shell2_ = &bs2_->shell(sh2);
  int ctr1 = bs1_->shell_to_center(sh1);
  int ctr2 = bs2_->shell_to_center(sh2);
  doublet_info_.AB2 = 0.0;
  for(int i=0; i<3; i++) {
    doublet_info_.A[i] = bs1_->r(ctr1,i);
    doublet_info_.B[i] = bs2_->r(ctr2,i);
    doublet_info_.AB2 += (doublet_info_.A[i] - doublet_info_.B[i])*
      (doublet_info_.A[i] - doublet_info_.B[i]);
  }
}

void Int1eLibint2::sort_contrdoublets_to_shelldoublet_(double *source, double *target)
{
  /*--- sort to the target ordering ---*/
  double *source_ints_buf = source;
  double *target_ints_buf = target;
  int target_bf1_offset = 0;
  int nbf2 = int_shell2_->nfunction();
  for (int gc1=0; gc1<int_shell1_->ncontraction(); gc1++) {
    int am1 = int_shell1_->am(gc1);
    int tsize1 = int_shell1_->nfunction(gc1);
    int target_bf2_offset = 0;
    for (int gc2=0; gc2<int_shell2_->ncontraction(); gc2++) {
      int am2 = int_shell2_->am(gc2);
      int tsize2 = int_shell2_->nfunction(gc2);
      
      for(int bf1=0;bf1<tsize1;bf1++) {
	double *target_ints_buf = target_ints_buffer_ + (target_bf1_offset+bf1)*nbf2 +
	                          target_bf2_offset;
	for(int bf2=0;bf2<tsize2;bf2++) {
  	  *(target_ints_buf++) = *(source_ints_buf++);
	}
      }
      target_bf2_offset += tsize2;
    }
    target_bf1_offset += tsize1;
  }
}

void Int1eLibint2::zero_buffers_vec_(const int ntypes)
{
  double *buf1 = cart_ints_;
  const int nlibint2 = ntypes * max_cart_doublet_size_;
  for(int i=0; i<nlibint2; i++,buf1++) {
    *buf1 = 0.0;
  }
  buf1 = target_ints_buffer_;
  const int nints = ntypes * max_doublet_size_;
  for(int i=0; i<nints; i++,buf1++) {
    *buf1 = 0.0;
  }
}

void Int1eLibint2::sort_contrdoublets_to_shelldoublet_vec_(const int ntypes, double *source, double *target)
{
  /*--- sort to the target ordering ---*/
  double *source_ints_buf = source;
  double *target_ints_buf = target;
  int target_bf1_offset = 0;
  int nbf2 = int_shell2_->nfunction();
  for (int gc1=0; gc1<int_shell1_->ncontraction(); gc1++) {
    int am1 = int_shell1_->am(gc1);
    int tsize1 = int_shell1_->nfunction(gc1);
    int target_bf2_offset = 0;
    for (int gc2=0; gc2<int_shell2_->ncontraction(); gc2++) {
      int am2 = int_shell2_->am(gc2);
      int tsize2 = int_shell2_->nfunction(gc2);
      
      for(int bf1=0;bf1<tsize1;bf1++) {
	double *target_ints_buf = target_ints_buffer_ + ((target_bf1_offset+bf1)*nbf2 +
	                          target_bf2_offset)*ntypes;
	for(int bf2=0;bf2<tsize2;bf2++) {
	  for(int type=0; type<ntypes; type++)
	    *(target_ints_buf++) = *(source_ints_buf++);
	}
      }
      target_bf2_offset += tsize2;
    }
    target_bf1_offset += tsize1;
  }
}

double **Int1eLibint2::init_block_(int a, int b)
{
  double **block = new double*[a];
  block[0] = new double[a*b];
  for(int i=1; i<a; i++)
    block[i] = block[i-1] + b;

  return block;
}

void Int1eLibint2::free_block_(double **block)
{
  delete[] block[0];
  delete[] block;
}

double ***Int1eLibint2::init_box_(int a, int b, int c)
{
  int i,j;

  double ***box = new double**[a];
  box[0] = new double*[a*b];
  for(i=1; i<a; i++)
    box[i] = box[i-1] + b;

  box[0][0] = new double[a*b*c];
  for(j=1; j<b; j++)
    box[0][j] = box[0][j-1] + c;
  for(i=1; i<a; i++) {
    box[i][0] = box[i-1][b-1] + c;
    for(j=1; j<b; j++)
      box[i][j] = box[i][j-1] + c;
  }

  return box;
}

void Int1eLibint2::free_box_(double ***box)
{
  delete[] box[0][0];
  delete[] box[0];
  delete[] box;
}


/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:

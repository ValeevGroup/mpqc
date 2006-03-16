//
// int1e.cc
//
// Copyright (C) 2004 Sandia National Laboratories.
//
// Author: Joseph Kenny <jpkenny@sandia.gov>
// Maintainer: Joseph Kenny
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

#ifdef __GNUG__
#pragma implementation
#endif

#include <chemistry/qc/intcca/int1e.h>
#include <util/class/scexception.h>
#include <Chemistry_Chemistry_QC_GaussianBasis_DerivCenters.hh>
#include <Chemistry_QC_GaussianBasis_ObIntEvalType.hh>

using namespace std;
using namespace sc;
using namespace Chemistry;
using namespace Chemistry::QC::GaussianBasis;

Int1eCCA::Int1eCCA(Integral *integral,
		   const Ref<GaussianBasisSet>&b1,
		   const Ref<GaussianBasisSet>&b2,
		   int order, IntegralEvaluatorFactory eval_factory, 
                   ObIntEvalType int_type, bool use_opaque,
                   Chemistry::QC::GaussianBasis::DerivCenters cca_dc ):
  bs1_(b1), bs2_(b2),
  overlap_ptr_(0), kinetic_ptr_(0),
  nuclear_ptr_(0), hcore_ptr_(0),
  integral_(integral), eval_factory_(eval_factory), use_opaque_(use_opaque),
  cca_dc_(cca_dc)
{

  int scratchsize=0,nshell2;
  
  /* The efield routines look like derivatives so bump up order if
   * it is zero to allow efield integrals to be computed.
   */
  deriv_lvl_ = order;
  if (order == 0) order = 1;

  nshell2 = bs1_->max_ncartesian_in_shell()*bs2_->max_ncartesian_in_shell();

  if (order == 0)
    scratchsize = nshell2;
  else if (order == 1) 
    scratchsize = nshell2*3;
  else
    throw InputError("invalid derivative level",
                     __FILE__,__LINE__);

  if( !use_opaque_ )
    buff_ = new double[scratchsize];
  if( deriv_lvl_ )
    temp_buff_ = new double[scratchsize];

  // create cca basis sets
  cca_bs1_ = GaussianBasis_Molecular::_create();
  cca_bs1_.initialize( bs1_.pointer(), bs1_->name() );
  if( bs1_.pointer() != bs2_.pointer() ) {
    cca_bs2_ = GaussianBasis_Molecular::_create();
    cca_bs2_.initialize( bs2_.pointer(), bs2_->name() );
  }
  else
    cca_bs2_ = cca_bs1_;

  if( int_type == ObIntEvalType_OVERLAP && deriv_lvl_ == 0 ) {
    overlap_ = eval_factory_.get_obint_evaluator2( int_type, 0, 
						   cca_bs1_, cca_bs2_,
						   cca_dc_ );
    overlap_ptr_ = &overlap_;
    if( use_opaque_ )
      buff_ = static_cast<double*>( overlap_ptr_->get_obint_buffer() );
  }

  else if( int_type == ObIntEvalType_OVERLAP && deriv_lvl_ == 1 ) {
    overlap_1der_ = eval_factory_.get_obint_evaluator2( int_type, 1,
							cca_bs1_, cca_bs2_,
							cca_dc_ );
    overlap_1der_ptr_ = &overlap_1der_;
    if( use_opaque_ )
      buff_ = static_cast<double*>( overlap_1der_ptr_->get_obint_buffer() );
  }

  else if( int_type == ObIntEvalType_KINETIC && deriv_lvl_ == 0 ) {
    kinetic_ = eval_factory_.get_obint_evaluator2( int_type, 0,
                                                      cca_bs1_, cca_bs2_,
                                                      cca_dc_ );
    kinetic_ptr_ = &kinetic_;
    if( use_opaque_ )
      buff_ = static_cast<double*>( kinetic_ptr_->get_obint_buffer() );
  }

  else if( int_type == ObIntEvalType_KINETIC && deriv_lvl_ == 1 ) {
    kinetic_1der_ = eval_factory_.get_obint_evaluator2( int_type, 1,
                                                           cca_bs1_, cca_bs2_,
                                                           cca_dc_ );
    kinetic_1der_ptr_ = &kinetic_1der_;
    if( use_opaque_ )
      buff_ = static_cast<double*>( kinetic_1der_ptr_->get_obint_buffer() );
  }

  else if( int_type == ObIntEvalType_NUCLEAR && deriv_lvl_ == 0 ) {
    nuclear_ = eval_factory_.get_obint_evaluator2( int_type, 0,
                                                      cca_bs1_, cca_bs2_,
                                                      cca_dc_ );
    nuclear_ptr_ = &nuclear_;
    if( use_opaque_ )
      buff_ = static_cast<double*>( nuclear_ptr_->get_obint_buffer() );
  }

  else if( int_type == ObIntEvalType_NUCLEAR && deriv_lvl_ == 1 ) {
    nuclear_1der_ = eval_factory_.get_obint_evaluator2( int_type, 1,
                                                           cca_bs1_, cca_bs2_,
                                                           cca_dc_ );
    nuclear_1der_ptr_ = &nuclear_1der_;
    if( use_opaque_ )
      buff_ = static_cast<double*>( nuclear_1der_ptr_->get_obint_buffer() );
  }

  else if( int_type == ObIntEvalType_OEHAM && deriv_lvl_ == 0 ) {
    hcore_ = eval_factory_.get_obint_evaluator2( int_type, 0,
                                                    cca_bs1_, cca_bs2_,
                                                    cca_dc_ );
    hcore_ptr_ = &hcore_;
    if( use_opaque_ )
      buff_ = static_cast<double*>( hcore_ptr_->get_obint_buffer() );
  }

  else if( int_type == ObIntEvalType_OEHAM && deriv_lvl_ == 1 ) {
    hcore_1der_ = eval_factory_.get_obint_evaluator2( int_type, 1,
                                                         cca_bs1_, cca_bs2_,
                                                         cca_dc_ );
    hcore_1der_ptr_ = &hcore_1der_;
    if( use_opaque_ )
      buff_ = static_cast<double*>( hcore_1der_ptr_->get_obint_buffer() );
  }

}

Int1eCCA::~Int1eCCA()
{
}

void
Int1eCCA::overlap( int ish, int jsh )
{
  cca_dc_.clear();
  if( use_opaque_ )
    overlap_ptr_->compute( ish, jsh, 0, -1 );

  else {
    sidl_buffer_ = overlap_ptr_->compute_array( ish, jsh, 0, -1 ); 
    copy_buffer();
  }
}  

void
Int1eCCA::overlap_1der(int ish, int jsh)
{
  if( use_opaque_ ) 
    overlap_1der_ptr_->compute( ish, jsh, 1, -1 );
  else {
    sidl_buffer_ = overlap_1der_ptr_->compute_array( ish, jsh, 1, -1 );
    copy_buffer();
  }
#ifndef INTV3_ORDER
  reorder_deriv( &(bs1_->shell(ish)), &(bs2_->shell(jsh)) );
#endif
}

void
Int1eCCA::overlap_1der(int ish, int jsh, int c)
{
  cca_dc_.clear();
  if( use_opaque_ )
    overlap_1der_ptr_->compute( ish, jsh, 1, c );
  else {
    sidl_buffer_ = overlap_1der_ptr_->compute_array( ish, jsh, 1, c );
    copy_buffer();
  }
#ifndef INTV3_ORDER
  reorder_deriv( &(bs1_->shell(ish)), &(bs2_->shell(jsh)) );
#endif
}

void
Int1eCCA::kinetic( int ish, int jsh )
{
  cca_dc_.clear();
  if( use_opaque_ )
    kinetic_ptr_->compute( ish, jsh, 0, -1 );
  else {
    sidl_buffer_ = kinetic_ptr_->compute_array( ish, jsh, 0, -1 ); 
    copy_buffer();
  }
}

void
Int1eCCA::kinetic_1der(int ish, int jsh )
{
  if( use_opaque_ )
    kinetic_1der_ptr_->compute( ish, jsh, 1, -1 );
  else {
    sidl_buffer_ = kinetic_1der_ptr_->compute_array( ish, jsh, 1, -1 );
    copy_buffer();
  }
#ifndef INTV3_ORDER
  reorder_deriv( &(bs1_->shell(ish)), &(bs2_->shell(jsh)) );
#endif
}

void
Int1eCCA::kinetic_1der(int ish, int jsh, int c)
{
  cca_dc_.clear();
  if( use_opaque_ )
    kinetic_1der_ptr_->compute( ish, jsh, 1, c );
  else {
    sidl_buffer_ = kinetic_1der_ptr_->compute_array( ish, jsh, 1, c );
    copy_buffer();
  }
#ifndef INTV3_ORDER
  reorder_deriv( &(bs1_->shell(ish)), &(bs2_->shell(jsh)) );
#endif
}


void
Int1eCCA::nuclear( int ish, int jsh )
{
  cca_dc_.clear();
  if( use_opaque_ )
    nuclear_ptr_->compute( ish, jsh, 0, -1 );
  else {
    sidl_buffer_ = nuclear_ptr_->compute_array( ish, jsh, 0, -1 ); 
    copy_buffer();
  }
}

void
Int1eCCA::nuclear_1der(int ish, int jsh, int c)
{
  cca_dc_.clear();
  if( use_opaque_ )
    nuclear_1der_ptr_->compute( ish, jsh, 1, c );
  else {
    sidl_buffer_ = nuclear_1der_ptr_->compute_array( ish, jsh, 1, c );
    copy_buffer();
  }
#ifndef INTV3_ORDER
  reorder_deriv( &(bs1_->shell(ish)), &(bs2_->shell(jsh)) );
#endif
}

void
Int1eCCA::nuclear_1der(int ish, int jsh)
{
  if( use_opaque_ )
    nuclear_1der_ptr_->compute( ish, jsh, 1, -1 );
  else {
    sidl_buffer_ = nuclear_1der_ptr_->compute_array( ish, jsh, 1, -1 );
    copy_buffer();
  }
#ifndef INTV3_ORDER
  reorder_deriv( &(bs1_->shell(ish)), &(bs2_->shell(jsh)) );
#endif
}

void
Int1eCCA::hcore( int ish, int jsh )
{
  cca_dc_.clear();
  if( use_opaque_ )
    hcore_ptr_->compute( ish, jsh, 0, -1 );
  else {
    sidl_buffer_ = hcore_ptr_->compute_array( ish, jsh, 0, -1 ); 
    copy_buffer();
  }
}

void
Int1eCCA::hcore_1der(int ish, int jsh, int c)
{
  cca_dc_.clear();
  if( use_opaque_ )
    hcore_1der_ptr_->compute( ish, jsh, 1, c );
  else {
    sidl_buffer_ = hcore_1der_ptr_->compute_array( ish, jsh, 1, c );
    copy_buffer();
  }
#ifndef INTV3_ORDER
  reorder_deriv( &(bs1_->shell(ish)), &(bs2_->shell(jsh)) );
#endif
}

void
Int1eCCA::hcore_1der(int ish, int jsh)
{
  if( use_opaque_ )
    hcore_1der_ptr_->compute( ish, jsh, 1, -1 );
  else {
    sidl_buffer_ = hcore_1der_ptr_->compute_array( ish, jsh, 1, -1 );
    copy_buffer();
  }
#ifndef INTV3_ORDER
  reorder_deriv( &(bs1_->shell(ish)), &(bs2_->shell(jsh)) );
#endif
}

void 
Int1eCCA::copy_buffer() 
{
  int sidl_size = 1 + sidl_buffer_.upper(0) - sidl_buffer_.lower(0);
  for(int i=0; i<sidl_size; ++i) 
    buff_[i] = sidl_buffer_.get(i);
}

void
Int1eCCA::reorder_deriv(sc::GaussianShell* s1, sc::GaussianShell* s2)
{
  int nfunc = s1->nfunction() * s2->nfunction();

  for(int i=0; i<nfunc*3; ++i)
    temp_buff_[i] = buff_[i];

  for(int i=0; i<nfunc; ++i)
    for( int di=0; di<3; ++di)
      buff_[i*3+di] = temp_buff_[di*nfunc+i];

/*
  std::cerr << "Int1eCCA reordering" << std::endl;
  for( int i=0; i<nfunc; ++i) {
    std::cerr << "integral " << i << std::endl;
    std::cerr << " dx: " << buff_[i*3] << std::endl;
    std::cerr << " dy: " << buff_[i*3+1] << std::endl;
    std::cerr << " dz: " << buff_[i*3+2] << std::endl;
  }
*/

}
    
  
  
  


/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// End:

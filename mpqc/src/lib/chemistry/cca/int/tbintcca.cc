//
// tbintcca.cc
//
// Copyright (C) 2004 Sandia National Laboratories
//
// Author: Joe Kenny <jpkenny@sandia.gov>
// Maintainer: Joe Kenny
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

#include "tbintcca.h"
#include <util/class/scexception.h>
#include <limits.h>

#include <Chemistry_Eri4IntegralDescr.hh>
#include <Chemistry_R12IntegralDescr.hh>
#include <Chemistry_R12T1IntegralDescr.hh>
#include <Chemistry_R12T2IntegralDescr.hh>

using namespace std;
using namespace sc;
using namespace Chemistry;
using namespace Chemistry::QC::GaussianBasis;

////////////////////////////////////////////////////////////////////////////
// TwoBodyIntCCA

TwoBodyIntCCA::TwoBodyIntCCA(Integral* integral,
                             const Ref<GaussianBasisSet> &bs1,
                             const Ref<GaussianBasisSet> &bs2,
			     const Ref<GaussianBasisSet> &bs3,
			     const Ref<GaussianBasisSet> &bs4,
			     IntegralSuperFactory fac,
			     CompositeIntegralDescr cdesc,
                             bool use_opaque) :
  TwoBodyInt(integral,bs1,bs2,bs3,bs4), 
  bs1_(bs1), bs2_(bs2), bs3_(bs3), bs4_(bs4),
  eval_factory_(fac), cdesc_(cdesc),
  use_opaque_(use_opaque)
{
  ndesc_ = cdesc_.get_n_descr();
  for( int i=0; i<ndesc_; ++i ) {
    descriptors_.push_back( cdesc_.get_descr(i) );
    types_.push_back( descriptors_.back().get_type() );
  }

  tbtype_to_buf_ = new double*[ndesc_];

  IntegralDescr desc = Chemistry::Eri4IntegralDescr::_create();
  dtype_to_tbtype_[desc.get_type()] = sc::TwoBodyInt::eri;
  desc = Chemistry::R12IntegralDescr::_create();
  dtype_to_tbtype_[desc.get_type()] = sc::TwoBodyInt::r12;
  desc = Chemistry::R12T1IntegralDescr::_create();
  dtype_to_tbtype_[desc.get_type()] = sc::TwoBodyInt::r12t1;
  desc = Chemistry::R12T2IntegralDescr::_create();
  dtype_to_tbtype_[desc.get_type()] = sc::TwoBodyInt::r12t2;

  int_bound_min_ = SCHAR_MIN;
  tol_ = pow(2.0,double(int_bound_min_));
  loginv_ = 1.0/log(2.0);

  int scratchsize = bs1_->max_ncartesian_in_shell()
    * bs2_->max_ncartesian_in_shell()
    * bs3_->max_ncartesian_in_shell()
    * bs4_->max_ncartesian_in_shell();
    
  if( !use_opaque_ ) buff_ = new double[scratchsize];

  // create cca basis sets
  cca_bs1_ = MPQC::GaussianBasis_Molecular::_create();
  cca_bs1_.initialize( bs1_.pointer(), bs1_->label() );
  if( bs1_.pointer() != bs2_.pointer() ) {
    cca_bs2_ = MPQC::GaussianBasis_Molecular::_create();
    cca_bs2_.initialize( bs2_.pointer(), bs2_->label() );
  }
  else
    cca_bs2_ = cca_bs1_;
  if( bs2_.pointer() != bs3_.pointer() ) {
    cca_bs3_ = MPQC::GaussianBasis_Molecular::_create();
    cca_bs3_.initialize( bs3_.pointer(), bs3_->label() );
  }
  else
    cca_bs3_ = cca_bs2_;
  if( bs3_.pointer() != bs4_.pointer() ) {
    cca_bs4_ = MPQC::GaussianBasis_Molecular::_create();
    cca_bs4_.initialize( bs4_.pointer(), bs4_->label() );
  }
  else
    cca_bs4_ = cca_bs3_;

  eval_ = eval_factory_.get_evaluator4( cdesc_, cca_bs1_, cca_bs2_, 
					cca_bs3_, cca_bs4_ );
  for( int i=0; i<ndesc_; ++i ) {
    IntegralDescr desc = cdesc_.get_descr(i);
    tbtype_to_buf_[ dtype_to_tbtype_[desc.get_type()] ]
      = static_cast<double*>( eval_.get_buffer(desc) );
  }

}

const double*
TwoBodyIntCCA::buffer(tbint_type te_type) const
{
  return tbtype_to_buf_[ te_type ];
}

TwoBodyIntCCA::~TwoBodyIntCCA()
{
  if( !use_opaque_ )
    delete buff_;
}

void
TwoBodyIntCCA::compute_shell( int i, int j, int k, int l )
{
  for( int ii=0; ii<ndesc_; ++ii ) {
    buffer_ = tbtype_to_buf_[ dtype_to_tbtype_[ types_[ii] ] ];
    eval_.compute( i, j, k, l );
    if( !redundant_ )
      remove_redundant( i, j, k, l );
  }
}

int
TwoBodyIntCCA::log2_shell_bound( int i, int j, int k, int l )
{

  double dbound = eval_.compute_bounds(i,j,k,l);

  // mpqc shouldn't return close to 0, but other codes might
  if( dbound < tol_ )
    return int_bound_min_;

  return static_cast<int>( logb(dbound) );


// something like this could be required for benchmarking if rounding
// effects become problematic, seems ok so far, so I'll leave it out
// for performance
/*
  double value = eval_.compute_bounds(i,j,k,l);
  int upper;
  if (value > tol_) {
    double log_dbl = log(value)*loginv_;
    upper = static_cast<int>( ceil(log_dbl) );
    // nasty check for rounding effects
    if( upper > int_bound_min_ ) {
      int lower = upper - 1;
      double lower_dbl = static_cast<double>( lower ); 
      if( abs(log_dbl - lower_dbl) < 0.01 )
        upper = lower;
    }
  }
  else upper = int_bound_min_;
  return upper;
*/

}

unsigned int
TwoBodyIntCCA::num_tbint_types() const
{
  return ndesc_;
}

////////////////////////////////////////////////////////////////////////////
// TwoBodyDerivIntCCA

TwoBodyDerivIntCCA::TwoBodyDerivIntCCA( Integral* integral,
					const Ref<GaussianBasisSet> &bs1,
					const Ref<GaussianBasisSet> &bs2,
					const Ref<GaussianBasisSet> &bs3,
					const Ref<GaussianBasisSet> &bs4,
					IntegralSuperFactory fac,
					CompositeIntegralDescr cdesc,
					bool use_opaque ) :
  TwoBodyDerivInt(integral,bs1,bs2,bs3,bs4),
  bs1_(bs1), bs2_(bs2), bs3_(bs3), bs4_(bs4),
  eval_factory_(fac), cdesc_(cdesc),
  use_opaque_(use_opaque)
{
  ndesc_ = cdesc_.get_n_descr();
  for( int i=0; i<ndesc_; ++i )
    descriptors_.push_back( cdesc_.get_descr(i) );

  int_bound_min_ = SCHAR_MIN;
  tol_ = pow(2.0,double(int_bound_min_));
  loginv_ = 1.0/log(2.0);

  // find max deriv level
  max_deriv_lvl_=0;
  int n_descr = cdesc_.get_n_descr();
  for( int i=0; i<n_descr; ++i ) {
    int temp = cdesc_.get_descr(i).get_deriv_lvl();
    if( temp > max_deriv_lvl_ )
      max_deriv_lvl_ = temp; 
  }

  int scratchsize = bs1_->max_ncartesian_in_shell()
    * bs2_->max_ncartesian_in_shell()
    * bs3_->max_ncartesian_in_shell()
    * bs4_->max_ncartesian_in_shell();
  if( max_deriv_lvl_ == 1 )
    scratchsize *= 9;
  else if( max_deriv_lvl_ != 0 )
    throw FeatureNotImplemented("only first order derivatives are available",
                                __FILE__,__LINE__);
    
  if( !use_opaque_ ) buff_ = new double[scratchsize];


  // create cca basis sets
  cca_bs1_ = MPQC::GaussianBasis_Molecular::_create();
  cca_bs1_.initialize( bs1_.pointer(), bs1_->label() );
  if( bs1_.pointer() != bs2_.pointer() ) {
    cca_bs2_ = MPQC::GaussianBasis_Molecular::_create();
    cca_bs2_.initialize( bs2_.pointer(), bs2_->label() );
  }
  else
    cca_bs2_ = cca_bs1_;
  if( bs2_.pointer() != bs3_.pointer() ) {
    cca_bs3_ = MPQC::GaussianBasis_Molecular::_create();
    cca_bs3_.initialize( bs3_.pointer(), bs3_->label() );
  }
  else
    cca_bs3_ = cca_bs2_;
  if( bs3_.pointer() != bs4_.pointer() ) {
    cca_bs4_ = MPQC::GaussianBasis_Molecular::_create();
    cca_bs4_.initialize( bs4_.pointer(), bs4_->label() );
  }
  else
    cca_bs4_ = cca_bs3_;

  // set factory config
/*
  sidl::array<string> sidl_factories = sidl::array<string>::create1d(n_descr);
  for( int i=0; i<n_descr; ++i ) 
    sidl_factories.set( i, factories_[i] );
  eval_factory_.set_source_factories( sidl_factories );
*/

  IntegralDescr idesc = cdesc_.get_descr(0);
  cca_dc_ = idesc.get_deriv_centers();
  
  eval_ = eval_factory_.get_evaluator4( cdesc_, cca_bs1_, cca_bs2_, 
					cca_bs3_, cca_bs4_ );
  buffer_ = static_cast<double*>( eval_.get_buffer( cdesc_.get_descr(0) ) );
  // and what happens for multiple buffers???
}

TwoBodyDerivIntCCA::~TwoBodyDerivIntCCA() 
{
  //delete buff_;
}

void
TwoBodyDerivIntCCA::compute_shell( int i, int j, int k, int l,
                                  sc::DerivCenters &dc )
{
  cca_dc_.clear();
  if( use_opaque_ )
    eval_.compute( i, j, k, l );
  else {   
    sidl_buffer_ = eval_.compute_array( i, j, k, l );
    int nelem = bs1_->shell(i).nfunction() * bs2_->shell(j).nfunction() *
      bs3_->shell(k).nfunction() * bs4_->shell(l).nfunction() * 3;
    copy_buffer(nelem);
  }

  dc.clear();
  if( cca_dc_.has_omitted_center() )
    dc.add_omitted(cca_dc_.omitted_center(),cca_dc_.omitted_atom());
  for( int i=0; i<cca_dc_.n(); ++i)
    dc.add_center(cca_dc_.center(i),cca_dc_.atom(i));
}

int
TwoBodyDerivIntCCA::log2_shell_bound(int i, int j, int k, int l)
{
  double dbound = eval_.compute_bounds(i,j,k,l);

  // mpqc shouldn't return close to 0, but other codes might
  if( dbound < tol_ )
    return int_bound_min_;

  return static_cast<int>( logb(dbound) );


// something like this could be required for benchmarking if rounding
// effects become problematic, seems ok so far, so I'll leave it out
// for performance
/*
  double value = eval_.compute_bounds(i,j,k,l);
  int upper;
  if (value > tol_) {
    double log_dbl = log(value)*loginv_;
    upper = static_cast<int>( ceil(log_dbl) );
    // nasty check for rounding effects
    if( upper > int_bound_min_ ) {
      int lower = upper - 1;
      double lower_dbl = static_cast<double>( lower );
      if( abs(log_dbl - lower_dbl) < 0.01 )
        upper = lower;
    }
  }
  else upper = int_bound_min_;
  return upper;
*/

}

unsigned int
TwoBodyDerivIntCCA::num_tbint_types() const
{
  return 1;
}

void
TwoBodyDerivIntCCA::copy_buffer(int n) 
{
  for( int i=0; i<n; ++i)
    buff_[i] = sidl_buffer_.get(i);
}

#ifndef INTV3_ORDER

/////////////////////////////////////////////////////////////////////////////
// Code for removing redundant integrals
// copied liberally from cints

void get_nonredundant_ints_( double *source, double *target, 
			     int e13e24, int e12, int e34,
			     GaussianShell* int_shell1,
			     GaussianShell* int_shell2,
			     GaussianShell* int_shell3,
			     GaussianShell* int_shell4 )
{

  //std::cout << "\nremoving redundant integrals";

  int i,j,k,l;

  int nbf1 = int_shell1->nfunction();
  int nbf2 = int_shell2->nfunction();
  int nbf3 = int_shell3->nfunction();
  int nbf4 = int_shell4->nfunction();

  double *redundant_ptr = source;
  double *nonredundant_ptr = target;

  int nbf34 = nbf3*nbf4;
  for (i=0; i<nbf1; i++) {
    int jmax = e12 ? i : nbf2-1;
    for (j=0; j<=jmax; j++) {
      int kmax = e13e24 ? i : nbf3-1;
      for (k=0; k<=kmax; k++) {
        int lmax = e34 ? ( (e13e24&&(i==k)) ? j : k) : 
                         ( (e13e24&&(i==k)) ? j : nbf4-1);
        for (l=0; l<=lmax; l++) {
          *(nonredundant_ptr++) = redundant_ptr[l];
        }
        redundant_ptr += nbf4;
      }
      redundant_ptr += (nbf3-(kmax+1))*nbf4;
    }
    redundant_ptr += (nbf2-(jmax+1))*nbf34;
  }
}

void
TwoBodyIntCCA::remove_redundant(int sh1, int sh2, int sh3, int sh4) {

  GaussianShell* int_shell1(&bs1_->shell(sh1));
  GaussianShell* int_shell2(&bs2_->shell(sh2));
  GaussianShell* int_shell3(&bs3_->shell(sh3));
  GaussianShell* int_shell4(&bs4_->shell(sh4));

  bool need_unique_ints_only = false;
  int e12,e34,e13e24;
  e12 = 0;
  if (int_shell1 == int_shell2 && int_shell1->nfunction()>1)
    e12 = 1;
  e34 = 0;
  if (int_shell3 == int_shell4 && int_shell3->nfunction()>1)
    e34 = 1;
  e13e24 = 0;
  if (int_shell1 == int_shell3 && int_shell2 ==
      int_shell4 && int_shell1->nfunction()*int_shell2->nfunction()>1)
    e13e24 = 1;

  if ( e12 || e34 || e13e24 )
    need_unique_ints_only = true;

  if (need_unique_ints_only) {
    std::cout.flush();
    get_nonredundant_ints_( buffer_, buffer_, e13e24, e12, e34,
                            int_shell1, int_shell2, int_shell3, int_shell4 );
  }
}

#else

/////////////////////////////////////////////////////////////////////////////
// Code for removing redundant integrals
// copied liberally from intV3 

static int
shell_offset(Ref<GaussianBasisSet> cs, int off)
{
  return off + cs->nshell();
}

#define INT_MAX1(n1) ((n1)-1)
#define INT_MAX2(e12,i,n2) ((e12)?(i):((n2)-1))
#define INT_MAX3(e13e24,i,n3) ((e13e24)?(i):((n3)-1))
#define INT_MAX4(e13e24,e34,i,j,k,n4) \
  ((e34)?(((e13e24)&&((k)==(i)))?(j):(k)) \
        :((e13e24)&&((k)==(i)))?(j):(n4)-1)

void
nonredundant_erep(double *buffer, int e12, int e34, int e13e24,
		  int n1, int n2, int n3, int n4,
		  int *red_off, int *nonred_off)
{
  int nonredundant_index;
  int i,j,k,l;
  
  double *redundant_ptr = &buffer[*red_off];
  double *nonredundant_ptr = &buffer[*nonred_off];
  
  nonredundant_index = 0;
  int n34 = n3*n4;
  for (i=0; i<n1; i++) {
    int jmax = INT_MAX2(e12,i,n2);
    for (j=0; j<=jmax; j++) {
      int kmax = INT_MAX3(e13e24,i,n3);
      for (k=0; k<=kmax; k++) {
        int lmax = INT_MAX4(e13e24,e34,i,j,k,n4);
        for (l=0; l<=lmax; l++) {
          nonredundant_ptr[l] = redundant_ptr[l];
	}
        redundant_ptr += n4;
        nonredundant_index += lmax+1;
        nonredundant_ptr += lmax+1;
      }
      redundant_ptr += (n3-(kmax+1))*n4;
    }
    redundant_ptr += (n2-(jmax+1))*n34;
  }
  *red_off += n1*n2*n34;
  *nonred_off += nonredundant_index;
}

void
TwoBodyIntCCA::remove_redundant(int is, int js, int ks, int ls) {

  int bs1_shell_offset = 0;
  int bs2_shell_offset, bs3_shell_offset, bs4_shell_offset;
  
  int shell_offset1 = shell_offset(bs1_,0);
  int shell_offset2, shell_offset3;
  if (bs2_ == bs1_) {
    shell_offset2 = shell_offset1;
    bs2_shell_offset = bs1_shell_offset;
  }
  else {
    shell_offset2 = shell_offset(bs2_,shell_offset1);
    bs2_shell_offset = shell_offset1;
  }
  
  if (bs3_ == bs1_) {
    shell_offset3 = shell_offset2;
    bs3_shell_offset = bs1_shell_offset;
  }
  else if (bs3_ == bs2_) {
    shell_offset3 = shell_offset2;
    bs3_shell_offset = bs2_shell_offset;
  }
  else {
    shell_offset3 = shell_offset(bs3_,shell_offset2);
    bs3_shell_offset = shell_offset2;
  }
  
  if (bs4_ == bs1_) {
    bs4_shell_offset = bs1_shell_offset;
  }
  else if (bs4_ == bs2_) {
    bs4_shell_offset = bs2_shell_offset;
  }
  else if (bs4_ == bs3_) {
    bs4_shell_offset = bs3_shell_offset;
  }
  else {
    bs4_shell_offset = shell_offset3;
  }
  
  int int_unit2 = 0;
  int int_unit4 = 0;
  int osh1 = is + bs1_shell_offset;
  int osh2;
  if (!int_unit2) osh2 = js + bs2_shell_offset;
  int osh3 = ks + bs3_shell_offset;
  int osh4;
  if (!int_unit4) osh4 = ls + bs4_shell_offset;
  
  int sh1 = is;
  int sh2;
  if (!int_unit2) sh2 = js;
  int sh3 = ks;
  int sh4;
  if (!int_unit4) sh4 = ls;
  
  GaussianShell* int_unit_shell = 0;
  GaussianShell* int_shell1 = &bs1_->shell(sh1);
  GaussianShell* int_shell2;
  if (!int_unit2) int_shell2 = &bs2_->shell(sh2);
  else int_shell2 = int_unit_shell;
  GaussianShell* int_shell3 = &bs3_->shell(sh3);
  GaussianShell* int_shell4;
  if (!int_unit4) int_shell4 = &bs4_->shell(sh4);
  else int_shell4 = int_unit_shell;
  
  int redundant_offset = 0;
  int nonredundant_offset = 0;

  if ((osh1 == osh4)&&(osh2 == osh3)&&(osh1 != osh2)) {
    throw ProgrammingError("nonredundant integrals cannot be generated",
                           __FILE__,__LINE__);
  }

  int e12 = (int_unit2?0:(osh1 == osh2));
  int e13e24 = ((osh1 == osh3)
		&& ((int_unit2 && int_unit4)
		    || ((int_unit2||int_unit4)?0:(osh2 == osh4))));
  int e34 = (int_unit4?0:(osh3 == osh4));
  if (e12||e34||e13e24) {
    nonredundant_erep(buffer_,e12,e34,e13e24,
		      int_shell1->nfunction(),
		      int_shell2->nfunction(),
		      int_shell3->nfunction(),
		      int_shell4->nfunction(),
		      &redundant_offset,
		      &nonredundant_offset);
  }
}

#endif

/////////////////////////////////////////////////////////////////////////////
// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:

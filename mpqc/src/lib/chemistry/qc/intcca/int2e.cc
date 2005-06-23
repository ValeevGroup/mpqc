//
// int2e.cc
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

#include <chemistry/qc/intcca/int2e.h>
#include <Chemistry_Chemistry_QC_GaussianBasis_DerivCenters.hh>
#include <util/class/scexception.h>

using namespace std;
using namespace sc;
using namespace MPQC;
using namespace Chemistry;
using namespace Chemistry::QC::GaussianBasis;

Int2eCCA::Int2eCCA(Integral *integral,
		   const Ref<GaussianBasisSet>&b1,
		   const Ref<GaussianBasisSet>&b2,
		   const Ref<GaussianBasisSet>&b3,
		   const Ref<GaussianBasisSet>&b4,
		   int order, size_t storage,
		   IntegralEvaluatorFactory eval_factory, 
                   bool use_opaque, string eval_type ):
  bs1_(b1), bs2_(b2), bs3_(b3), bs4_(b4),
  erep_ptr_(0), integral_(integral), eval_factory_(eval_factory),
  use_opaque_(use_opaque), buffer_(0)
{

  /* Allocate storage for the integral buffer. */
  int maxsize = bs1_->max_ncartesian_in_shell()
                *bs2_->max_ncartesian_in_shell()
                *bs3_->max_ncartesian_in_shell()
                *bs4_->max_ncartesian_in_shell();
  if( order == 1 )
    maxsize *= 9;
  else if( order != 0 )
    throw FeatureNotImplemented("only first order derivatives are available",
                                __FILE__,__LINE__);
  if( !use_opaque ) buffer_ = new double[maxsize];

  cca_bs1_ = GaussianBasis_Molecular::_create();
  cca_bs2_ = GaussianBasis_Molecular::_create();
  cca_bs3_ = GaussianBasis_Molecular::_create();
  cca_bs4_ = GaussianBasis_Molecular::_create();
  cca_bs1_.initialize( bs1_.pointer(), bs1_->name() );
  cca_bs2_.initialize( bs2_.pointer(), bs2_->name() );
  cca_bs3_.initialize( bs3_.pointer(), bs3_->name() );
  cca_bs4_.initialize( bs4_.pointer(), bs4_->name() );
 
  if( eval_type == "eri" ) {
    erep_ = eval_factory_.get_integral_evaluator4( "eri2", 0,
                                                   cca_bs1_, cca_bs2_, 
                                                   cca_bs3_, cca_bs4_ );
    erep_ptr_ = &erep_;
    if( use_opaque_ )
      buffer_ = static_cast<double*>( erep_ptr_->get_buffer() );
  }
  else if( eval_type == "eri_deriv1") {
    erep_1der_ = eval_factory_.get_integral_evaluator4( "eri2", 1,
                                                        cca_bs1_, cca_bs2_,
                                                        cca_bs3_, cca_bs4_ );
    erep_1der_ptr_ = &erep_1der_;
    if( use_opaque_ )
      buffer_ = static_cast<double*>( erep_1der_ptr_->get_buffer() );
  }
  else
    throw InputError("unrecognized integral type",
                     __FILE__,__LINE__);
  if (!buffer_)
    throw ProgrammingError("buffer not assigned",
                           __FILE__,__LINE__);
}

void
Int2eCCA::compute_erep( int is, int js, int ks, int ls )
{
  Chemistry_QC_GaussianBasis_DerivCenters dc;
  dc = Chemistry_QC_GaussianBasis_DerivCenters::_create();
  if( use_opaque_ )
    erep_ptr_->compute( is, js, ks, ls, 0, dc );
  else {   
    sidl_buffer_ = erep_ptr_->compute_array( is, js, ks, ls, 0, dc );
    int nelem = bs1_->shell(is).nfunction() * bs2_->shell(js).nfunction() *
                bs3_->shell(ks).nfunction() * bs4_->shell(ls).nfunction();
    copy_buffer(nelem);
  }

  if(!redundant_) {
    remove_redundant(is,js,ks,ls);
  }
}  

void
Int2eCCA::compute_erep_1der( int is, int js, int ks, int ls, 
                             Chemistry::QC::GaussianBasis::DerivCenters &dc )
{

  if( use_opaque_ )
    erep_ptr_->compute( is, js, ks, ls, 1, dc );
  else {
    sidl_buffer_ = erep_ptr_->compute_array( is, js, ks, ls, 1, dc );
    int nelem = bs1_->shell(is).nfunction() * bs2_->shell(js).nfunction() *
                bs3_->shell(ks).nfunction() * bs4_->shell(ls).nfunction();
    copy_buffer(nelem);
  }

  if(!redundant_) {
    remove_redundant(is,js,ks,ls);
  }
}

void 
Int2eCCA::copy_buffer( int n ) 
{
  for( int i=0; i<n; ++i)
     buffer_[i] = sidl_buffer_.get(i);
}

#ifndef INTV3_ORDER

/////////////////////////////////////////////////////////////////////////////
// Code for removing redundant integrals
// copied liberally from cints

void get_nonredundant_ints_(double *source, double *target, 
                            int e13e24, int e12, int e34,
                            GaussianShell* int_shell1,GaussianShell* int_shell2,
                            GaussianShell* int_shell3,GaussianShell* int_shell4)
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
Int2eCCA::remove_redundant(int sh1, int sh2, int sh3, int sh4) {

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
Int2eCCA::remove_redundant(int is, int js, int ks, int ls) {

  std::cout << "\nREMOVING REDUNDANT INTEGRALS";

  int bs1_shell_offset = 0;
  int bs2_shell_offset, bs3_shell_offset, bs4_shell_offset;
  
  int shell_offset1 = shell_offset(bs1_,0);
  int shell_offset2, shell_offset3, shell_offset4;
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
// End:

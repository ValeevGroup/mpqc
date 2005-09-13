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

#include <chemistry/qc/intcca/tbintcca.h>
#include <Chemistry_Chemistry_QC_GaussianBasis_DerivCenters.hh>
#include <util/class/scexception.h>
#include <limits.h>

using namespace Chemistry::QC::GaussianBasis;
using namespace sc;

////////////////////////////////////////////////////////////////////////////
// TwoBodyIntCCA

TwoBodyIntCCA::TwoBodyIntCCA(Integral* integral,
                             const Ref<GaussianBasisSet> &bs1,
                             const Ref<GaussianBasisSet> &bs2,
			     const Ref<GaussianBasisSet> &bs3,
			     const Ref<GaussianBasisSet> &bs4,
			     size_t storage,
			     IntegralEvaluatorFactory eval_factory, 
                             bool use_opaque, string eval_type) :
  TwoBodyInt(integral,bs1,bs2,bs3,bs4)
{
  int2ecca_ = new Int2eCCA(integral,bs1,bs2,bs3,bs4,0,storage,
                           eval_factory,use_opaque,eval_type);
  buffer_ = int2ecca_->buffer();
  int2ecca_->set_redundant(redundant_);
  int_bound_min_ = SCHAR_MIN;
  tol_ = pow(2.0,double(int_bound_min_));
  loginv_ = 1.0/log(2.0);
}

void
TwoBodyIntCCA::compute_shell(int is, int js, int ks, int ls)
{
  int2ecca_->compute_erep(is,js,ks,ls);
}

int
TwoBodyIntCCA::log2_shell_bound(int is, int js, int ks, int ls)
{
  double value = int2ecca_->compute_bounds(is,js,ks,ls);
  int upper;
  if (value > tol_) {
    double log_dbl = log(value)*loginv_;
    upper = static_cast<int>( ceil(log_dbl) );
    // nasty check for rounding effects
    if( upper > SCHAR_MIN ) {
      int lower = upper - 1;
      double lower_dbl = static_cast<double>( lower ); 
      if( abs(log_dbl - lower_dbl) < 0.01 )
        upper = lower;
    }
  }
  else upper = int_bound_min_;
  return upper;
}

void
TwoBodyIntCCA::set_integral_storage(size_t storage)
{
//  throw FeatureNotImplemented("set_integral_storage needs to be implemented",
//                              __FILE__,__LINE__);
}


////////////////////////////////////////////////////////////////////////////
// TwoBodyDerivIntCCA

TwoBodyDerivIntCCA::TwoBodyDerivIntCCA(Integral* integral,
                             const Ref<GaussianBasisSet> &bs1,
                             const Ref<GaussianBasisSet> &bs2,
                             const Ref<GaussianBasisSet> &bs3,
                             const Ref<GaussianBasisSet> &bs4,
                             size_t storage,
                             IntegralEvaluatorFactory eval_factory,
                             bool use_opaque, string eval_type) :
  TwoBodyDerivInt(integral,bs1,bs2,bs3,bs4)
{
  int2ecca_ = new Int2eCCA(integral,bs1,bs2,bs3,bs4,1,storage,
                           eval_factory,use_opaque,eval_type);
  buffer_ = int2ecca_->buffer();
  int2ecca_->set_redundant(0);
  int_bound_min_ = SCHAR_MIN;
  tol_ = pow(2.0,double(int_bound_min_));
  loginv_ = 1.0/log(2.0);
}

void
TwoBodyDerivIntCCA::compute_shell(int is, int js, int ks, int ls,
                                  DerivCenters &dc )
{
  Chemistry::QC::GaussianBasis::DerivCenters cca_dc;
  cca_dc = Chemistry_QC_GaussianBasis_DerivCenters::_create();
  for( int id=0; id<cca_dc.n(); ++id ) {
    if( id == cca_dc.omitted_center() )
       dc.add_omitted(cca_dc.center(id),cca_dc.atom(id));
     else
       dc.add_center(cca_dc.center(id),cca_dc.atom(id));
  }

  int2ecca_->compute_erep_1der(is,js,ks,ls,cca_dc);

}

int
TwoBodyDerivIntCCA::log2_shell_bound(int is, int js, int ks, int ls)
{
  double value = int2ecca_->compute_bounds(is,js,ks,ls);
  int upper;
  if (value > tol_) {
    double log_dbl = log(value)*loginv_;
    upper = static_cast<int>( ceil(log_dbl) );
    // nasty check for rounding effects
    if( upper > SCHAR_MIN ) {
      int lower = upper - 1;
      double lower_dbl = static_cast<double>( lower );
      if( abs(log_dbl - lower_dbl) < 0.01 )
        upper = lower;
    }
  }
  else upper = int_bound_min_;
  return upper;
}

/////////////////////////////////////////////////////////////////////////////
// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:

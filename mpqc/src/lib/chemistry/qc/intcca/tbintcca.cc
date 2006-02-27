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
  cca_dc_ = Chemistry_QC_GaussianBasis_DerivCenters::_create();
  int2ecca_ = new Int2eCCA(integral,bs1,bs2,bs3,bs4,0,storage,
                           eval_factory,use_opaque,eval_type,cca_dc_);
  buffer_ = int2ecca_->buffer();
  int2ecca_->set_redundant(redundant_);
  int_bound_min_ = SCHAR_MIN;
  tol_ = pow(2.0,double(int_bound_min_));
  loginv_ = 1.0/log(2.0);
  /* need a way to query this through CCA interfaces */
  n_types_ = 1;
}

void
TwoBodyIntCCA::compute_shell(int is, int js, int ks, int ls)
{
  int2ecca_->compute_erep(is,js,ks,ls);
}

int
TwoBodyIntCCA::log2_shell_bound(int is, int js, int ks, int ls)
{
/*
  double value = int2ecca_->compute_bounds(is,js,ks,ls);
  int upper;
  if (value > tol_) {
    double log_dbl = log(value)*loginv_;
    upper = static_cast<int>( ceil(log_dbl) );
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

  return static_cast<int>( logb( int2ecca_->compute_bounds(is,js,ks,ls) ) );
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
  cca_dc_ = Chemistry_QC_GaussianBasis_DerivCenters::_create();
  int2ecca_ = new Int2eCCA(integral,bs1,bs2,bs3,bs4,1,storage,
                           eval_factory,use_opaque,eval_type,cca_dc_);
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
  int2ecca_->compute_erep_1der(is,js,ks,ls);

  dc.clear();
  if( cca_dc_.has_omitted_center() )
    dc.add_omitted(cca_dc_.omitted_center(),cca_dc_.omitted_atom());
  for( int i=0; i<cca_dc_.n(); ++i) 
    dc.add_center(cca_dc_.center(i),cca_dc_.atom(i));

/* 
  for( int i=0; i<cca_dc.n(); ++i)
     std::cerr << "intcca: cca_dc: n " << i << " center " << cca_dc.center(i) << " atom " << cca_dc.atom(i) << std
::endl;
  if( cca_dc.has_omitted_center() ) {
    std::cerr << "intcca: cca_dc: omitted center is " << cca_dc.omitted_center() << std::endl;
    std::cerr << "intcca: cca_dc: omitted atom is " << cca_dc.omitted_atom() << std::endl;
  }

  for( int i=0; i<dc.n(); ++i)
     std::cerr << "intcca: dc: n " << i << " center " << dc.center(i) << " atom " << dc.atom(i) << std::endl;
  if( dc.has_omitted_center() ) {
    std::cerr << "intcca: dc: omitted center is " << dc.omitted_center() << std::endl;
    std::cerr << "intcca: dc: omitted atom is " << dc.omitted_atom() << std::endl;
  } 

    if( is==8 && js==0 && ks==0 && ls==0 ) {
      std::cerr << "cca_dc.n(): " << cca_dc_.n() << std::endl;
      GaussianShell* s1 = &( bs1_->shell(is) );
      GaussianShell* s2 = &( bs2_->shell(js) );
      GaussianShell* s3 = &( bs3_->shell(ks) );
      GaussianShell* s4 = &( bs4_->shell(ls) );
      int nfunc = s1->nfunction() * s2->nfunction() * s3->nfunction() * s4->nfunction();
      std::cerr << "\nintcca: computing shell " << is << " " <<  js << " " << ks << " " << ls << std::endl;
      std::cerr << "intcca buffer for shell quartet:\n";
      int nc1 = s1->ncontraction();
      int nc2 = s2->ncontraction();
      int nc3 = s3->ncontraction();
      int nc4 = s4->ncontraction();
      std::cerr << "shellnum1: " << is << std::endl;
      for (int i=0; i<nc1; ++i)
        std::cerr << "am: " << s1->am(i) << std::endl;
      std::cerr << "shellnum2: " << js << std::endl;
      for (int i=0; i<nc2; ++i)
        std::cerr << "am: " << s2->am(i) << std::endl;
      std::cerr << "shellnum3: " << ks << std::endl;
      for (int i=0; i<nc3; ++i)
        std::cerr << "am: " << s3->am(i) << std::endl;
      std::cerr << "shellnum4: " << ls << std::endl;
      for (int i=0; i<nc4; ++i)
        std::cerr << "am: " << s4->am(i) << std::endl;

      for( int i=0; i<dc.n(); ++i) {
        std::cerr << "n " << i << " center " << dc.center(i) << " atom " << dc.atom(i) << std::endl;
        std::cerr << "  dx\n";
        for( int j=0; j<nfunc; ++j)
          std::cerr << "  " << buffer_[i*nfunc*3+j] << std::endl;
        std::cerr << "  dy\n";
        for( int j=nfunc; j<nfunc*2; ++j)
          std::cerr << "  " << buffer_[i*nfunc*3+j] << std::endl;
        std::cerr << "  dz\n";
        for( int j=nfunc*2; j<nfunc*3; ++j)
          std::cerr << "  " << buffer_[i*nfunc*3+j] << std::endl;
      }
    }
*/

}

int
TwoBodyDerivIntCCA::log2_shell_bound(int is, int js, int ks, int ls)
{
/*
  double value = int2ecca_->compute_bounds_1der(is,js,ks,ls);
  int upper;
  if (value > tol_) {
    double log_dbl = log(value)*loginv_;
    upper = static_cast<int>( ceil(log_dbl) );
    if( upper > SCHAR_MIN ) {
      int lower = upper - 1;
      double lower_dbl = static_cast<double>( lower );
      if( abs(log_dbl - lower_dbl) < 0.01 )
        upper = lower;
    }
  }
  else upper = int_bound_min_;
  return upper;
*/

  return static_cast<int>( logb( int2ecca_->compute_bounds_1der(is,js,ks,ls) ) );
}

/////////////////////////////////////////////////////////////////////////////
// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:

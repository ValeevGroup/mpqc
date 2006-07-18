//
// obintcca.cc
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

#include <chemistry/cca/int/obintcca.h>
#include <util/class/scexception.h>

using namespace std;
using namespace sc;
using namespace Chemistry;
using namespace Chemistry::QC::GaussianBasis;

////////////////////////////////////////////////////////////////////////////
// OneBodyIntCCA

OneBodyIntCCA::OneBodyIntCCA( Integral* integral,
			      const Ref<GaussianBasisSet>& bs1, 
			      const Ref<GaussianBasisSet>& bs2,
			      IntegralSuperFactory fac,
			      CompositeIntegralDescr cdesc,
			      bool  use_opaque,
                              bool reorder
                              ):
  OneBodyInt(integral,bs1,bs2), bs1_(bs1), bs2_(bs2),
  eval_factory_(fac), cdesc_(cdesc), use_opaque_(use_opaque),
  reorder_(reorder)
{
  // create cca basis sets
  cca_bs1_ = MPQC::GaussianBasis_Molecular::_create();
  cca_bs1_.initialize( bs1_.pointer(), bs1_->label() );
  if( bs1_.pointer() != bs2_.pointer() ) {
    cca_bs2_ = MPQC::GaussianBasis_Molecular::_create();
    cca_bs2_.initialize( bs2_.pointer(), bs2_->label() );
  }
  else
    cca_bs2_ = cca_bs1_;

  // there are no onebody evaluators currently in mpqc
  // that handle multiple types, so CompositeDescr contains exactly 1 Descr

  IntegralDescr desc = cdesc_.get_descr(0);
  n_segment_ = desc.get_n_segment();
  type_ = desc.get_type();
  int scratchsize = bs1_->max_ncartesian_in_shell()
                      *  bs2_->max_ncartesian_in_shell()
                      *  n_segment_;

  temp_buffer_ = new double[scratchsize];

  eval_ = eval_factory_.get_evaluator2( cdesc_, cca_bs1_, cca_bs2_ );
  if( use_opaque_ )
    buffer_ = static_cast<double*>( eval_.get_buffer(desc) );
  else
    buffer_ = new double[scratchsize];
}

OneBodyIntCCA::~OneBodyIntCCA()
{
}

void
OneBodyIntCCA::compute_shell(int i, int j)
{
  if( use_opaque_ )
    eval_.compute( i, j );
  else {
    sidl_buffer_ = eval_.compute_array( type_, i, j );
    int sidl_size = 1 + sidl_buffer_.upper(0) - sidl_buffer_.lower(0);
    for(int ii=0; ii<sidl_size; ++ii)
      buffer_[ii] = sidl_buffer_.get(ii);
  }

  // reorder for mpqc's wacky 1-body multi-segment ordering
  if( n_segment_ > 1 && reorder_ ) {
    GaussianShell* s1 = &( bs1_->shell(i) );
    GaussianShell* s2 = &( bs2_->shell(j) );
    int nfunc = s1->nfunction() * s2->nfunction();
    for(int i=0; i<nfunc*n_segment_; ++i)
      temp_buffer_[i] = buffer_[i];
    for(int i=0; i<nfunc; ++i)
      for( int si=0; si<n_segment_; ++si)
        buffer_[i*n_segment_+si] = temp_buffer_[si*nfunc+i];
  }

  // temporary debugging stuff for cca integrals comparison
/*
  if( 1 ) {
    std::cerr << "CCA buffer for shell doublet:\n";
    std::cerr << "shellnum1: " << i << std::endl;
    GaussianShell* s1 = &( bs1_->shell(i) );
    int nc1 = s1->ncontraction();
    for (int ii=0; ii<nc1; ++ii)
      std::cerr << "am: " << s1->am(ii) << std::endl;
    std::cerr << "shellnum2: " << j << std::endl;
    GaussianShell* s2 = &( bs2_->shell(j) );
    int nc2 = s2->ncontraction();
    for (int ii=0; ii<nc2; ++ii)
      std::cerr << "am: " << s2->am(ii) << std::endl;

    int nfunc = s1->max_cartesian() * s2->max_cartesian();
    for( int ii=0; ii<nfunc; ++ii)
      std::cerr << buffer_[ii] << std::endl;
    std::cerr << std::endl;
  }
*/

}

bool
OneBodyIntCCA::cloneable()
{
  return true;
}

Ref<OneBodyInt>
OneBodyIntCCA::clone()
{
  return new OneBodyIntCCA( integral_, bs1_, bs2_, 
                            eval_factory_, cdesc_, use_opaque_, reorder_ );
}

// ////////////////////////////////////////////////////////////////////////////
// // OneBodyDerivIntCCA

OneBodyDerivIntCCA::OneBodyDerivIntCCA(Integral *integral,
                                       const Ref<GaussianBasisSet>&bs1,
                                       const Ref<GaussianBasisSet>&bs2,
                                       IntegralSuperFactory eval_factory,
				       CompositeIntegralDescr cdesc,
                                       bool use_opaque,
                                       bool reorder
                                       ):
  OneBodyDerivInt(integral,bs1,bs2), bs1_(bs1), bs2_(bs2),
  eval_factory_(eval_factory), cdesc_(cdesc), use_opaque_(use_opaque),
  reorder_(reorder)
{
  // create cca basis sets
  cca_bs1_ = MPQC::GaussianBasis_Molecular::_create();
  cca_bs1_.initialize( bs1_.pointer(), bs1_->label() );
  if( bs1_.pointer() != bs2_.pointer() ) {
    cca_bs2_ = MPQC::GaussianBasis_Molecular::_create();
    cca_bs2_.initialize( bs2_.pointer(), bs2_->label() );
  }
  else
    cca_bs2_ = cca_bs1_;

  // there are no onebody evaluators currently in mpqc
  // that handle multiple types, so CompositeDescr contains exactly 1 Descr

  IntegralDescr desc = cdesc_.get_descr(0);
  n_segment_ = desc.get_n_segment();
  int scratchsize = bs1_->max_ncartesian_in_shell()
                      *  bs2_->max_ncartesian_in_shell()
                      *  n_segment_ * 3;

  temp_buffer_ = new double[scratchsize];

  cca_dc_ = desc.get_deriv_centers();
  eval_ = eval_factory_.get_evaluator2( cdesc_, cca_bs1_, cca_bs2_ );
  if( use_opaque_ )
    buffer_ = static_cast<double*>( eval_.get_buffer( desc ) );
  else
    buffer_ = new double[scratchsize];
}

OneBodyDerivIntCCA::~OneBodyDerivIntCCA()
{
}

void
OneBodyDerivIntCCA::compute_shell(int i, int j, DerivCenters& c)
{
 throw SCException("I thought this was never called for one-body evals",
                   __FILE__,__LINE__);
}

void 
OneBodyDerivIntCCA::compute_shell(int i, int j, int atom) 
{
  cca_dc_.set_deriv_atom( atom );
  eval_.compute(i,j);

  // reorder for mpqc's wacky 1-body derivative ordering
  if( n_segment_ > 1 ) 
    throw SCException("Deritative, multi-segment buffer structure unclear",
                      __FILE__,__LINE__);

  if( reorder_ ) {
    GaussianShell* s1 = &( bs1_->shell(i) );
    GaussianShell* s2 = &( bs2_->shell(j) );
    int nfunc = s1->nfunction() * s2->nfunction();

    for(int i=0; i<nfunc*3; ++i)
      temp_buffer_[i] = buffer_[i];
    for(int i=0; i<nfunc; ++i)
      for( int di=0; di<3; ++di)
        buffer_[i*3+di] = temp_buffer_[di*nfunc+i];
  }

// temporary debuging stuff for cca integrals comparison
/*
  if( 1 ) {
    std::cerr << "\ncca buffer:\n";
    int nc1 = s1->ncontraction();
    int nc2 = s2->ncontraction();

    int nfunc = s1->max_cartesian() * s2->max_cartesian();
    std::cerr << " dx: " << std::endl;
    for( int i=0; i<nfunc; ++i)
      std::cerr << "     " << buffer_[i] << std::endl;
    std::cerr << " dy: " << std::endl;
    for( int i=0; i<nfunc; ++i)
      std::cerr << "     " << buffer_[nfunc+i] << std::endl;
    std::cerr << " dz: " << std::endl;
    for( int i=0; i<nfunc; ++i)
      std::cerr << "     " << buffer_[2*nfunc+i] << std::endl;
  }
*/

}

bool
OneBodyDerivIntCCA::cloneable()
{
  return true;
}

Ref<OneBodyDerivInt>
OneBodyDerivIntCCA::clone()
{
  return new OneBodyDerivIntCCA( integral_, bs1_, bs2_,
                                 eval_factory_, cdesc_, 
                                 use_opaque_, reorder_ );
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:

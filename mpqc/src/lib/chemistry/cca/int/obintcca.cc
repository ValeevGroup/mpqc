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

////////////////////////////////////////////////////////////////////////////
// OneBodyIntCCA

OneBodyIntCCA::OneBodyIntCCA( Integral* integral,
			      const Ref<GaussianBasisSet>& bs1, 
			      const Ref<GaussianBasisSet>& bs2,
			      IntegralSuperFactory fac,
			      CompositeIntegralDescr cdesc,
			      vector<string> factories, 
			      bool  use_opaque ):
  OneBodyInt(integral,bs1,bs2), bs1_(bs1), bs2_(bs2),
  eval_factory_(fac), cdesc_(cdesc), factories_(factories), 
   use_opaque_(use_opaque) 
{
  
  /* The efield routines look like derivatives so nshell*3 */
  int scratchsize=0,nshell2;
  nshell2 = bs1_->max_ncartesian_in_shell()*bs2_->max_ncartesian_in_shell();
  scratchsize = nshell2*3;
  if( !use_opaque_ )
    buffer_ = new double[scratchsize];

  // create cca basis sets
  cca_bs1_ = MPQC::GaussianBasis_Molecular::_create();
  cca_bs1_.initialize( bs1_.pointer(), bs1_->name() );
  if( bs1_.pointer() != bs2_.pointer() ) {
    cca_bs2_ = MPQC::GaussianBasis_Molecular::_create();
    cca_bs2_.initialize( bs2_.pointer(), bs2_->name() );
  }
  else
    cca_bs2_ = cca_bs1_;

  // set factory config, CompositeDescr should contain exactly 1 Descr
  sidl::array<string> sidl_factories = sidl::array<string>::create1d(1);
  sidl_factories.set( 0, factories_[0] );
  eval_factory_.set_source_factories( sidl_factories );

  eval_ = eval_factory_.get_evaluator2( cdesc_, cca_bs1_, cca_bs2_ );
  buffer_ = static_cast<double*>( eval_.get_buffer( cdesc_.get_descr(0) ) );
  std::cerr << "IntCCA set buffer to " << buffer_ << std::endl;
  
}

OneBodyIntCCA::~OneBodyIntCCA()
{
}

void
OneBodyIntCCA::compute_shell(int i, int j)
{
  eval_.compute( i, j );
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
                            eval_factory_, cdesc_, factories_, use_opaque_ );
}

// ////////////////////////////////////////////////////////////////////////////
// // OneBodyDerivIntCCA

OneBodyDerivIntCCA::OneBodyDerivIntCCA(Integral *integral,
                                       const Ref<GaussianBasisSet>&bs1,
                                       const Ref<GaussianBasisSet>&bs2,
                                       IntegralSuperFactory eval_factory,
				       CompositeIntegralDescr cdesc,
				       vector<string> factories,
                                       bool use_opaque ):
  OneBodyDerivInt(integral,bs1,bs2), bs1_(bs1), bs2_(bs2),
  eval_factory_(eval_factory), cdesc_(cdesc), factories_(factories),
  use_opaque_(use_opaque)
{
}

OneBodyDerivIntCCA::~OneBodyDerivIntCCA()
{
}

void
OneBodyDerivIntCCA::compute_shell(int i, int j, DerivCenters& c)
{
}

void 
OneBodyDerivIntCCA::compute_shell(int i, int j, int c) 
{
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:

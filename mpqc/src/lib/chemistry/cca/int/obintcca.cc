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
			      bool  use_opaque ):
  OneBodyInt(integral,bs1,bs2), bs1_(bs1), bs2_(bs2),
  eval_factory_(fac), cdesc_(cdesc), use_opaque_(use_opaque) 
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

  // set factory config, there are no onebody evaluators currently in mpqc
  // that handle multiple types, so CompositeDescr contains exactly 1 Descr
/*
  sidl::array<string> sidl_factories = sidl::array<string>::create1d(1);
  sidl_factories.set( 0, factories_[0] );
  eval_factory_.set_source_factories( sidl_factories );
*/

  eval_ = eval_factory_.get_evaluator2( cdesc_, cca_bs1_, cca_bs2_ );
  buffer_ = static_cast<double*>( eval_.get_buffer( cdesc_.get_descr(0) ) );
}

OneBodyIntCCA::~OneBodyIntCCA()
{
}

void
OneBodyIntCCA::compute_shell(int i, int j)
{
  eval_.compute( i, j );

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
                            eval_factory_, cdesc_, use_opaque_ );
}

// ////////////////////////////////////////////////////////////////////////////
// // OneBodyDerivIntCCA

OneBodyDerivIntCCA::OneBodyDerivIntCCA(Integral *integral,
                                       const Ref<GaussianBasisSet>&bs1,
                                       const Ref<GaussianBasisSet>&bs2,
                                       IntegralSuperFactory eval_factory,
				       CompositeIntegralDescr cdesc,
                                       bool use_opaque ):
  OneBodyDerivInt(integral,bs1,bs2), bs1_(bs1), bs2_(bs2),
  eval_factory_(eval_factory), cdesc_(cdesc), use_opaque_(use_opaque)
{
  int scratchsize=0,nshell2;
  nshell2 = bs1_->max_ncartesian_in_shell()*bs2_->max_ncartesian_in_shell();
  scratchsize = nshell2*3;
  if( !use_opaque_ )
    buffer_ = new double[scratchsize];
  temp_buffer_ = new double[scratchsize];

  // create cca basis sets
  cca_bs1_ = MPQC::GaussianBasis_Molecular::_create();
  cca_bs1_.initialize( bs1_.pointer(), bs1_->name() );
  if( bs1_.pointer() != bs2_.pointer() ) {
    cca_bs2_ = MPQC::GaussianBasis_Molecular::_create();
    cca_bs2_.initialize( bs2_.pointer(), bs2_->name() );
  }
  else
    cca_bs2_ = cca_bs1_;

  // set factory config, there are no onebody evaluators currently in mpqc
  // that handle multiple types, so CompositeDescr contains exactly 1 Descr
/*
  sidl::array<string> sidl_factories = sidl::array<string>::create1d(1);
  sidl_factories.set( 0, factories_[0] );
  eval_factory_.set_source_factories( sidl_factories );
*/

  IntegralDescr idesc = cdesc_.get_descr(0);
  cca_dc_ = idesc.get_deriv_centers();

  eval_ = eval_factory_.get_evaluator2( cdesc_, cca_bs1_, cca_bs2_ );
  buffer_ = static_cast<double*>( eval_.get_buffer( idesc ) );
}

OneBodyDerivIntCCA::~OneBodyDerivIntCCA()
{
}

void
OneBodyDerivIntCCA::compute_shell(int i, int j, DerivCenters& c)
{

 throw SCException("I thought this was never called for one-body evals",
                   __FILE__,__LINE__);
/*
  c.clear();
  c.add_center(0,basis1(),i);
  c.add_omitted(1,basis2(),j);
  for( int id=0; id<c.n(); ++id ) {
     if( id == c.omitted_center() )
       cca_dc_.add_omitted(c.center(id),c.atom(id));
     else
       cca_dc_.add_center(c.center(id),c.atom(id));
  }
*/
}

void 
OneBodyDerivIntCCA::compute_shell(int i, int j, int atom) 
{
  cca_dc_.set_deriv_atom( atom );
  eval_.compute(i,j);

  GaussianShell* s1 = &( bs1_->shell(i) );
  GaussianShell* s2 = &( bs2_->shell(j) );
  int nfunc = s1->nfunction() * s2->nfunction();


  // reorder for mpqc's wacky 1-body derivative ordering
  for(int i=0; i<nfunc*3; ++i)
    temp_buffer_[i] = buffer_[i];
  for(int i=0; i<nfunc; ++i)
    for( int di=0; di<3; ++di)
      buffer_[i*3+di] = temp_buffer_[di*nfunc+i];

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
                                 eval_factory_, cdesc_, use_opaque_ );
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:

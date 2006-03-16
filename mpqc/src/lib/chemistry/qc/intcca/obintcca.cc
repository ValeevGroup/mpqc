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

#include <chemistry/qc/intcca/obintcca.h>
#include <util/class/scexception.h>
#include <Chemistry_Chemistry_QC_GaussianBasis_DerivCenters.hh>

using namespace std;
using namespace Chemistry;
using namespace Chemistry::QC::GaussianBasis;
using namespace sc;

////////////////////////////////////////////////////////////////////////////
// OneBodyIntCCA

OneBodyIntCCA::OneBodyIntCCA(Integral* integral,
                           const Ref<GaussianBasisSet>&bs1,
                           const Ref<GaussianBasisSet>&bs2,
			   IntegralEvaluatorFactory eval_factory,
                           IntegralFunction ifunc,
                           bool use_opaque ):
  OneBodyInt(integral,bs1,bs2), intfunc_(ifunc), eval_factory_(eval_factory), 
  use_opaque_(use_opaque) 
{
  ObIntEvalType int_type;
  if( ifunc == &Int1eCCA::overlap ) int_type = ObIntEvalType_OVERLAP;
  else if (ifunc == &Int1eCCA::kinetic ) int_type = ObIntEvalType_KINETIC;
  else if (ifunc == &Int1eCCA::nuclear ) int_type = ObIntEvalType_NUCLEAR;
  else if (ifunc == &Int1eCCA::hcore ) int_type = ObIntEvalType_OEHAM;
  cca_dc_ = Chemistry_QC_GaussianBasis_DerivCenters::_create();
  int1ecca_ = new Int1eCCA(integral,bs1,bs2,0,eval_factory,
                           int_type,use_opaque,cca_dc_);
  buffer_ = int1ecca_->buffer();
}

OneBodyIntCCA::~OneBodyIntCCA()
{
}

void
OneBodyIntCCA::compute_shell(int i, int j)
{
  (int1ecca_.pointer()->*intfunc_)(i, j);
}

bool
OneBodyIntCCA::cloneable()
{
  return true;
}

Ref<OneBodyInt>
OneBodyIntCCA::clone()
{
  return new OneBodyIntCCA(integral_, bs1_, bs2_, 
                           eval_factory_, intfunc_, use_opaque_ );
}

// ////////////////////////////////////////////////////////////////////////////
// // OneBodyDerivIntCCA

OneBodyDerivIntCCA::OneBodyDerivIntCCA(Integral *integral,
                                       const Ref<GaussianBasisSet>&bs1,
                                       const Ref<GaussianBasisSet>&bs2,
                                       IntegralEvaluatorFactory eval_factory,
                                       bool use_opaque,
                                       ObIntEvalType int_type ):
  OneBodyDerivInt(integral,bs1,bs2), eval_factory_(eval_factory),
  use_opaque_(use_opaque), eval_type_(int_type)
{
  cca_dc_ = Chemistry_QC_GaussianBasis_DerivCenters::_create();
  int1ecca_ = new Int1eCCA(integral,bs1,bs2,1,eval_factory,
                           eval_type_,use_opaque,cca_dc_);
  buffer_ = int1ecca_->buffer();
}

OneBodyDerivIntCCA::~OneBodyDerivIntCCA()
{
}

void
OneBodyDerivIntCCA::compute_shell(int i, int j, DerivCenters& c)
{

  std::cerr << "THIS IS NEVER USED\n";

  cca_dc_.clear();
  
  if( eval_type_ == ObIntEvalType_OVERLAP )
    int1ecca_->overlap_1der(i,j);
  else if( eval_type_ == ObIntEvalType_KINETIC )  
    int1ecca_->kinetic_1der(i,j);
  else if( eval_type_ == ObIntEvalType_NUCLEAR )
    int1ecca_->nuclear_1der(i,j);
  else if( eval_type_ == ObIntEvalType_OEHAM )
    int1ecca_->hcore_1der(i,j);

  c.clear();
  for( int id=0; id<cca_dc_.n(); ++id ) {
    if( id == cca_dc_.omitted_center() )
      c.add_omitted(cca_dc_.center(id),cca_dc_.atom(id));
     else
       c.add_center(cca_dc_.center(id),cca_dc_.atom(id));
  }

}

void 
OneBodyDerivIntCCA::compute_shell(int i, int j, int c) {

  if( eval_type_ == ObIntEvalType_OVERLAP )
    int1ecca_->overlap_1der(i,j,c);
  else if( eval_type_ == ObIntEvalType_KINETIC )
    int1ecca_->kinetic_1der(i,j,c);
  else if( eval_type_ == ObIntEvalType_NUCLEAR )
    int1ecca_->nuclear_1der(i,j,c);
  else if( eval_type_ == ObIntEvalType_OEHAM )
    int1ecca_->hcore_1der(i,j,c);
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:

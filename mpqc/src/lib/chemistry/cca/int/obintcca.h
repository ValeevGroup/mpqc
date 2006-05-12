//
// obintcca.h
//
// Copyright (C) 2004, Sandia National Laboratories
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
#pragma interface
#endif

#ifndef _chemistry_cca_int_obintcca_h
#define _chemistry_cca_int_obintcca_h

#include <vector>
#include <chemistry/qc/basis/obint.h>
#include <Chemistry_QC_GaussianBasis_IntegralSuperFactory.hh>
#include <Chemistry_QC_GaussianBasis_CompositeIntegralDescr.hh>
#include <Chemistry_QC_GaussianBasis_DerivCenters.hh>
#include <MPQC_GaussianBasis_Molecular.hh>

using namespace std;
using namespace sc;
using namespace Chemistry;
using namespace Chemistry::QC::GaussianBasis;

namespace sc {

// /////////////////////////////////////////////////////////////////////////

/** This implements one body integrals through the CCA interface. */

  class OneBodyIntCCA : public OneBodyInt {
    
  private:
    Integral* integral_;
    Ref<GaussianBasisSet> bs1_, bs2_;
    IntegralSuperFactory eval_factory_;
    CompositeIntegralDescr cdesc_;
    vector<string> factories_;
    bool use_opaque_;
    MPQC::GaussianBasis_Molecular cca_bs1_, cca_bs2_;
    double* buff_;
    IntegralEvaluator2 eval_;

  protected:
    Chemistry::QC::GaussianBasis::DerivCenters cca_dc_;
    
  public:
    OneBodyIntCCA( Integral* integral,
		   const Ref<GaussianBasisSet>&, 
		   const Ref<GaussianBasisSet>&,
		   IntegralSuperFactory,
		   CompositeIntegralDescr,
                   vector<string> factories,
		   bool );
    ~OneBodyIntCCA();
    void compute_shell(int,int);
    bool cloneable();
    Ref<OneBodyInt> clone();
};

///////////////////////////////////////////////////////////////////////////////

 /** This implements one body derivative integrals through the
     CCA interface */

  class OneBodyDerivIntCCA : public OneBodyDerivInt {
    
  private:
    Integral* integral_;
    Ref<GaussianBasisSet> bs1_, bs2_;
    IntegralSuperFactory eval_factory_;
    CompositeIntegralDescr cdesc_;
    vector<string> factories_;
    bool use_opaque_;
    MPQC::GaussianBasis_Molecular cca_bs1_, cca_bs2_;
    double* buff_;
    IntegralEvaluator2 eval_;

  protected:
    Chemistry::QC::GaussianBasis::DerivCenters cca_dc_;
 
  public:
    OneBodyDerivIntCCA( Integral* integral,
		        const Ref<GaussianBasisSet>&,
		        const Ref<GaussianBasisSet>&,
			IntegralSuperFactory,
		        CompositeIntegralDescr,
                        vector<string> factories,
			bool );
    ~OneBodyDerivIntCCA();
    void compute_shell(int, int, DerivCenters&);
    void compute_shell(int, int, int);
    bool cloneable();
    Ref<OneBodyDerivInt> clone();
  };

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:

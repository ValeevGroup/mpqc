//
// tbintcca.h
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

#ifndef _chemistry_qc_intcca_tbintcca_h
#define _chemistry_qc_intcca_tbintcca_h

#include <vector>
#include <chemistry/qc/basis/tbint.h>
#include <Chemistry_QC_GaussianBasis_IntegralSuperFactory.hh>
#include <Chemistry_QC_GaussianBasis_DerivCenters.hh>
#include <MPQC_GaussianBasis_Molecular.hh>

using namespace std;
using namespace sc;
using namespace Chemistry;
using namespace Chemistry::QC::GaussianBasis;

namespace sc {

  //////////////////////////////////////////////////////////////////////////

/** This implements two body integrals through the CCA interface. */

  class TwoBodyIntCCA : public TwoBodyInt {
  
  private:

    Integral* integral_;
    Ref<GaussianBasisSet> bs1_, bs2_, bs3_, bs4_;
    IntegralSuperFactory eval_factory_;
    CompositeIntegralDescr cdesc_;
    vector<string> factories_;
    bool use_opaque_;
    MPQC::GaussianBasis_Molecular cca_bs1_, cca_bs2_, cca_bs3_, cca_bs4_;
    double* buff_;
    IntegralEvaluator4 eval_;
    bool redundant_;
    int int_bound_min_;
    double tol_;
    double loginv_;

  public:
    TwoBodyIntCCA( Integral* integral,
		   const Ref<GaussianBasisSet>&b1, 
		   const Ref<GaussianBasisSet>&b2,
		   const Ref<GaussianBasisSet>&b3, 
		   const Ref<GaussianBasisSet>&b4,	   
		   size_t,
		   IntegralSuperFactory,
		   CompositeIntegralDescr,
		   vector<string> factories,
		   bool );
    ~TwoBodyIntCCA();

    void compute_shell(int,int,int,int);
    int log2_shell_bound(int,int,int,int);
    int redundant() const { return redundant_; }
    void set_redundant(int i) { redundant_ = i; }
    unsigned int num_tbint_types() const;
    void remove_redundant(int,int,int,int);
};


  ////////////////////////////////////////////////////////////////////////////

  /** This implements two body derivative integrals through the 
      CCA interface. */
  class TwoBodyDerivIntCCA : public TwoBodyDerivInt {

  private:
    Integral* integral_;
    Ref<GaussianBasisSet> bs1_, bs2_, bs3_, bs4_;
    IntegralSuperFactory eval_factory_;
    CompositeIntegralDescr cdesc_;
    vector<string> factories_;
    bool use_opaque_;
    MPQC::GaussianBasis_Molecular cca_bs1_, cca_bs2_, cca_bs3_, cca_bs4_;
    double* buff_;
    IntegralEvaluator4 eval_;
    bool redundant_;
    double tol_;
    double loginv_;
    int int_bound_min_;
    int max_deriv_lvl_;
    Chemistry::QC::GaussianBasis::DerivCenters cca_dc_;
    sidl::array<int> sidl_buffer_;

  public:
    TwoBodyDerivIntCCA( Integral* integral,
			const Ref<GaussianBasisSet>&b1, 
			const Ref<GaussianBasisSet>&b2,
			const Ref<GaussianBasisSet>&b3, 
			const Ref<GaussianBasisSet>&b4,	   
			size_t,
			IntegralSuperFactory,
			CompositeIntegralDescr,
			vector<string> factories,
			bool );
    ~TwoBodyDerivIntCCA();

    void compute_shell(int,int,int,int,DerivCenters&);
    int log2_shell_bound(int,int,int,int);
    int redundant() const { return redundant_; }
    void set_redundant(int i) { redundant_ = i; }
    unsigned int num_tbint_types() const;
    void copy_buffer(int);
  };
  
}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:

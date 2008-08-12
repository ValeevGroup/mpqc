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
#include <sidl_cxx.hxx>
#include <chemistry/qc/basis/obint.h>
#include <Chemistry_QC_GaussianBasis_IntegralEvaluatorFactoryInterface.hxx>
#include <Chemistry_QC_GaussianBasis_CompositeDescrInterface.hxx>
#include <Chemistry_QC_GaussianBasis_DerivCentersInterface.hxx>
#include <MPQC_GaussianBasisMolecular.hxx>

namespace sc {

// /////////////////////////////////////////////////////////////////////////

/** This implements one body integrals through the CCA interface. */

  class OneBodyIntCCA : public OneBodyInt {
    
  private:
    Integral* integral_;
    Ref<GaussianBasisSet> bs1_, bs2_;
    Chemistry::QC::GaussianBasis::IntegralEvaluatorFactoryInterface eval_factory_;
    Chemistry::QC::GaussianBasis::CompositeDescrInterface cdesc_;
    bool reorder_;
    MPQC::GaussianBasisMolecular cca_bs1_, cca_bs2_;
    Chemistry::QC::GaussianBasis::IntegralEvaluator2Interface eval_;
    double* temp_buffer_;
    int n_segment_;
    std::string type_;
    sidl::array<double> sidl_buffer_;

  protected:
    
  public:
    OneBodyIntCCA( Integral* integral,
		   const Ref<GaussianBasisSet>&, 
		   const Ref<GaussianBasisSet>&,
		   Chemistry::QC::GaussianBasis::IntegralEvaluatorFactoryInterface,
		   Chemistry::QC::GaussianBasis::CompositeDescrInterface,
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
    Chemistry::QC::GaussianBasis::IntegralEvaluatorFactoryInterface eval_factory_;
    Chemistry::QC::GaussianBasis::CompositeDescrInterface cdesc_;
    bool reorder_;
    MPQC::GaussianBasisMolecular cca_bs1_, cca_bs2_;
    double* buff_;
    double* temp_buffer_;
    Chemistry::QC::GaussianBasis::IntegralEvaluator2Interface eval_;
    Chemistry::QC::GaussianBasis::DerivCentersInterface cca_dc_;
    int n_segment_;
    std::string type_;
    sidl::array<double> sidl_buffer_;

  public:
    OneBodyDerivIntCCA( Integral* integral,
		        const Ref<GaussianBasisSet>&,
		        const Ref<GaussianBasisSet>&,
			Chemistry::QC::GaussianBasis::IntegralEvaluatorFactoryInterface,
		        Chemistry::QC::GaussianBasis::CompositeDescrInterface,
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

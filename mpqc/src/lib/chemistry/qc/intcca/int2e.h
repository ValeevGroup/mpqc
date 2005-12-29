//
// int1e.h
//
// Copyright (C) 2004 Sandia National Laboratories.
//
// Author: Joseph Kenny <jpkenny@sandia.gov>
// Maintainer: JPK
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

#ifndef _chemistry_qc_intcca_int2e_h
#define _chemistry_qc_intcca_int2e_h

#include <sidl_cxx.hh>
#include <Chemistry_QC_GaussianBasis_IntegralEvaluatorFactory.hh>
#include <Chemistry_QC_GaussianBasis_IntegralEvaluator4.hh>
#include <Chemistry_QC_GaussianBasis_DerivCenters.hh>
#include <MPQC_GaussianBasis_Molecular.hh>
#include <chemistry/qc/basis/integral.h>

using namespace std;
using namespace MPQC;
using namespace Chemistry;
using namespace Chemistry::QC::GaussianBasis;

namespace sc {

class Integral;

/** Int2eCCA adapts CCA integrals components for use within SC.  It is
    used by TwoBodyIntCCA and TwoBodyDerivIntCCA to implement the 
    IntegralCCA class. */
class Int2eCCA: public RefCount {

  private:
    IntegralEvaluatorFactory eval_factory_;
    Ref<GaussianBasisSet> bs1_;
    Ref<GaussianBasisSet> bs2_;
    Ref<GaussianBasisSet> bs3_;
    Ref<GaussianBasisSet> bs4_;
    GaussianBasis_Molecular cca_bs1_;
    GaussianBasis_Molecular cca_bs2_;
    GaussianBasis_Molecular cca_bs3_;
    GaussianBasis_Molecular cca_bs4_;
    sidl::array<double> sidl_buffer_;
    double *buffer_;
    bool use_opaque_;
    void copy_buffer(int);
    IntegralEvaluator4 erep_;
    IntegralEvaluator4 erep_1der_;
    IntegralEvaluator4 *erep_ptr_;
    IntegralEvaluator4 *erep_1der_ptr_;
    int redundant_;
    void remove_redundant(int,int,int,int);

  protected:
    Integral *integral_;

  public:
    Int2eCCA(Integral *integral,
             const Ref<GaussianBasisSet>&b1,
             const Ref<GaussianBasisSet>&b2,
             const Ref<GaussianBasisSet>&b3,
             const Ref<GaussianBasisSet>&b4,
             int order, size_t storage, IntegralEvaluatorFactory, 
             bool, string );
    ~Int2eCCA() {};
    double *buffer() { return buffer_; }
    void compute_erep( int is, int js, int ks, int ls );
    void compute_erep_1der( int is, int js, int ks, int ls,
                            Chemistry::QC::GaussianBasis::DerivCenters &dc);
    int redundant() const { return redundant_; }
    void set_redundant(int i) { redundant_ = i; }

};

}

#endif

// Local Variables:
// mode: c++
// End:

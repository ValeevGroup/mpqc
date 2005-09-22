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

#ifndef _chemistry_qc_intcca_int1e_h
#define _chemistry_qc_intcca_int1e_h

#include <sidl_cxx.hh>
#include <Chemistry_QC_GaussianBasis_IntegralEvaluatorFactory.hh>
#include <Chemistry_QC_GaussianBasis_IntegralEvaluator2.hh>
#include <Chemistry_QC_GaussianBasis_Molecular.hh>
#include <Chemistry_Chemistry_QC_GaussianBasis_DerivCenters.hh>
#include <Chemistry_Chemistry_QC_GaussianBasis_DerivCenters.hh>
#include <MPQC_GaussianBasis_Molecular.hh>
#include <chemistry/qc/basis/integral.h>

using namespace std;
using namespace Chemistry;
using namespace Chemistry::QC::GaussianBasis;
using namespace MPQC;
namespace sc {

class Integral;

/** Int1eCCA adapts CCA integrals components for use within SC.  It is
    used by OneBodyIntCCA and OneBodyDerivIntCCA to implement the 
    IntegralCCA class. */
class Int1eCCA: public RefCount {

  private:
    IntegralEvaluatorFactory eval_factory_;
    Ref<GaussianBasisSet> bs1_;
    Ref<GaussianBasisSet> bs2_;
    GaussianBasis_Molecular cca_bs1_;
    GaussianBasis_Molecular cca_bs2_;
    sidl::array<double> sidl_buffer_;
    double *buff_;
    bool use_opaque_;
    void copy_buffer();
    IntegralEvaluator2 overlap_;
    IntegralEvaluator2 overlap_1der_;
    IntegralEvaluator2 kinetic_;
    IntegralEvaluator2 kinetic_1der_;
    IntegralEvaluator2 nuclear_;
    IntegralEvaluator2 nuclear_1der_;
    IntegralEvaluator2 hcore_;
    IntegralEvaluator2 hcore_1der_;
    IntegralEvaluator2 *overlap_ptr_;
    IntegralEvaluator2 *overlap_1der_ptr_;
    IntegralEvaluator2 *kinetic_ptr_;
    IntegralEvaluator2 *kinetic_1der_ptr_;
    IntegralEvaluator2 *nuclear_ptr_;
    IntegralEvaluator2 *nuclear_1der_ptr_;
    IntegralEvaluator2 *hcore_ptr_;
    IntegralEvaluator2 *hcore_1der_ptr_;
    Chemistry_QC_GaussianBasis_DerivCenters cca_dc_;

  protected:
    Integral *integral_;

  public:
    Int1eCCA(Integral *integral,
             const Ref<GaussianBasisSet>&b1,
	     const Ref<GaussianBasisSet>&b2,
	     int order, IntegralEvaluatorFactory, std::string, bool);
    ~Int1eCCA();

    double *buffer() { return buff_; }
    void overlap(int ish, int jsh);
    void overlap_1der(int ish, int jsh, 
                      Chemistry_QC_GaussianBasis_DerivCenters &dc);
    void overlap_1der(int ish, int jsh, int c);
    void kinetic(int ish, int jsh);
    void kinetic_1der(int ish, int jsh, 
                      Chemistry_QC_GaussianBasis_DerivCenters &dc);
    void kinetic_1der(int ish, int jsh, int c);
    void nuclear(int ish, int jsh);
    void nuclear_1der(int ish, int jsh,
                      Chemistry_QC_GaussianBasis_DerivCenters &dc);
    void nuclear_1der(int ish, int jsh, int c);
    void hcore(int ish, int jsh);
    void hcore_1der(int ish, int jsh,
                    Chemistry_QC_GaussianBasis_DerivCenters &dc);
    void hcore_1der(int ish, int jsh, int c);
                    
};

}

#endif

// Local Variables:
// mode: c++
// End:

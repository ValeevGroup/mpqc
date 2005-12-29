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

#include <chemistry/qc/basis/tbint.h>
#include <chemistry/qc/intcca/int2e.h>
#include <Chemistry_QC_GaussianBasis_IntegralEvaluatorFactory.hh>

using namespace std;
using namespace Chemistry::QC::GaussianBasis;

namespace sc {

/** This implements two body integrals through the CCA interface. */
class TwoBodyIntCCA : public TwoBodyInt {
  protected:
    Ref<Int2eCCA> int2ecca_;

  public:
    TwoBodyIntCCA(Integral*,
		  const Ref<GaussianBasisSet>&b1, 
		  const Ref<GaussianBasisSet>&b2,
		  const Ref<GaussianBasisSet>&b3, 
		  const Ref<GaussianBasisSet>&b4,
                  size_t storage, IntegralEvaluatorFactory, 
                  bool, string );
    ~TwoBodyIntCCA() {};

    int log2_shell_bound(int,int,int,int);
    void compute_shell(int,int,int,int);
    
    size_t storage_used();
    void set_integral_storage(size_t storage);
    int redundant() const { return int2ecca_->redundant(); }
    void set_redundant(int i) { int2ecca_->set_redundant(i); }
};

/** This implements two body derivative integrals 
    through the CCA interface. */
class TwoBodyDerivIntCCA : public TwoBodyDerivInt {
  protected:
    Ref<Int2eCCA> int2ecca_;

  public:
    TwoBodyDerivIntCCA(Integral*,
                  const Ref<GaussianBasisSet>&b1,
                  const Ref<GaussianBasisSet>&b2,
                  const Ref<GaussianBasisSet>&b3,
                  const Ref<GaussianBasisSet>&b4,
                  size_t storage, IntegralEvaluatorFactory, bool, string);
    ~TwoBodyDerivIntCCA() {};

    int log2_shell_bound(int,int,int,int);
    void compute_shell(int,int,int,int,DerivCenters&);

    size_t storage_used();
    void set_integral_storage(size_t storage);
    int redundant() const { return int2ecca_->redundant(); }
    void set_redundant(int i) { int2ecca_->set_redundant(i); }
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:

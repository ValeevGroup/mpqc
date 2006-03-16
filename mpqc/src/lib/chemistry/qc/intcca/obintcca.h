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

#ifndef _chemistry_qc_intcca_obintcca_h
#define _chemistry_qc_intcca_obintcca_h

#include <chemistry/qc/basis/obint.h>
#include <chemistry/qc/intcca/int1e.h>
#include <Chemistry_QC_GaussianBasis_IntegralEvaluatorFactory.hh>
#include <Chemistry_QC_GaussianBasis_DerivCenters.hh>

using namespace Chemistry::QC::GaussianBasis;

namespace sc {

// /////////////////////////////////////////////////////////////////////////

/** This implements one body integrals through the CCA interface. It is
    given a function pointer to the IntCCA member that computes the
    particular integral of interest. */
class OneBodyIntCCA : public OneBodyInt {
  private:
    IntegralEvaluatorFactory eval_factory_;
    bool use_opaque_;
  protected:
    Ref<sc::Int1eCCA> int1ecca_;
    typedef void (sc::Int1eCCA::*IntegralFunction)(int,int);
    IntegralFunction intfunc_;
    Chemistry::QC::GaussianBasis::DerivCenters cca_dc_;
  public:
    OneBodyIntCCA(Integral*,
                 const Ref<GaussianBasisSet>&, const Ref<GaussianBasisSet>&,
                 IntegralEvaluatorFactory, IntegralFunction, bool );
    ~OneBodyIntCCA();
    void compute_shell(int,int);
    bool cloneable();
    Ref<OneBodyInt> clone();
};

///////////////////////////////////////////////////////////////////////////////

 /** This implements one body derivative integrals. It
     is given a function pointer to the Int1eCCA member that computes the
     particular integral of interest. */
class OneBodyDerivIntCCA : public OneBodyDerivInt {
  private:
    IntegralEvaluatorFactory eval_factory_;
    bool use_opaque_;
    ObIntEvalType eval_type_;
  protected:
    Ref<Int1eCCA> int1ecca_;
    typedef void (Int1eCCA::*IntegralFunction)(int, int, DerivCenters&);
    Chemistry::QC::GaussianBasis::DerivCenters cca_dc_;
  public:
    OneBodyDerivIntCCA(Integral*,
                      const Ref<GaussianBasisSet>&,
                      const Ref<GaussianBasisSet>&,
                      IntegralEvaluatorFactory,
                      bool, ObIntEvalType);
    ~OneBodyDerivIntCCA();
    void compute_shell(int, int, DerivCenters&);
    void compute_shell(int, int, int);
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:

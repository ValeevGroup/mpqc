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

#ifndef _chemistry_qc_intcca_obintcca_h
#define _chemistry_qc_intcca_obintcca_h

#include <chemistry/qc/basis/obint.h>
#include <chemistry/qc/intcca/int1e.h>
#include <Chemistry_QC_GaussianBasis_IntegralEvaluatorFactory.hh>

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
    Ref<Int1eCCA> int1ecca_;
    typedef void (Int1eCCA::*IntegralFunction)(int,int);
    IntegralFunction intfunc_;
  public:
    OneBodyIntCCA(Integral*,
                 const Ref<GaussianBasisSet>&, const Ref<GaussianBasisSet>&,
                 IntegralFunction, IntegralEvaluatorFactory, bool );
    ~OneBodyIntCCA();
    void compute_shell(int,int);
    bool cloneable();
    Ref<OneBodyInt> clone();
    void set_buffer( double* );
};

// class PointChargeIntV3 : public OneBodyInt
// {
//   protected:
//     Ref<Int1eV3> int1ev3_;
//     Ref<PointChargeData> data_;
//   public:
//     PointChargeIntV3(Integral*,
//                      const Ref<GaussianBasisSet>&,
//                      const Ref<GaussianBasisSet>&,
//                      const Ref<PointChargeData>&);
//     ~PointChargeIntV3();
//     void compute_shell(int,int);
// };

// class EfieldDotVectorIntV3: public OneBodyInt
// {
//   protected:
//     Ref<Int1eV3> int1ev3_;
//     Ref<EfieldDotVectorData> data_;
//   public:
//     EfieldDotVectorIntV3(Integral*,
//                          const Ref<GaussianBasisSet>&,
//                          const Ref<GaussianBasisSet>&,
//                          const Ref<EfieldDotVectorData>&);
//     ~EfieldDotVectorIntV3();
//     void compute_shell(int,int);
// };

// class DipoleIntV3: public OneBodyInt
// {
//   protected:
//     Ref<Int1eV3> int1ev3_;
//     Ref<DipoleData> data_;
//   public:
//     DipoleIntV3(Integral*,
//                 const Ref<GaussianBasisSet>&,
//                 const Ref<GaussianBasisSet>&,
//                 const Ref<DipoleData>&);
//     ~DipoleIntV3();
//     void compute_shell(int,int);
// };

// // /////////////////////////////////////////////////////////////////////////

// /** This implements one body derivative integrals in the IntV3 library. It
//     is given a function pointer to the Int1eV3 member that computes the
//     particular integral of interest. */
// class OneBodyDerivIntV3 : public OneBodyDerivInt {
//   protected:
//     Ref<Int1eV3> int1ev3_;
//     typedef void (Int1eV3::*IntegralFunction)(int,int,int,int);
//     IntegralFunction intfunc_;
//   public:
//     OneBodyDerivIntV3(Integral*,
//                       const Ref<GaussianBasisSet>&,
//                       const Ref<GaussianBasisSet>&,
//                       IntegralFunction);
//     ~OneBodyDerivIntV3();
//     void compute_shell(int,int,DerivCenters&);
//     void compute_shell(int,int,int);
// };

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:

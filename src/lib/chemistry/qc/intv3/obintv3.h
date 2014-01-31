//
// obintv3.h
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@limitpt.com>
// Maintainer: LPS
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

#ifndef _chemistry_qc_intv3_obintv3_h
#define _chemistry_qc_intv3_obintv3_h

#include <chemistry/qc/basis/obint.h>
#include <chemistry/qc/intv3/int1e.h>

namespace sc {

// /////////////////////////////////////////////////////////////////////////

/** This implements most one body integrals in the IntV3 library. It is
    given a function pointer to the Int1eV3 member that computes the
    particular integral of interest. */
class OneBodyIntV3 : public OneBodyInt {
  protected:
    Ref<Int1eV3> int1ev3_;
    typedef void (Int1eV3::*IntegralFunction)(int,int);
    IntegralFunction intfunc_;
  public:
    OneBodyIntV3(Integral*,
                 const Ref<GaussianBasisSet>&, const Ref<GaussianBasisSet>&,
                 IntegralFunction);
    ~OneBodyIntV3();
    void compute_shell(int,int);
    bool cloneable() const;
    Ref<OneBodyInt> clone();
};

class PointChargeIntV3 : public OneBodyInt
{
  protected:
    Ref<Int1eV3> int1ev3_;
    Ref<PointChargeData> data_;
  public:
    PointChargeIntV3(Integral*,
                     const Ref<GaussianBasisSet>&,
                     const Ref<GaussianBasisSet>&,
                     const Ref<PointChargeData>&);
    ~PointChargeIntV3();
    void compute_shell(int,int);
};

class EfieldIntV3: public OneBodyInt
{
  protected:
    Ref<Int1eV3> int1ev3_;
    Ref<IntParamsOrigin> data_;
  public:
    EfieldIntV3(Integral*,
                const Ref<GaussianBasisSet>&,
                const Ref<GaussianBasisSet>&,
                const Ref<IntParamsOrigin>&);
    ~EfieldIntV3();
    void compute_shell(int,int);
};


class EfieldDotVectorIntV3: public OneBodyInt
{
  protected:
    Ref<Int1eV3> int1ev3_;
    Ref<EfieldDotVectorData> data_;
  public:
    EfieldDotVectorIntV3(Integral*,
                         const Ref<GaussianBasisSet>&,
                         const Ref<GaussianBasisSet>&,
                         const Ref<EfieldDotVectorData>&);
    ~EfieldDotVectorIntV3();
    void compute_shell(int,int);
};

class DipoleIntV3: public OneBodyInt
{
  protected:
    Ref<Int1eV3> int1ev3_;
    Ref<DipoleData> data_;
  public:
    DipoleIntV3(Integral*,
                const Ref<GaussianBasisSet>&,
                const Ref<GaussianBasisSet>&,
                const Ref<DipoleData>&);
    ~DipoleIntV3();
    void compute_shell(int,int);
};

// /////////////////////////////////////////////////////////////////////////

/** This implements one body derivative integrals in the IntV3 library. It
    is given a function pointer to the Int1eV3 member that computes the
    particular integral of interest. */
class OneBodyDerivIntV3 : public OneBodyDerivInt {
  protected:
    Ref<Int1eV3> int1ev3_;
    typedef void (Int1eV3::*IntegralFunction)(int,int,int,int);
    IntegralFunction intfunc_;
  public:
    OneBodyDerivIntV3(Integral*,
                      const Ref<GaussianBasisSet>&,
                      const Ref<GaussianBasisSet>&,
                      IntegralFunction);
    ~OneBodyDerivIntV3();
    void compute_shell(int,int,int);
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:

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

///////////////////////////////////////////////////////////////////////////

class OneBodyIntV3 : public OneBodyInt {
  protected:
    RefInt1eV3 int1ev3_;
    typedef void (Int1eV3::*IntegralFunction)(int,int);
    IntegralFunction intfunc_;
  public:
    OneBodyIntV3(Integral*,
                 const RefGaussianBasisSet&, const RefGaussianBasisSet&,
                 IntegralFunction);
    ~OneBodyIntV3();
    void compute_shell(int,int);
};

class PointChargeIntV3 : public OneBodyInt
{
  protected:
    RefInt1eV3 int1ev3_;
    RefPointChargeData data_;
  public:
    PointChargeIntV3(Integral*,
                     const RefGaussianBasisSet&,
                     const RefGaussianBasisSet&,
                     const RefPointChargeData&);
    ~PointChargeIntV3();
    void compute_shell(int,int);
};

class EfieldDotVectorIntV3: public OneBodyInt
{
  protected:
    RefInt1eV3 int1ev3_;
    RefEfieldDotVectorData data_;
  public:
    EfieldDotVectorIntV3(Integral*,
                         const RefGaussianBasisSet&,
                         const RefGaussianBasisSet&,
                         const RefEfieldDotVectorData&);
    ~EfieldDotVectorIntV3();
    void compute_shell(int,int);
};

class DipoleIntV3: public OneBodyInt
{
  protected:
    RefInt1eV3 int1ev3_;
    RefDipoleData data_;
  public:
    DipoleIntV3(Integral*,
                const RefGaussianBasisSet&,
                const RefGaussianBasisSet&,
                const RefDipoleData&);
    ~DipoleIntV3();
    void compute_shell(int,int);
};

///////////////////////////////////////////////////////////////////////////

class OneBodyDerivIntV3 : public OneBodyDerivInt {
  protected:
    RefInt1eV3 int1ev3_;
    typedef void (Int1eV3::*IntegralFunction)(int,int,int,int);
    IntegralFunction intfunc_;
  public:
    OneBodyDerivIntV3(Integral*,
                      const RefGaussianBasisSet&,
                      const RefGaussianBasisSet&,
                      IntegralFunction);
    ~OneBodyDerivIntV3();
    void compute_shell(int,int,DerivCenters&);
    void compute_shell(int,int,int);
};

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:

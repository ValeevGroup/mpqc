//
// tbintv3.h
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

#ifndef _chemistry_qc_intv3_tbintv3_h
#define _chemistry_qc_intv3_tbintv3_h

#include <chemistry/qc/basis/tbint.h>
#include <chemistry/qc/intv3/int2e.h>

namespace sc {

/** This implements electron repulsion integrals in the IntV3 library. */
class TwoBodyIntV3 : public TwoBodyInt {
  protected:
    Ref<Int2eV3> int2ev3_;
    Ref<TwoBodyOperSetDescr> descr_;

  public:
    TwoBodyIntV3(Integral*integral,
                 const Ref<GaussianBasisSet>&b1,
                 const Ref<GaussianBasisSet>&b2,
                 const Ref<GaussianBasisSet>&b3,
                 const Ref<GaussianBasisSet>&b4,
                 size_t storage);
    ~TwoBodyIntV3();

    TwoBodyOperSet::type type() const { return TwoBodyOperSet::ERI; }
    const Ref<TwoBodyOperSetDescr>& descr() const { return descr_; }

    int log2_shell_bound(int,int,int,int);
    void compute_shell(int,int,int,int);

    size_t storage_used() { return int2ev3_->storage_used(); }
    void set_integral_storage(size_t storage);
};

/** This implements electron repulsion integrals involving three centers in
 * the IntV3 library. */
class TwoBodyThreeCenterIntV3 : public TwoBodyThreeCenterInt {
  protected:
    Ref<Int2eV3> int2ev3_;
    Ref<TwoBodyOperSetDescr> descr_;

  public:
    TwoBodyThreeCenterIntV3(Integral*integral,
                            const Ref<GaussianBasisSet>&b1,
                            const Ref<GaussianBasisSet>&b2,
                            const Ref<GaussianBasisSet>&b3,
                            size_t storage);
    ~TwoBodyThreeCenterIntV3();

    TwoBodyOperSet::type type() const { return TwoBodyOperSet::ERI; }
    const Ref<TwoBodyOperSetDescr>& descr() const { return descr_; }

    int log2_shell_bound(int,int,int);
    void compute_shell(int,int,int);

    size_t storage_used() { return int2ev3_->storage_used(); }
    void set_integral_storage(size_t storage);
};

/** This implements electron repulsion integrals involving two centers in
 * the IntV3 library. */
class TwoBodyTwoCenterIntV3 : public TwoBodyTwoCenterInt {
  protected:
    Ref<Int2eV3> int2ev3_;
    Ref<TwoBodyOperSetDescr> descr_;

  public:
    TwoBodyTwoCenterIntV3(Integral*integral,
                          const Ref<GaussianBasisSet>&b1,
                          const Ref<GaussianBasisSet>&b2,
                          size_t storage);
    ~TwoBodyTwoCenterIntV3();

    TwoBodyOperSet::type type() const { return TwoBodyOperSet::ERI; }
    const Ref<TwoBodyOperSetDescr>& descr() const { return descr_; }

    int log2_shell_bound(int,int);
    void compute_shell(int,int);

    size_t storage_used() { return int2ev3_->storage_used(); }
    void set_integral_storage(size_t storage);
};

/** This implements electron repulsion derivative integrals in the IntV3
    library. */
class TwoBodyDerivIntV3 : public TwoBodyDerivInt {
  protected:
    Ref<Int2eV3> int2ev3_;

  public:
    TwoBodyDerivIntV3(Integral*integral,
                      const Ref<GaussianBasisSet>&b1,
                      const Ref<GaussianBasisSet>&b2,
                      const Ref<GaussianBasisSet>&b3,
                      const Ref<GaussianBasisSet>&b4,
                      size_t storage);
    ~TwoBodyDerivIntV3();

    int log2_shell_bound(int,int,int,int);
    void compute_shell(int,int,int,int,DerivCenters&);

    size_t storage_used() { return int2ev3_->storage_used(); }
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:

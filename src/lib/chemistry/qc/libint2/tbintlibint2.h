//
// tbintlibint2.h
//
// Copyright (C) 2001 Edward Valeev
//
// Author: Edward Valeev <evaleev@vt.edu>
// Maintainer: EV
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

#ifndef _chemistry_qc_libint2_tbint_h
#define _chemistry_qc_libint2_tbint_h

#include <chemistry/qc/basis/tbint.h>
#include <chemistry/qc/libint2/int2e.h>
#include <chemistry/qc/libint2/tbosar.h>

namespace sc {

/** This implements 4-center two-electron integrals in the IntLibint2 library. */
class TwoBodyIntLibint2 : public TwoBodyInt {

    TwoBodyOperSet::type int2etype_;
    Ref<TwoBodyOperSetDescr> descr_;
    Ref<IntParams> params_;

  protected:
    Ref<Int2eLibint2> int2elibint2_;

  public:
    typedef IntParamsG12::ContractedGeminal ContractedGeminal;
    TwoBodyIntLibint2(Integral*integral,
                 const Ref<GaussianBasisSet>&b1,
                 const Ref<GaussianBasisSet>&b2,
                 const Ref<GaussianBasisSet>&b3,
                 const Ref<GaussianBasisSet>&b4,
                 size_t storage, TwoBodyOperSet::type int2etype,
                 const Ref<IntParams>& params);
    virtual ~TwoBodyIntLibint2();

    TwoBodyOperSet::type type() const { return int2etype_; }
    const Ref<TwoBodyOperSetDescr>& descr() const { return descr_; }

    int log2_shell_bound(int,int,int,int);
    void compute_shell(int,int,int,int);

    size_t used_storage() const { return int2elibint2_->storage_used(); }
    void set_integral_storage(size_t storage);

    const double *buffer(TwoBodyOper::type te_type) const {
      return int2elibint2_->buffer( descr_->opertype(te_type) );
    }

    bool cloneable() const;
    Ref<TwoBodyInt> clone();

  private:
    /// similar to the standard constructor, but saves some work (and pain)
    /// for clone() method
    TwoBodyIntLibint2(Integral*integral,
                      const Ref<GaussianBasisSet>&b1,
                      const Ref<GaussianBasisSet>&b2,
                      const Ref<GaussianBasisSet>&b3,
                      const Ref<GaussianBasisSet>&b4,
                      size_t storage, TwoBodyOperSet::type int2etype,
                      const Ref<IntParams>& params,
                      const Ref<Log2Bounds>& bounds);
};

/** This implements 3-center 2-body integrals in the IntLibint2 library. */
class TwoBodyThreeCenterIntLibint2 : public TwoBodyThreeCenterInt {

  TwoBodyOperSet::type int2etype_;
  Ref<TwoBodyOperSetDescr> descr_;
  Ref<IntParams> params_;

  protected:
    Ref<Int2eLibint2> int2elibint2_;

  public:
    typedef IntParamsG12::ContractedGeminal ContractedGeminal;
    TwoBodyThreeCenterIntLibint2(Integral*integral,
                 const Ref<GaussianBasisSet>&b1,
                 const Ref<GaussianBasisSet>&b2,
                 const Ref<GaussianBasisSet>&b3,
                 size_t storage, TwoBodyOperSet::type int2etype,
                 const Ref<IntParams>& params);
    virtual ~TwoBodyThreeCenterIntLibint2();

    TwoBodyOperSet::type type() const { return int2etype_; }
    const Ref<TwoBodyOperSetDescr>& descr() const { return descr_; }

    int log2_shell_bound(int,int,int);
    void compute_shell(int,int,int);

    size_t used_storage() const { return int2elibint2_->storage_used(); }
    void set_integral_storage(size_t storage);

    const double *buffer(TwoBodyOper::type te_type) const {
      return int2elibint2_->buffer( descr_->opertype(te_type) );
    }

    bool cloneable() const;
    Ref<TwoBodyThreeCenterInt> clone();

  private:
    /// similar to the standard constructor, but saves some work (and pain)
    /// for clone() method
    TwoBodyThreeCenterIntLibint2(Integral*integral,
                 const Ref<GaussianBasisSet>&b1,
                 const Ref<GaussianBasisSet>&b2,
                 const Ref<GaussianBasisSet>&b3,
                 const Ref<GaussianBasisSet>&bunit,
                 size_t storage, TwoBodyOperSet::type int2etype,
                 const Ref<IntParams>& params,
                 const Ref<Log2Bounds>& bounds);

};

/** This implements 2-center 2-body integrals in the IntLibint2 library. */
class TwoBodyTwoCenterIntLibint2 : public TwoBodyTwoCenterInt {

  TwoBodyOperSet::type int2etype_;
  Ref<TwoBodyOperSetDescr> descr_;
  Ref<IntParams> params_;

  protected:
    Ref<Int2eLibint2> int2elibint2_;

  public:
    typedef IntParamsG12::ContractedGeminal ContractedGeminal;
    TwoBodyTwoCenterIntLibint2(Integral*integral,
                 const Ref<GaussianBasisSet>&b1,
                 const Ref<GaussianBasisSet>&b2,
                 size_t storage, TwoBodyOperSet::type int2etype,
                 const Ref<IntParams>& params);
    virtual ~TwoBodyTwoCenterIntLibint2();

    TwoBodyOperSet::type type() const { return int2etype_; }
    const Ref<TwoBodyOperSetDescr>& descr() const { return descr_; }

    bool cloneable() const;
    Ref<TwoBodyTwoCenterInt> clone();

    int log2_shell_bound(int,int);
    void compute_shell(int,int);

    size_t used_storage() const { return int2elibint2_->storage_used(); }
    void set_integral_storage(size_t storage);

    const double *buffer(TwoBodyOper::type te_type) const {
      return int2elibint2_->buffer( descr_->opertype(te_type) );
    }

  private:
    /// similar to the standard constructor, but saves some work (and pain)
    /// for clone() method
    TwoBodyTwoCenterIntLibint2(Integral*integral,
                 const Ref<GaussianBasisSet>&b1,
                 const Ref<GaussianBasisSet>&b2,
                 const Ref<GaussianBasisSet>&bunit,
                 size_t storage, TwoBodyOperSet::type int2etype,
                 const Ref<IntParams>& params,
                 const Ref<Log2Bounds>& bounds);
};

/** This implements electron repulsion derivative integrals in the IntV3
    library. */
class TwoBodyDerivIntLibint2 : public TwoBodyDerivInt {
  protected:
    Ref<Int2eLibint2> int2elibint2_;

  public:
    TwoBodyDerivIntLibint2(Integral*integral,
                      const Ref<GaussianBasisSet>&b1,
                      const Ref<GaussianBasisSet>&b2,
                      const Ref<GaussianBasisSet>&b3,
                      const Ref<GaussianBasisSet>&b4,
                      size_t storage, TwoBodyOperSet::type int2etype);
    virtual ~TwoBodyDerivIntLibint2();

    int log2_shell_bound(int,int,int,int);
    void compute_shell(int,int,int,int,DerivCenters&);

    size_t used_storage() const { return int2elibint2_->storage_used(); }
};

    namespace libint2 {

      template <class Int2e>
      struct Int2eCreator {
          Ref<Int2e> operator()(Integral*integral,
		     const Ref<GaussianBasisSet>& b1,
		     const Ref<GaussianBasisSet>& b2,
		     const Ref<GaussianBasisSet>& b3,
		     const Ref<GaussianBasisSet>& b4,
		     size_t storage,
		     const Ref<IntParams>& params) {
            return new Int2e(integral,b1,b2,b3,b4,storage,params);
          }
      };

    }

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:

//
// intdescr.h
//
// Copyright (C) 2005 Edward Valeev
//
// Author: Edward Valeev <edward.valeev@chemistry.gatech.edu>
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

#ifdef __GNUG__
#pragma interface
#endif

#include <chemistry/qc/basis/integral.h>
#include <chemistry/qc/basis/intparams.h>

#ifndef _chemistry_qc_basis_intdescr_h
#define _chemistry_qc_basis_intdescr_h

namespace sc {
  
  /** For a set of integrals (e.g. 1/r12 and exp(-gamma*r12)), IntegralSetDescr describes properties of
      integral types, their parameters (e.g. gamma), and how to construct an evaluator object.
      IntEval is the base evaluator type (e.g. TwoBodyInt) for all integrals in the integral set. */
  template <typename IntEval>
  class IntegralSetDescr : public RefCount {
    public:
    IntegralSetDescr() {}

    /// call appropriate method to produce TwoBodyInt corresponding to this Set
    virtual Ref<IntEval> inteval() const =0;
    /// how many integral sets
    virtual unsigned int num_sets() const =0;
    /// iterator over Property classes for each Integral in the set (Property describes symmetry, labels, etc.)
    //PropertyIterator propiter();
  };
  
  /// TwoBodyIntDescr describes a TwoBodyInt
  typedef IntegralSetDescr<TwoBodyInt> TwoBodyIntDescr;
  
  /** TwoBodyIntDescrERI describes single set of electron repulsion integrals
    */
  class TwoBodyIntDescrERI : public TwoBodyIntDescr {
    public:
    TwoBodyIntDescrERI(const Ref<Integral>& IF);
    
    /// which factory is used
    Ref<Integral> factory() const { return factory_; }
    /// implementation of Integral::inteval()
    Ref<TwoBodyInt> inteval() const;
    /// implementation of Integral::num_sets()
    unsigned int num_sets() const { return 1; }

    private:
    /// which factory is used
    Ref<Integral> factory_;
  };
  
  /** TwoBodyIntDescrR12 describes a set of integrals used in MP2-F12 theories
      using linear r12 correlation factor. The following integrals are computed:
      1) 1/r_{12} 2) r_{12} 3) [T_1,r_{12}] 4) [T_2,r_{12}]
    */
  class TwoBodyIntDescrR12 : public TwoBodyIntDescr {
    public:
    TwoBodyIntDescrR12(const Ref<Integral>& IF);
    
    /// which factory is used
    Ref<Integral> factory() const { return factory_; }
    /// implementation of Integral::inteval()
    Ref<TwoBodyInt> inteval() const;
    /// implementation of Integral::num_sets()
    unsigned int num_sets() const { return 4; }

    private:
    /// which factory is used
    Ref<Integral> factory_;
  };
  
  /** TwoBodyIntDescrG12 describes a set of integrals used in MP2-F12 theories
      using Gaussian geminal correlation factors. The following integrals are computed:
      1) 1/r_{12} 2) g_{12} = exp(-gamma * r_{12}^2) 3) [T_1,g_{12}]
      4) [T_2,g_{12}] 5) g_{12}/r_{12} 6) [g_{12}',[T_1,g_{12}]] */
  class TwoBodyIntDescrG12 : public TwoBodyIntDescr {
    public:
    /// Compute integrals using geminal parameters params
    TwoBodyIntDescrG12(const Ref<Integral>& IF,
                        const Ref<IntParamsG12>& g12params);
    
    /// which factory is used
    Ref<Integral> factory() const { return factory_; }
    /// implementation of Integral::inteval()
    Ref<TwoBodyInt> inteval() const;
    /// implementation of Integral::num_sets()
    unsigned int num_sets() const { return 6; }

    private:
    /// geminal parameters
    Ref<IntParamsG12> params_;
    /// which factory is used
    Ref<Integral> factory_;
  };

  /** TwoBodyIntDescrG12NC describes a set of integrals used in MP2-F12 theories
      using Gaussian geminal correlation factors (without kinetic energy commutators).
      The following integrals are computed:
      1) 1/r_{12} 2) g_{12} = exp(-gamma * r_{12}^2)
      3) g_{12}/r_{12} 4) [g_{12}',[T_1,g_{12}]]
      5) (exp(g12')-exp(g12))/(exp(g12)+exp(g12')) g12 * g12'
    */
  class TwoBodyIntDescrG12NC : public TwoBodyIntDescr {
    public:
    /// Compute integrals using geminal parameters params
    TwoBodyIntDescrG12NC(const Ref<Integral>& IF,
                         const Ref<IntParamsG12>& g12params);

    /// which factory is used
    Ref<Integral> factory() const { return factory_; }
    /// implementation of Integral::inteval()
    Ref<TwoBodyInt> inteval() const;
    /// implementation of Integral::num_sets()
    unsigned int num_sets() const { return 5; }

    private:
    /// geminal parameters
    Ref<IntParamsG12> params_;
    /// which factory is used
    Ref<Integral> factory_;
  };

  /** TwoBodyIntDescrGenG12 describes a set of integrals used in MP2-F12 theories
      using general Gaussian geminal correlation factors. The following integrals are computed:
      1) 1/r_{12} 2) g_{12} = exp(-gamma * r_{12}^2)
      3) g_{12}/r_{12} 4) r1.r1 g12*g12' 5) r1.r2 g12*g12'
      6) [g_{12}',[T_1,g_{12}]] */
  class TwoBodyIntDescrGenG12 : public TwoBodyIntDescr {
    public:
    /// Compute integrals using geminal parameters params
    TwoBodyIntDescrGenG12(const Ref<Integral>& IF,
                          const Ref<IntParamsGenG12>& gg12params);

    /// which factory is used
    Ref<Integral> factory() const { return factory_; }
    /// implementation of Integral::inteval()
    Ref<TwoBodyInt> inteval() const;
    /// implementation of Integral::num_sets()
    unsigned int num_sets() const { return 6; }

    private:
    /// geminal parameters
    Ref<IntParamsGenG12> params_;
    /// which factory is used
    Ref<Integral> factory_;
  };

}

#endif


//
// intparams.h
//
// Copyright (C) 2005 Edward Valeev
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

#ifndef _chemistry_qc_basis_intparams_h
#define _chemistry_qc_basis_intparams_h

#include <vector>
#include <utility>
#include <util/state/state.h>
#include <util/state/statein.h>
#include <util/state/stateout.h>

namespace sc {

  /** This class passes optional operator parameters. These parameters can be
      passed to
      1) Factory methods which initialize XXXBodyInt objects
      2) compute_shell methods of XXXBodyInt objects
  */
  class IntParams : virtual public SavableState {
    public:
      IntParams(unsigned int nparams=0);
      IntParams(StateIn&);
      virtual ~IntParams();
      void save_data_state(StateOut&);

      unsigned int nparams() const;

    protected:
      template <typename T> const T* downcast(const IntParams& p) const {
        try {
          const T& voidref = dynamic_cast<const T&>(p);
          return &voidref;
        }
        catch (std::bad_cast&) { return 0; }
      }

    private:
      unsigned int nparams_;

      friend bool operator==(const IntParams& p1, const IntParams& p2);
      virtual bool equiv(const IntParams& other) const =0;
  };
  inline bool operator==(const IntParams& p1, const IntParams& p2) {
    return p1.equiv(p2);
  }

  /** Passes params to Integral::electron_repulsion() and other factory methods which do not need parameters */
  class IntParamsVoid : public IntParams {
    public:
      IntParamsVoid();
      IntParamsVoid(StateIn&);
      ~IntParamsVoid();
      void save_data_state(StateOut&);
    private:
      bool equiv(const IntParams& other) const;
  };

  /** Passes params to Integral::g12() */
  class IntParamsG12 : public IntParams {
    public:
      ///     std::pair<  g,     c   >  as in c * exp( - g*r12)
      typedef std::pair<double,double> PrimitiveGeminal;
      typedef std::vector<PrimitiveGeminal> ContractedGeminal;
      /// 1 = e^(-0.0 * r_{12})
      static ContractedGeminal zero_exponent_geminal;
      /// null (i.e., invalid) geminal
      static ContractedGeminal null_geminal;

      /// Request integrals with only 1 geminal (g12, g12/r12, [Ti,g12])
      IntParamsG12(const ContractedGeminal& bra);
      /// Request integrals with 2 geminals (g12*g12', [g12,[T1,g12']], [Ti,g12*g12'])
      IntParamsG12(const ContractedGeminal& bra,
                   const ContractedGeminal& ket);
      IntParamsG12(StateIn&);
      ~IntParamsG12();

      void save_data_state(StateOut&);

      const ContractedGeminal& bra() const;
      const ContractedGeminal& ket() const;

      static PrimitiveGeminal product(const PrimitiveGeminal& A,
                                      const PrimitiveGeminal& B);
      static ContractedGeminal product(const ContractedGeminal& A,
                                       const ContractedGeminal& B);

    private:
      bool equiv(const IntParams& other) const;

      /// An invalid exponent (exponent of null_geminal)
      static double null_exponent;

      ContractedGeminal bra_;
      ContractedGeminal ket_;
  };

  /** Passes params to Integral::geng12() */
  class IntParamsGenG12 : public IntParams {
    public:
      ///     std::pair<std::pair<  a,     g   >   c,   >  as in c * exp( - a*(r1^2+r2^2) - g*r12^2)
      typedef std::pair<std::pair<double,double>, double> PrimitiveGeminal;
      typedef std::vector<PrimitiveGeminal> ContractedGeminal;
      /// 1 = e^(- 0*(r1+r2) - 0*r12)
      static ContractedGeminal zero_exponent_geminal;
      /// null (i.e., invalid) geminal
      static ContractedGeminal null_geminal;

      /// Request integrals with only 1 geminal (g12, g12/r12, [Ti,g12])
      IntParamsGenG12(const ContractedGeminal& bra);
      /// Request integrals with 2 geminals (g12*g12', [g12,[T1,g12']], [Ti,g12*g12'])
      /// IntParamsGenG12 can only be used to compute integrals with the same integral
      /// in bra and ket, hence this will throw if bra != ket.
      IntParamsGenG12(const ContractedGeminal& bra,
                   const ContractedGeminal& ket);
      IntParamsGenG12(StateIn&);
      ~IntParamsGenG12();
      void save_data_state(StateOut&);

      const ContractedGeminal& bra() const;
      const ContractedGeminal& ket() const;

    private:
      bool equiv(const IntParams& other) const;

      /// An invalid exponent (exponent of null_geminal)
      static double null_exponent;

      ContractedGeminal bra_;
      ContractedGeminal ket_;
  };


}

#endif

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:

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

  class EfieldDotVectorData: public RefCount
  {
    public:
      EfieldDotVectorData() {};
      ~EfieldDotVectorData();

      double position[3];
      double vector[3];

      void set_position(double*);
      void set_vector(double*);
  };


  class PointChargeData: public RefCount
  {
    private:
      int ncharges_;
      const double *charges_;
      const double *const*positions_;
      double *alloced_charges_;
      double **alloced_positions_;

    public:
      // If copy_data is 0, the passed positions and charges will
      // be stored (but not freed).
      PointChargeData(int ncharge,
                      const double *const*positions, const double *charges,
                      int copy_data = 0);
      ~PointChargeData();

      int ncharges() const { return ncharges_; }
      const double *charges() const { return charges_; }
      const double *const*positions() const { return positions_; }
  };

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

      /// Some IntParams implementations can have a variable number of params; otherwise this will return 0
      unsigned int nparams() const;

    protected:
      template <typename T> const T* downcast(const IntParams& p) const {
        const T* castptr = dynamic_cast<const T*>(&p);
        return castptr;
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

  /** Passes params to Integral::dipole() and other factory methods which need r information */
  class IntParamsOrigin : public IntParams {
    public:
      IntParamsOrigin();
      IntParamsOrigin(const double (&O)[3]);
      template <typename Real> explicit IntParamsOrigin(const Real* O) : O_(3) {
        std::copy(O, O+3, O_.begin());
      }
      IntParamsOrigin(StateIn&);
      ~IntParamsOrigin();
      void save_data_state(StateOut&);

      const double* r() const;
      double r(unsigned int xyz) const;

    private:
      std::vector<double> O_;
      bool equiv(const IntParams& other) const;
  };

  /** Used to pass params to Integral::g12().
      Represents (a pair of) linear combinations of exponentials:
      \f$
         \sum_i c_i \exp( - g_i x )
      \f$
      where \f$ x \f$ is \f$ r_{12}^2 \f$ for Gaussian geminals.
    */
  class IntParamsG12 : public IntParams {
    public:
      ///     std::pair<  g,     c   >  as in c * exp( - g*r12^2)
      typedef std::pair<double,double> PrimitiveGeminal;
      typedef std::vector<PrimitiveGeminal> ContractedGeminal;
      /// 1 = e^(-0.0 * r_{12}^2)
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

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:

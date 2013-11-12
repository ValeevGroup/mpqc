//
// gaussianfit.h
//
// Copyright (C) 2007 Edward Valeev
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

#ifndef _math_optimize_gaussianfit_h
#define _math_optimize_gaussianfit_h

#include <stdexcept>
#include <vector>
#include <cmath>

namespace sc {

  /** GaussianFit<Function> is a fit of Function(x)*Weight(x) to N Gaussians on range [left,right]
   Valid Function and Weight are Unary Functions which take and return a double.
   */
  template<typename Function, typename Weight>
  class GaussianFit {
    public:
      // If 1, "weigh" the function, else weigh the square of the error
      static const bool weigh_F = 0;

      typedef double Exp;
      typedef double Coef;
      typedef std::pair<Exp, Coef> Gaussian;
      typedef std::vector<Gaussian> Gaussians;
      /// Will fit a function F with weight W evenly sampled on NP points in an interval [left,right] to N Gaussians
      GaussianFit(unsigned int N, const Weight& W, double left, double right,
                  unsigned int NP);
      ~GaussianFit();

      /// fit F and return the parameters of the fit.
      const Gaussians& operator()(const Function& F) const;

    private:
      /// expand Function in terms of Gaussians of form x^k_ * exp(-a*x^2)
      static const int k_ = 0;

      /// weight function
      Weight weight_;
      /// (exponent,coefficient) pairs. One of 2 representations, hence must be mutable
      mutable std::vector<std::pair<Exp, Coef> > gaussians_;
      double left_;
      double right_;
      unsigned int npts_;

      /// parameter values in format suitable for lmfit. One of 2 representations, hence must be mutable
      mutable double* p_;
      /// scratch space, holds zeros
      double* scratch_;

      /// Change to nonzero to debug
      static const int classdebug_ = 0;

      /// only changes representation, hence logically const
      void extract_params() const;
      /// only changes representation, hence logically const
      void assign_params() const;
  };

  namespace math {
    /// Slater1D(k,x) = \f$ c x^k \exp(-a*x) \f$
    class Slater1D {
      public:
        Slater1D(double a, int k = 0, double c = 1.0) :
          a_(a), k_(k), c_(c) {
        }
        double operator()(double x) const {
          return c_ * std::pow(x, k_) * std::exp(-a_ * x);
        }
      private:
        int k_;
        double a_;
        double c_;
    };
    /// Gaussian1D(k,x) = c x^k exp(-a*x^2)
    class Gaussian1D {
      public:
        Gaussian1D(double a, int k = 0, double c = 1.0) :
          a_(a), k_(k), c_(c) {
        }
        double operator()(double x) const {
          return c_ * std::pow(x, k_) * std::exp(-a_ * x * x);
        }
      private:
        int k_;
        double a_;
        double c_;
    };
    /// PowerExponential1D(k,l,x) = c x^k exp(-a*x^l)
    class PowerExponential1D {
      public:
        PowerExponential1D(double a, int l = 2, int k = 0, double c = 1.0) :
          a_(a), k_(k), l_(l), c_(c) {
        }
        double operator()(double x) const {
          return c_ * std::pow(x, k_) * std::exp(-a_ * std::pow(x, l_));
        }
      private:
        int k_;
        int l_;
        double a_;
        double c_;
    };
  }

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:

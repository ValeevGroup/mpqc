//
// gaussianfit.timpl.h
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

#ifndef _math_optimize_gaussianfittimpl_h
#define _math_optimize_gaussianfittimpl_h

#include <cmath>
#include <algorithm>
#include <math/optimize/levmar/lm.h>
#include <math/optimize/gaussianfit.h>
#include <util/misc/scexception.h>

extern "C" {
  void __eval_slater(double* params, double* f, int nparam, int np,
                     void *extraparams);
  void __eval_slater_dfdp(double* params, double* dfdp, int nparam, int np,
                          void *extraparams);
  void __eval_slater_pgauss(double* params, double* f, int nparam, int np,
                            void *extraparams);
  void __eval_slater_dfdp_pgauss(double* params, double* dfdp, int nparam,
                                 int np, void *extraparams);
}
;

namespace {
  extern "C" {
    typedef void
        (*eval_f_ptr)(double *p, double *hx, int m, int n, void *adata);
    typedef void (*eval_dfdp_ptr)(double *p, double *hx, int m, int n,
                                  void *adata);
  }
  ;
}

namespace sc {
  
  template <typename Function, typename Weight> GaussianFit<Function, Weight>::GaussianFit(
                                                                                           unsigned int N,
                                                                                           const Weight& W,
                                                                                           double left,
                                                                                           double right,
                                                                                           unsigned int NP) :
    weight_(W), gaussians_(), left_(left), right_(right), npts_(NP), p_(new double[2*N]), scratch_(new double[NP])
        {
    if (left < 0.0 || right < 0.0 || left >= right || N < 1 || NP < 2)
      throw sc::ProgrammingError("GaussianFit::GaussianFit() -- invalid parameters",__FILE__,__LINE__);
    
    // Set initial parameters
    gaussians_.resize(N);
    // center exponents at 3 and use exp. ratio of 3
    if (N%2) {
      // Odd N
      double gamma = 3.0;
      for (int i=N/2; i>=0; --i, gamma/=3.0)
        gaussians_[i] = std::make_pair(gamma, 1.0);
      gamma = 9.0;
      for (int i=N/2+1; i<N; ++i, gamma*=3.0)
        gaussians_[i] = std::make_pair(gamma, 1.0);
    } else {
      // Even N
      double gamma = 3.0/sqrt(3.0);
      for (int i=N/2-1; i>=0; --i, gamma/=3.0)
        gaussians_[i] = std::make_pair(gamma, 1.0);
      gamma = 3.0*sqrt(3.0);
      for (int i=N/2; i<N; ++i, gamma*=3.0)
        gaussians_[i] = std::make_pair(gamma, 1.0);
    }
    
    // zero out scratch
    for (unsigned int p=0; p<npts_; ++p)
      scratch_[p] = 0.0;
  }
  
  template <typename Function, typename Weight> GaussianFit<Function, Weight>::~GaussianFit() {
    delete[] p_;
    delete[] scratch_;
  }
  
  namespace {
    
    template <class Function, class Weight> struct ExtraParams {
        const Function* f; // functor
        const Weight* w; // weight
        unsigned int npts; // number of points
        double lb; // lower bound
        double ub; // upper bound
        int k; // fit to x^k exp(-a*x^2) functions
    };
    
    template <class Function, class Weight> void eval_f(double* params,
                                                        double* f, int nparam,
                                                        int np,
                                                        void *extraparams) {
      typedef ExtraParams<Function,Weight> XPS;
      typedef GaussianFit<Function,Weight> GaussFit;
      XPS* XParams = static_cast<XPS*>(extraparams);
      const double lb = XParams->lb;
      const double ub = XParams->ub;
      const double dx = (ub - lb)/(np - 1);
      const double npts = XParams->npts;
      const Function& F = *(XParams->f);
      const Weight& W = *(XParams->w);
      const int k = XParams->k;
      
      double x = 0.0;
      for (unsigned int p=0; p<npts; ++p, x+=dx) {
        double fitvalue = 0.0;
        for (unsigned int i=0; i<nparam;) {
          const double gamma = params[i++];
          const double coef = params[i++];
          fitvalue += coef * std::pow(x, k) * std::exp(-gamma*x*x);
        }
        
        if (GaussFit::weigh_F) {
          const double df = fitvalue - W(x) * F(x);
          f[p] = df;
        } else {
          const double df = fitvalue - F(x);
          // lm will minimize the difference between the target function (identically 0.0) and the fitting function (Function - fit)
          // to weigh the sum of squares by W need to multiply this by the square root of W
          f[p] = df * std::sqrt(W(x));
        }
      }
    }
    template <class Function, class Weight> void eval_dfdp(double* params,
                                                           double* dfdp,
                                                           int nparam, int np,
                                                           void *extraparams) {
      typedef ExtraParams<Function,Weight> XPS;
      typedef GaussianFit<Function,Weight> GaussFit;
      XPS* XParams = static_cast<XPS*>(extraparams);
      const double lb = XParams->lb;
      const double ub = XParams->ub;
      const double dx = (ub - lb)/(np - 1);
      const double npts = XParams->npts;
      const Function& F = *(XParams->f);
      const Weight& W = *(XParams->w);
      const int k = XParams->k;
      
      double x = 0.0;
      unsigned int j = 0;
      for (unsigned int p=0, j=0; p<npts; ++p, x+=dx) {
        const double sqrt_W = std::sqrt(W(x));
        for (unsigned int i=0; i<nparam;) {
          const unsigned int p1 = i++;
          const unsigned int p2 = i++;
          const double gamma = params[p1];
          const double coef = params[p2];
          const double expgxx = std::pow(x, k) * std::exp(-gamma*x*x);
          if (GaussFit::weigh_F) {
            dfdp[j++] = - (x * x * coef * expgxx);
            dfdp[j++] = expgxx;
          } else {
            dfdp[j++] = -sqrt_W * (x * x * coef * expgxx);
            dfdp[j++] = sqrt_W * expgxx;
          }
        }
      }
    }
  } // namespace

  namespace detail {
    template <class Function, class Weight> struct __to_extern_C_eval {
        static eval_f_ptr f_ptr;
        static eval_dfdp_ptr dfdp_ptr;
    };
  
  }
  
  template <typename Function, typename Weight> const typename GaussianFit<Function,Weight>::Gaussians& GaussianFit<
      Function, Weight>::operator()(const Function& F) const
  {
    ExtraParams<Function,Weight> xp;
    xp.f = &F;
    xp.w = &weight_;
    xp.npts = npts_;
    xp.lb = left_;
    xp.ub = right_;
    xp.k = k_;
    const int ngaussians = gaussians_.size();

    extract_params();

    const int niter = dlevmar_der(detail::__to_extern_C_eval<Function,Weight>::f_ptr,
                                  detail::__to_extern_C_eval<Function,Weight>::dfdp_ptr,
        p_, scratch_, 2*ngaussians, npts_, 100000, NULL, NULL, NULL, NULL, static_cast<void*>(&xp));
    if (niter> 0) {
      if (classdebug_) {
        std::cout << "GaussianFit() -- converged in " << niter << " iterations" << std::endl;
        for(unsigned int g=0, p=0; g<ngaussians; ++g) {
          std::cout << " gamma = " << p_[p++];
          std::cout << " coef = " << p_[p++] << std::endl;
        }
      }

      assign_params();
    }
    else {
      throw AlgorithmException("GaussianFit() -- fit failed",__FILE__,__LINE__);
    }

    return gaussians_;
  }

  template <typename Function,
  typename Weight>
  void
  GaussianFit<Function,Weight>::extract_params() const
  {
    const unsigned int ngaussians = gaussians_.size();
    unsigned int pp;
    for(unsigned int p=0, pp=0; p<ngaussians; ++p) {
      p_[pp++] = gaussians_[p].first;
      p_[pp++] = gaussians_[p].second;
    }
  }

  template <typename Function,
  typename Weight>
  void
  GaussianFit<Function,Weight>::assign_params() const
  {
    const unsigned int ngaussians = gaussians_.size();
    unsigned int pp;
    for(unsigned int p=0, pp=0; p<ngaussians; ++p) {
      gaussians_[p].first = p_[pp++];
      gaussians_[p].second = p_[pp++];
    }
  }
};

#endif // include guards
/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:

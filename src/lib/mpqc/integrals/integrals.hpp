//
// integrals.hpp
//
// Copyright (C) 2013 Drew Lewis and Andrey Asadchev
//
// Authors: Drew Lewis
// Maintainer: Drew Lewis and Edward Valeev
//
// This file is part of the MPQC Toolkit.
//
// The MPQC Toolkit is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The MPQC Toolkit is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the MPQC Toolkit; see the file COPYING.LIB.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//
// The U.S. Government is granted a limited license as per AL 91-7.
//

#ifndef MPQC_INTEGRALS_INTEGRALS_HPP
#define MPQC_INTEGRALS_INTEGRALS_HPP

#include <vector>

#include <mpqc/math/tensor.hpp>
#include <mpqc/utility/foreach.hpp>

#include <chemistry/molecule/molecule.h>
#include <chemistry/qc/basis/integral.h>

#include "shell.hpp"

namespace mpqc {
namespace integrals {

    /// @addtogroup ChemistryBasisIntegralRange
    /// @{
    /**
     * Wraps an MPQC integral engine (e.g. sc::TwoBodyInt)
     */
    template<class RefEngine>
    class Integrals {
    public:
        typedef TensorRef<const double,2, TensorRowMajor > Tensor2;
        typedef TensorRef<const double,3, TensorRowMajor > Tensor3;
        typedef TensorRef<const double,4, TensorRowMajor > Tensor4;

        /**
         * Constructor for Integrals
         * @param engineptr is a sc::Ref to an MPQC integral engine
         */
        explicit Integrals(RefEngine engineptr)
            : engine_(engineptr) {}

        RefEngine& engine() {
            return engine_;
        }

        /**
         * Calls the MPQC integral object on shells p and q and returns a
         * TensorRef holding the integral buffer.
         */
        Tensor2 operator()(Shell p, Shell q) {
            size_t dims[] = {size_t(p.size()), size_t(q.size()) };
            engine_->compute_shell(p.index(), q.index());
            return Tensor2(engine_->buffer(), dims);
        }

        Tensor3 operator()(Shell p, Shell q, Shell r) {
            size_t dims[] = {size_t( p.size()), size_t(q.size()),
                             size_t(r.size()) };
            engine_->compute_shell(p.index(), q.index(), r.index());
            return Tensor3(engine_->buffer(), dims);
        }

        Tensor4 operator()(Shell p, Shell q, Shell r, Shell s) {
            size_t dims[] = { size_t(p.size()), size_t(q.size()),
                              size_t(r.size()), size_t(s.size()) };
            engine_->compute_shell(p.index(), q.index(),
                                     r.index(), s.index());
            return Tensor4(engine_->buffer(), dims);
        }

    private:
        RefEngine engine_;
    };

    /// @} //ChemistryBasisIntegralRange

  namespace detail {

    /// Determines the number of functions in a given shell and returns a range of Shell objects
    /// for mpqc integrals.
    inline std::vector<Shell> pack(sc::Ref<sc::GaussianBasisSet> basis,
                                   const std::vector<int> &S) {
        int f = 0;
        std::vector<Shell> shells;
        foreach(int s, S){
            int n = basis->shell(s).nfunction();
            shells.push_back(Shell(s,range(f,f+n)));
            f += n;
        }
        return shells;
    }


    //inline size_t extent(const std::vector<Shell> &S) {
    //    size_t extent = 0;
    //    foreach(Shell s, S){
    //        extent = std::max(extent, *s.end());
    //    }
    //    return extent;
    //}



    /// Evaluates shell integrals and places them in a TensorRef.
    template<class Engine>
    void evaluate(Integrals<Engine> integral,
                  const std::vector<int> &P,
                  const std::vector<int> &Q,
                  TensorRef<double,2, TensorRowMajor > &ints) {
        std::vector<Shell> shells[] = {
            pack(integral.engine()->basis1(), P),
            pack(integral.engine()->basis2(), Q)
        };

        foreach(Shell p, shells[0]){
            foreach(Shell q, shells[1]){
                ints(p,q) = integral(p,q);
            }
        }

    }

    template<class Engine>
    void evaluate(Integrals<Engine> integral,
                  const std::vector<int> &P,
                  const std::vector<int> &Q,
                  const std::vector<int> &R,
                  TensorRef<double,3, TensorRowMajor > &ints) {
        std::vector<Shell> shells[] = {
            pack(integral.engine()->basis1(), P),
            pack(integral.engine()->basis2(), Q),
            pack(integral.engine()->basis3(), R)
        };
        foreach(Shell p, shells[0]){
            foreach(Shell q, shells[1]){
                foreach(Shell r, shells[2]){
                    ints(p,q,r) = integral(p,q,r);
                }
            }
        }
    }

    template<class Engine>
    void evaluate(Integrals<Engine> integral,
                  const std::vector<int> &P,
                  const std::vector<int> &Q,
                  const std::vector<int> &R,
                  const std::vector<int> &S,
                  TensorRef<double,4, TensorRowMajor > &ints) {
        std::vector<Shell> shells[] = {
            pack(integral.engine()->basis1(), P),
            pack(integral.engine()->basis2(), Q),
            pack(integral.engine()->basis3(), R),
            pack(integral.engine()->basis4(), S),
        };
        foreach(Shell p, shells[0]){
            foreach(Shell q, shells[1]){
                foreach(Shell r, shells[2]){
                    foreach(Shell s, shells[3]){
                        ints(p,q,r,s) = integral(p,q,r,s);
                    }
                }
            }
        }
    }

} // namespace detail
} // namespace integrals
} // namespace mpqc


namespace mpqc {
namespace integrals {

    /// @addtogroup ChemistryBasisIntegralRange
    /// @{

    /**
       Evaluate set of shell blocks of integrals (p|O|q)
       @param[in] engine integral engine
       @param[in] P list of p shell indices
       @param[in] Q list of q shell indices
       @param[out] (p|O|q) integrals
     */
    inline void evaluate(sc::Ref<sc::OneBodyInt> &engine,
                         const std::vector<int> &P,
                         const std::vector<int> &Q,
                         TensorRef<double,2, TensorRowMajor > &ints) {
      detail::evaluate(Integrals<sc::Ref<sc::OneBodyInt> >(engine), P, Q,
                       ints);
    }

    inline void evaluate(sc::Ref<sc::TwoBodyTwoCenterInt> &engine,
                         const std::vector<int> &P,
                         const std::vector<int> &Q,
                         TensorRef<double,2, TensorRowMajor > &ints) {
      detail::evaluate(Integrals<sc::Ref<sc::TwoBodyTwoCenterInt> >(engine),
                       P, Q, ints);
    }

    inline void evaluate(sc::Ref<sc::TwoBodyThreeCenterInt> &engine,
                         const std::vector<int> &P,
                         const std::vector<int> &Q,
                         const std::vector<int> &R,
                         TensorRef<double,3, TensorRowMajor > &ints) {
      detail::evaluate(Integrals<sc::Ref<sc::TwoBodyThreeCenterInt> >(engine),
                       P, Q, R, ints);
    }

    inline void evaluate(sc::Ref<sc::TwoBodyInt> engine,
                         const std::vector<int> &P,
                         const std::vector<int> &Q,
                         const std::vector<int> &R,
                         const std::vector<int> &S,
                         TensorRef<double,4, TensorRowMajor > &ints) {
      detail::evaluate(Integrals<sc::Ref<sc::TwoBodyInt> >(engine),
                       P, Q, R, S, ints);
    }


    /// @} // ChemistryBasisIntegralRange
} // namespace integrals
} // namespace mpqc


#endif /* MPQC_INTEGRALS_INTEGRALS_HPP */

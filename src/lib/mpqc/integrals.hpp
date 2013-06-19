#ifndef MPQC_INTEGRALS_HPP
#define MPQC_INTEGRALS_HPP

#include <vector>
#include <boost/foreach.hpp>

#include "mpqc/range.hpp"
#include "mpqc/tensor.hpp"

#include <chemistry/molecule/molecule.h>
#include <chemistry/qc/basis/integral.h>

namespace mpqc {
namespace integrals {
namespace detail {


    struct Shell : range {
        Shell(int index, range r) : range(r), index_(index) {}
        int index() const {
            return index_;
        }
    private:
        int index_;
    };

//     struct Screening {
//         Screening()
//         bool test(const Shell &p, const Shell &q,
//                   const Shell &r, const Shell &s) {
//         }
//     private:
//         if (tbi.log2_shell_bound(i,j,k,l)+pmaxijkl < tol)
// //                 continue;
//         std::vector<char> pmax_;
//     };

    inline std::vector<Shell> pack(sc::Ref<sc::GaussianBasisSet> basis,
                                   const std::vector<int> &S) {
        int f = 0;
        std::vector<Shell> shells;
        BOOST_FOREACH (int s, S) {
            int n = basis->shell(s).nfunction();
            shells.push_back(Shell(s, range(f, f+n)));
            f += n;
        }
        return shells;
    }


    inline size_t extent(const std::vector<Shell> &S) {
        int extent = 0;
        BOOST_FOREACH (Shell s, S) {
            extent = std::max(extent, *s.end());
        }
        return extent;
    }


    template<class Engine>
    struct Integral {

        typedef integrals::detail::Shell Shell;
        typedef tensor<const double,2> Tensor2;
        typedef tensor<const double,3> Tensor3;
        typedef tensor<const double,4> Tensor4;

        explicit Integral(sc::Ref<Engine> integral)
            : integral_(integral) {}

        sc::Ref<Engine>& engine() {
            return integral_;
        }

        int max2(Shell p, Shell q, Shell r, Shell s) const {
            int i = p.index();
            int j = q.index();
            int k = r.index();
            int l = s.index();
            return integral_->log2_shell_bound(i,j,k,l);
        }

        Tensor2 operator()(Shell p, Shell q) {
            size_t dims[] = { p.size(), q.size() };
            integral_->compute_shell(p.index(), q.index());
            return Tensor2(integral_->buffer(), dims);
        }

        Tensor3 operator()(Shell p, Shell q, Shell r) {
            size_t dims[] = { p.size(), q.size(), r.size() };
            integral_->compute_shell(p.index(), q.index(), r.index());
            return Tensor3(integral_->buffer(), dims);
        }

        Tensor4 operator()(Shell p, Shell q, Shell r, Shell s) {
            size_t dims[] = { p.size(), q.size(), r.size(), s.size() };
            integral_->compute_shell(p.index(), q.index(),
                                     r.index(), s.index());
            return Tensor4(integral_->buffer(), dims);
        }

    private:
        sc::Ref<Engine> integral_;
    };

    template<class Engine>
    void evaluate(Integral<Engine> integral,
                  const std::vector<int> &P,
                  const std::vector<int> &Q,
                  tensor<double,2> &ints) {
        std::vector<Shell> shells[] = {
            pack(integral.engine()->basis1(), P),
            pack(integral.engine()->basis2(), Q)
        };
        BOOST_FOREACH (Shell p, shells[0]) {
            BOOST_FOREACH (Shell q, shells[1]) {
                ints(p,q) = integral(p,q);
            }
        }
    }

    template<class Engine>
    void evaluate(Integral<Engine> integral,
                  const std::vector<int> &P,
                  const std::vector<int> &Q,
                  const std::vector<int> &R,
                  tensor<double,3> &ints) {
        std::vector<Shell> shells[] = {
            pack(integral.engine()->basis1(), P),
            pack(integral.engine()->basis2(), Q),
            pack(integral.engine()->basis3(), R)
        };
        BOOST_FOREACH (Shell p, shells[0]) {
            BOOST_FOREACH (Shell q, shells[1]) {
                BOOST_FOREACH (Shell r, shells[2]) {
                    ints(p,q,r) = integral(p,q,r);
                }
            }
        }
    }

    template<class Engine>
    void evaluate(Integral<Engine> integral,
                  const std::vector<int> &P,
                  const std::vector<int> &Q,
                  const std::vector<int> &R,
                  const std::vector<int> &S,
                  tensor<double,4> &ints) {
        std::vector<Shell> shells[] = {
            pack(integral.engine()->basis1(), P),
            pack(integral.engine()->basis2(), Q),
            pack(integral.engine()->basis3(), R),
            pack(integral.engine()->basis4(), S),
        };
        BOOST_FOREACH (Shell p, shells[0]) {
            BOOST_FOREACH (Shell q, shells[1]) {
                BOOST_FOREACH (Shell r, shells[2]) {
                    BOOST_FOREACH (Shell s, shells[3]) {
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

    /**
       Evaluate list of shell integrals (p|O|q)
       @param engine integral engine
       @param P list of p shell indices
       @param Q list of q shell indices
       @param[out] (p|O|q) integrals
     */
    void evaluate(sc::Ref<sc::OneBodyInt> engine,
                  const std::vector<int> &P,
                  const std::vector<int> &Q,
                  tensor<double,2> &ints) {
        detail::evaluate(detail::Integral<sc::OneBodyInt>(engine), P, Q, ints);
    }

    void evaluate(sc::Ref<sc::TwoBodyTwoCenterInt> engine,
                  const std::vector<int> &P,
                  const std::vector<int> &Q,
                  tensor<double,2> &ints) {
        detail::evaluate(detail::Integral<sc::TwoBodyTwoCenterInt>(engine),
                         P, Q, ints);
    }

    void evaluate(sc::Ref<sc::TwoBodyThreeCenterInt> engine,
                  const std::vector<int> &P,
                  const std::vector<int> &Q,
                  const std::vector<int> &R,
                  tensor<double,3> &ints) {
        detail::evaluate(detail::Integral<sc::TwoBodyThreeCenterInt>(engine),
                         P, Q, R, ints);
    }

    void evaluate(sc::Ref<sc::TwoBodyInt> engine,
                  const std::vector<int> &P,
                  const std::vector<int> &Q,
                  const std::vector<int> &R,
                  const std::vector<int> &S,
                  tensor<double,4> &ints) {
        detail::evaluate(detail::Integral<sc::TwoBodyInt>(engine),
                         P, Q, R, S, ints);
    }

} // namespace integrals
} // namespace mpqc


#endif /* MPQC_INTEGRALS_HPP */

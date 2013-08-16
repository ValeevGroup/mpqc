#ifndef MPQC_INTEGRALS_INTEGRALS_HPP
#define MPQC_INTEGRALS_INTEGRALS_HPP

#include <vector>

#include <mpqc/math/tensor.hpp>
#include <mpqc/utility/foreach.hpp>

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


    inline size_t extent(const std::vector<Shell> &S) {
        int extent = 0;
        foreach(Shell s, S){
            extent = std::max(extent, *s.end());
        }
        return extent;
    }


    template<class RefEngine>
    class Integrals {

    public:
        typedef TensorRef<const double,2, TensorRowMajor > Tensor2;
        typedef TensorRef<const double,3, TensorRowMajor > Tensor3;
        typedef TensorRef<const double,4, TensorRowMajor > Tensor4;

        explicit Integrals(RefEngine integral)
            : integral_(integral) {}

        RefEngine& engine() {
            return integral_;
        }

        Tensor2 operator()(Shell p, Shell q) {
            size_t dims[] = {size_t(p.size()), size_t(q.size()) };
            integral_->compute_shell(p.index(), q.index());
            return Tensor2(integral_->buffer(), dims);
        }

        Tensor3 operator()(Shell p, Shell q, Shell r) {
            size_t dims[] = {size_t( p.size()), size_t(q.size()),
                             size_t(r.size()) };
            integral_->compute_shell(p.index(), q.index(), r.index());
            return Tensor3(integral_->buffer(), dims);
        }

        Tensor4 operator()(Shell p, Shell q, Shell r, Shell s) {
            size_t dims[] = { size_t(p.size()), size_t(q.size()),
                              size_t(r.size()), size_t(s.size()) };
            integral_->compute_shell(p.index(), q.index(),
                                     r.index(), s.index());
            return Tensor4(integral_->buffer(), dims);
        }

    private:
        RefEngine integral_;
    };

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

    /**
       Evaluate list of shell integrals (p|O|q)
       @param engine integral engine
       @param P list of p shell indices
       @param Q list of q shell indices
       @param[out] (p|O|q) integrals
     */
    inline void evaluate(sc::Ref<sc::OneBodyInt> &engine,
                  const std::vector<int> &P,
                  const std::vector<int> &Q,
                  TensorRef<double,2, TensorRowMajor > &ints) {
        detail::evaluate(detail::Integrals<sc::Ref<sc::OneBodyInt> >(engine), P, Q, ints);
    }

    inline void evaluate(sc::Ref<sc::TwoBodyTwoCenterInt> &engine,
                  const std::vector<int> &P,
                  const std::vector<int> &Q,
                  TensorRef<double,2, TensorRowMajor > &ints) {
        detail::evaluate(detail::Integrals<sc::Ref<sc::TwoBodyTwoCenterInt> >(engine),
                         P, Q, ints);
    }

    inline void evaluate(sc::Ref<sc::TwoBodyThreeCenterInt> &engine,
                  const std::vector<int> &P,
                  const std::vector<int> &Q,
                  const std::vector<int> &R,
                  TensorRef<double,3, TensorRowMajor > &ints) {
        detail::evaluate(detail::Integrals<sc::Ref<sc::TwoBodyThreeCenterInt> >(engine),
                         P, Q, R, ints);
    }

    inline void evaluate(sc::Ref<sc::TwoBodyInt> engine,
                  const std::vector<int> &P,
                  const std::vector<int> &Q,
                  const std::vector<int> &R,
                  const std::vector<int> &S,
                  TensorRef<double,4, TensorRowMajor > &ints) {
        detail::evaluate(detail::Integrals<sc::Ref<sc::TwoBodyInt> >(engine),
                         P, Q, R, S, ints);
    }


} // namespace integrals
} // namespace mpqc


#endif /* MPQC_INTEGRALS_INTEGRALS_HPP */

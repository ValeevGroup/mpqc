#ifndef MPQC_CI_SIGMA3_HPP
#define MPQC_CI_SIGMA3_HPP

#include "mpqc/ci/ci.hpp"
#include "mpqc/ci/string.hpp"
#include "mpqc/ci/subspace.hpp"

#include "mpqc/utility/timer.hpp"
#include "mpqc/range.hpp"
#include "mpqc/math/matrix.hpp"
#include "mpqc/omp.hpp"

#include "mpqc/array.hpp"
#include "mpqc/array/functions.hpp"

#include <boost/type_traits/is_same.hpp>

namespace mpqc {
namespace ci {

    /// One-particle excitation from string I to string J, {sign, ij, I, J}
    struct Excitation {
        float sgn; /** replacement parity, -1 or 1 */
        int integral; /** integral index, (i**2+i)/2 + j */
        int I; /** ref. string index */
        int J; /** exc. string index */
    };

    // inline bool operator<(const Excitation &a, const Excitation &b) {
    //     // sort by I or by J if I's are the same
    //     return ((a.I != b.I) ? (a.I < b.I) : (a.J < b.J));
    // }

    inline std::ostream& operator<<(std::ostream &os, const Excitation &r) {
        os << r.I << " -> " << r.J
           << ", phase = " << r.sgn
           << ", ij = " << r.integral;
        return os;
    }

    /// Evaluate valid single excitations from string I into subspace S
    /// and append results to excitations vector.
    template<class CI, class Spin>
    void append(const CI &ci,
                const String &I, const Subspace<Spin> &S,
                std::vector<Excitation> &V) {
        // this particular traversal is used to generate semi-ordered list
        int idx = (int)ci.template strings<Spin>()[I];
        size_t i = I.size();
        while (i--) {
            if (!I[i]) continue; // empty orbital
            for (size_t j = 0; j < I.size(); ++j) {
                if (I[j] && (i != j)) continue; // not an empty orbital 
                String J = I.swap(i, j);
                if (!ci.test(J)) continue;
                int jdx = (int)ci.template strings<Spin>()[J];
                if (!S.test(jdx)) continue; // string not in S
                Excitation t;
                t.sgn = (float)sgn(I, i, j);
                t.integral = (int)index(i,j);
                t.I = idx;
                t.J = jdx;
                V.push_back(t);
            }
        }
    }

    /// Vector of single particle excitations from I to J subspace
    template<class Spin>
    struct Excitations {
        template<class CI>
        Excitations(const CI &ci, const Subspace<Spin> &I, const Subspace<Spin> &J)
            : I_(I), J_(J)
        {
            BOOST_FOREACH (int i, I) {
                append(ci, ci.template strings<Spin>()[i], J, data_);
            }
        }
        typedef std::vector<Excitation>::const_iterator const_iterator;
        const_iterator begin() const { return data_.begin(); }
        const_iterator end() const { return data_.end(); }
        size_t size() const { return data_.size(); }
        const Excitation& operator[](int i) const { return data_.at(i); }
        const Subspace<Spin>& I() const { return this->I_; }
        const Subspace<Spin>& J() const { return this->J_; }
    private:
        std::vector<Excitation> data_;
        Subspace<Spin> I_, J_;
    };

#define MPQC_CI_SIGMA3_NAIVE
#ifdef MPQC_CI_SIGMA3_NAIVE

    /// Naive non-vectorizable sigma3 kernel
    void sigma3(const Excitations<Alpha> &alpha, const Excitations<Beta> &beta,
                const Matrix &V, const Matrix &C, Matrix &S) {                
        // beta->beta excitations
        BOOST_FOREACH (auto b, beta) {
            // alpha->alpha excitations
            BOOST_FOREACH (auto a, alpha) {
                int Ia = a.I;
                int Ja = a.J;
                int Ib = b.I;
                int Jb = b.J;
                double c = C(Ja - *alpha.J().begin(), Jb - *beta.J().begin());
                double s = a.sgn*b.sgn*V(a.integral,b.integral)*c;
                S(Ia - *alpha.I().begin(), Ib - *beta.I().begin()) += s;
            }
        }
    }

#else // MPQC_CI_SIGMA3_NAIVE

    /// Compute s = phase*v*c
    void sigma3(int n, double phase,
                const double * RESTRICT v,
                const double * RESTRICT c,
                double * RESTRICT s) {
        for (int i = 0; i < n; ++i) {
            s[i] += phase*c[i]*v[i];
        }
    }

    /// Compute s = phase*v*c
    template<class V, class C, class S>
    void sigma3(int n, double phase, const V &v, const C &c, S &s) {
        sigma3(n, phase, v.data(), c.data(), s.data());
        // for (int i = 0; i < n; ++i) {
        //     s[i] += phase*c[i]*v[i];
        // }
    }

    /// Vectorized sigma3 kernel
    /// @param A Ia->Ja excitations
    /// @param B Ib->Jb excitations
    /// @param V V(ij,kl)
    /// @param C C(Ja,Jb) block
    /// @param S S(Ia,Ib) block
    /// @param mutex mutex to lock S updates
    template<class SpinA, class SpinB>
    void sigma3(const Excitations<SpinA> &A, const Excitations<SpinB> &B,
                const Matrix &V, const Matrix &C, Matrix &S) {
        static_assert(!boost::is_same<SpinA,SpinB>::value, "spins must not be the same");
        const int BLOCK = 128;
        mpqc::Matrix c, v;
        mpqc::Vector s;
        std::vector<int> Ia;
        // iterate over A excitations in blocks
        for (int ba = 0; ba < A.size(); ba += BLOCK) {
            int na = std::min<int>(A.size() - ba, BLOCK);
            c.resize(na, C.cols());
            v.resize(na, V.cols());
            Ia.clear();
            /// pack the Ia->Ja dimension
            for (int a = 0; a < na; ++a) {
                auto aa = A[a+ba];
                c.row(a) = C.row(aa.J - *A.J().begin());
                v.row(a) = V.col(aa.integral);
                v.row(a) *= aa.sgn;
                Ia.push_back(aa.I - *A.I().begin()); // record Ia index
            }
            s = mpqc::Vector::Zero(na);
            for (auto b = B.begin(); b != B.end(); ++b) {
                int Ib = b->I - *B.I().begin();
                int Jb = b->J - *B.J().begin();
                // vectorized s += sgn(Ib,Jb)*v*c
                sigma3(na, b->sgn, v.col(b->integral), c.col(Jb), s);
                // next iteration updates same Ib
                if ((b+1) != B.end() && (b+1)->I == (b)->I) continue;
                // scatter
                //mutex.lock();
                for (int a = 0; a < Ia.size(); ++a) {
                    S(Ia[a], Ib) += s(a);
                    s(a) = 0;
                }
                //mutex.unlock();
            }
        }
    }

#endif // MPQC_CI_SIGMA3_NAIVE

}
}

#endif /* MPQC_CI_SIGMA3_HPP */


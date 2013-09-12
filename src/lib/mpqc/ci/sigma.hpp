#include "mpqc/ci/ci.hpp"
#include "mpqc/ci/string.hpp"

#include "mpqc/utility/timer.hpp"
#include "mpqc/range.hpp"
#include "mpqc/math/matrix.hpp"
#include "mpqc/omp.hpp"

#include "mpqc/array.hpp"
#include "mpqc/array/functions.hpp"

#define MPQC_PROFILE_ENABLE
#include "mpqc/utility/profile.hpp"

namespace mpqc {
namespace ci {

    /// Computes parity, 1 for even parity, -1 for odd
    int sgn(size_t ij) {
        return (ij % 2) ? -1 : 1;
    }

    /// Compute parity of exciting i to j (or j to i).
    /// Counts number of set bits in open interval (i:j)
    /// and returns parity of the pop count.
    inline int sgn(const String &I, int i, int j) {
        uint64_t b = I.to_ulong();
        int n = std::abs(i - j);
        if (j < i)
            std::swap(i, j);
        b = b << (63 - j);
        b = b << 1;
        b = b >> 2;
        b = b >> (63 - n);
        size_t p = String::bitset(b).count();
        //size_t p = _mm_popcnt_u64(b);
#ifndef NDEBUG
        assert(p < I.count());
        size_t q = 0;
        for (int k = std::min(i,j)+1; k < std::max(i,j); ++k) {
            q += I[k];
        }
        // printf("string %s(%lu) [%i,%i] b=%lu, p=%lu, q=%lu\n",
        //        std::string(I).c_str(), I.to_ulong(), i,j, b, p, q);
        assert(p == q);
#endif
        return sgn(p);
    }

    /**
     * a set of 1-particle replacements (replacements are particle-number
     *  conserving operators that change occupation numbers).
     */
    struct Replacement {
        /**
         * a one-particle replacement = {sign, ij, string}, where for replacement I->J
         * IJ = (I > J) ? I*(I+1)/2 + J : J*(J+1)/2 + I.
         */
        struct Tuple {
            float sgn; /** replacement parity, -1 or 1 */
            int integral; /** integral index, (i**2+i)/2 + j */
            int I; /** ref. string index */
            int J; /** exc. string index */
            bool operator<(const Tuple &t) const {
                return (this->J < t.J);
            }
        };
        typedef std::vector<Tuple> Single;
    };

    /**
     * produces all single replacements from argument String list[s]
     * @param[in] list the set of Strings
     * @param[in] s the index of the argument String in list
     * @param[out] r set of single replacements
     */
    template<class L, class CI>
    void replacements(const String::List<L> &list, size_t s,
                      Replacement::Single &r, const CI &ci) {
        const String &I = list.at(s);
        // this particular traversal is used to generate semi-ordered list
        size_t i = I.size();
        size_t pos = r.size();
        while (i--) {
            if (!I[i])
                continue;
            for (size_t j = 0; j < I.size(); ++j) {
                if (I[j] && (i != j))
                    continue;
                String J = I.swap(i, j);
                if (!ci.test(J))
                    continue;
                Replacement::Tuple t;
                t.sgn = (float) sgn(I, i, j);
                t.integral = (int) index(j, i);
                t.I = (int) s;
                t.J = (int) list[J];
                r.push_back(t);
                //std::cout << J << "@" << t.J << std::endl;
            }
            //std::cout << "---" << std::endl;
        }
    }

    template<class L, class CI>
    Replacement::Single replacements(const String::List<L> &list, size_t s,
                                     const CI &ci) {
        Replacement::Single r;
        replacements(list, s, r, ci);
        return r;
    }

    template<class L, class CI>
    Replacement::Single replacements(const String::List<L> &list, size_t begin,
                                     size_t end, const CI &ci) {
        Replacement::Single S, M;
        for (size_t s = begin; s < end; ++s) {
            size_t pos = S.size();
            replacements(list, s, S, ci);
        }
        return S;
    }

    template<class L, class F_, class Type>
    size_t sigma2(const String &I, const String::List<L> &list,
                  const mpqc::Vector &h,
                  const mpqc::Matrix &V,
                  const Eigen::MatrixBase<F_> &f,
                  const CI<Type> &ci) {
        auto &F = const_cast<Eigen::MatrixBase<F_>&>(f);
        size_t ops = 0;
        size_t count = I.count();
        std::vector<int> O, E(1); // reserve k->k excitation
        for (int l = 0; l < I.size(); ++l) {
            if (I[l])
                O.push_back(l); // occ. orbs
            if (!I[l])
                E.push_back(l); // exc. orbs
        }
        //asm("#andrey");
        for (auto o = O.begin(); o < O.end(); ++o) {
            int k = *o;
            E[0] = k; // k->k
            for (auto e = E.begin(); e < E.end(); ++e) {
                int l = *e;
                String J = I.swap(k, l);
                if (!ci.test(J))
                    continue;
                double sgn_kl = sgn(I, k, l);
                int kl = index(k, l);
                F(list[J]) += sgn_kl * h(kl);
                std::swap(*e, *o); // l->k
                foreach (int i, O) {
                    E[0] = i;
                    //for (int j = 0; j < J.size(); ++j) {
                    //if (J[j] && j != i) continue;
                    foreach (int j, E) {
                        String K = J.swap(i,j);
                        if (!ci.test(K)) continue;
                        F(list[K]) += 0.5*sgn_kl*sgn(J,i,j)*V(index(i,j),kl);
                    }
                }
                // restore original vectors
                E[0] = k;
                std::swap(*e, *o);
                ops += O.size() * E.size() * 3;
            }
        }
        //for (int i = 0; i < F.size(); ++i) 
        //printf("F(%i) = %e\n", i, F(i));
        return ops;
    }

    template<class Vector>
    size_t sigma3(const Replacement::Single &A, const std::vector<int> &index,
                  const Matrix &V, const Vector &C, Matrix &S) {
        size_t ops = 0;
        foreach (auto const &ij, A) {
            double cij = ij.sgn*C(ij.J);
            if (fabs(cij) < 1e-5) continue;
            auto s = S.row(ij.I);
            auto v = V.col(ij.integral);
            //#pragma vector aligned
            for (int i = 0; i < index.size(); ++i) {
                s(index[i]) += cij*v(i);
            }
            ++ops;
        }
        return 2 * ops * index.size();
    }

    void sigma3(int n, const double *v, const double *c, int phase,
                double* __restrict__ s) {
        assert((phase == 1) || (phase == -1));
        if (phase == 1) {
            for (int i = 0; i < n; ++i) {
                s[i] += v[i] * c[i];
            }
        }
        if (phase == -1) {
            for (int i = 0; i < n; ++i) {
                s[i] -= v[i] * c[i];
            }
        }
    }

    size_t sigma3(const Replacement::Single &A, const Matrix &V,
                  const Matrix &C, Matrix &S, size_t offset = 0) {
        foreach (auto const &ij, A) {
            sigma3(S.rows(),
                   V.col(ij.integral).data(),
                   C.col(ij.J).data(),
                   ij.sgn,
                   S.col(ij.I-offset).data());
        }
        return 2 * A.size() * S.rows();
    }

    /**
     * computes sigma3 contribution to a block of sigma:
     *
     * \f$ \sigma_3(I_a, I_b) =
     *     \sum_{J_b \in R^1(I_b)}
     *     \mathrm{sign}(I_b,J_b) \sum_{J_a \in R^1(I_a)} \mathrm{sign}(I_a,J_a)
     *     V(\mathrm{ind}(I_a,J_a), \mathrm{ind}(I_b,J_b) ) C(J_a, J_b )
     * \f$
     * <blockquote>
     * the algorithm:\n
     *   foreach \f$ I_b \in \{ I_b \} \f$ \n
     *    compute \f$ \{ R^1(I_b) \} = \{ J_b \} \f$ \n
     *    foreach block of \f$ J_b \f$ \n
     *      copy block of \f$ C(J_a, J_b)  \f$ \n
     *      compute \f$ \tilde{V} (..., J_b) = V ( ... , \mathrm{ind}(I_b,J_b) )
     *                       \mathrm{sign}(I_b,J_b) \f$ \n
     *      foreach block of \f$ I_a \in \{ I_a \} \f$ \n
     *        compute \f$ \{ R^1(I_a) \} = \{ J_a \} \f$ foreach \f$ I_a \f$ in the block \n
     *        *** calling sigma3() *** \n
     *        foreach \f$ I_a \f$ in the \f$ I_a \f$ block \n
     *          foreach \f$ J_a \f$ in \f$ \{ R^1(I_a) \} \f$ \n
     *            s = 0 \n
     *            foreach \f$ J_b \f$ in the \f$ J_b \f$ block \n
     *              s += \f$ C(J_a, J_b) * \tilde{V} ( \mathrm{ind}(I_a,J_a) , J_b) \f$ \n
     *            endfor \n
     *            s *= \f$ \mathrm{sign}(I_a,J_a) \f$ \n
     *          endfor \n
     *          \f$ \sigma_3(I_a,I_b) += s \f$ \n
     *        endfor \n
     *      endfor \n
     *    endfor \n
     *  endfor
     *  </blockquote>
     *
     * @param[in] alpha the list of alpha-spin Strings in the target block
     * @param[in] beta the list of beta-spin Strings in the target block
     * @param[in] h' core Hamiltonian matrix not used right now,
     \f$ h'_{kl} = h_{kl} + 0.5 \sum_j (kj|jl) \f$
     * @param[in] V Coulomb integrals matrix
     * @param[in] C CI coefficients
     * @param[out] S the result matrix; the contribution is accumulated
     */
    template<class Index>
    void sigma(const CI< Full<Index> > &ci,
               const mpqc::Vector &h, const Matrix &V,
               const mpqc::Array<double> &C,
               mpqc::Array<double> &S) {

        const auto &alpha = ci.alpha;
        const auto &beta = ci.beta;
        auto &comm = ci.comm;

        size_t ne = alpha[0].count();
        size_t no = alpha[0].size();

        struct {
            double sigma1, sigma2, kernel12;
            double sigma3, kernel3;
            size_t ops;
        } time = { };

        struct {
            std::vector<range> ranges;
            std::vector<Replacement::Single> list;
        } AA;

        // generate and sort aa lists
        {
            size_t block = ci.block2; //64;
            for (size_t A = 0; A < alpha.size(); A += block) {
                size_t na = std::min(block, alpha.size() - A);
                AA.ranges.push_back(range(A, A + na));
                AA.list.push_back(replacements(alpha, A, A + na, ci));
                std::sort(AA.list.back().begin(), AA.list.back().end());
            }
        }

        {
            timer t;
            mpqc::Vector H = h;
            for (int j = 0; j < no; ++j) {
                for (int i = 0; i <= j; ++i) {
                    double v = 0;
                    for (int k = 0; k < no; ++k) {
                        v += V(index(i, k), index(k, j));
                    }
                    H(index(i, j)) -= 0.5 * v;
                }
            }

            foreach (auto rj, range(ci.local.beta).block(ci.block)) {
                MPQC_PROFILE_LINE;

                mpqc::Matrix s, c;
                {
                    MPQC_PROFILE_LINE;
                    c.resize(alpha.size(), rj.size());
                    C(alpha,rj) >> c;
                    c = Matrix(c).transpose();
                }

                s.resize(rj.size(), alpha.size());
                s.fill(0);

                {
                    MPQC_PROFILE_LINE;
                    timer t;
#pragma omp parallel for schedule(static,1)
                    for (int i = 0; i < alpha.size(); ++i) {
                        mpqc::Vector F(alpha.size());
                        F.fill(0);
                        size_t ops = ci::sigma2(alpha[i], alpha, H, V, F, ci);
                        for (int k = 0; k < alpha.size(); ++k) {
                            double f = F(k);
                            if (fabs(f) < 1e-14) continue;
                            s.col(i) += f*c.col(k);
                        }
                    }
                    time.kernel12 += t;
                }

                {
                    S(alpha,rj) << Matrix(s.transpose());
                    //S(alpha,rj) << s;
                }

            }

            S.sync();

            time.sigma1 = t;
            t.reset();

            int ms = 0;
            if (ms == 0) {
                size_t BLOCK = sqrt(ci.block * alpha.size());
                symmetrize(S, ci.comm, BLOCK);
            }

            S.sync();

            time.sigma2 = t;
            // std::cout << "sigma1+2 time: " << time.sigma1+time.sigma2 << std::endl;
            // std::cout << "sigma1+2 kernel time: " << time.kernel12 << std::endl;
        }

        MPI::Task task(ci.comm);

        while (true) {
            timer t;

            range rb = task.next(beta, ci.block);
            if (!rb.size())
                break;

            Matrix s(alpha.size(), rb.size());
            S(alpha, rb) >> s;

            // ci.comm.printf("rb = [%i,%i]\n", rb.begin(), rb.end());

            // all beta->beta single replacements in rb, sorted
            auto BB = replacements(beta, *rb.begin(), *rb.end(), ci);
            std::sort(BB.begin(), BB.end());

            // unique J indices, sorted
            std::vector<int> Jb;
            {
                std::set<int> J;
                foreach (auto bb, BB) {
                    J.insert(bb.J);
                }
                Jb.resize(J.size());
                std::copy(J.begin(), J.end(), Jb.begin());
                std::sort(Jb.begin(), Jb.end());
            }

            int BLOCK = 128;

            // loop over unique J indices, block at a time
#pragma omp parallel for schedule(dynamic,1)
            for (size_t jb = 0; jb < Jb.size(); jb += BLOCK) {
                size_t nb = std::min<int>(BLOCK, Jb.size() - jb);

                Matrix c(alpha.size(), nb);
                Matrix ct, vt;

                timer t;

                Replacement::Single IJ;
                std::map<int, int> Jb_index;
                std::vector<int> J_index;

                // load J block
                {
                    for (size_t j = 0; j < nb; ++j) {
                        int J = Jb.at(j + jb);
                        foreach (const auto &bb, BB) {
                            if (bb.J == J) IJ.push_back(bb);
                        }
                        Jb_index[J] = j;
                        J_index.push_back(J);
                        //std::cout << "J=" << J << std::endl;
                        // MPQC_PROFILE_LINE;
                        // c.col(j) = (const Vector&)C(alpha, J);
                    }
                }

                // aggregate read
                {
                    for (size_t j = 0; j < nb;) {
                        int J = J_index[j];
                        int n = 1;
                        while (true) {
                            if (j + n == nb)
                                break;
                            if (J + n != J_index[j + n])
                                break;
                            ++n;
                        }
                        //std::cout << "RJ = " << range(J,J+n) << std::endl;
                        const Matrix &cj = C(alpha, range(J, J + n));
                        c(alpha, range(j, j + n)) = cj;
                        //C(alpha,range(J,J+n)) >> c(alpha,range(j,j+n));
                        j += n;
                    }
                }

                //ci.comm.printf("C reading took %f s\n", (double)t);

                t.reset();

                // loop over excitations into J block
                for (size_t b = 0; b < IJ.size(); b += BLOCK) {
                    size_t nij = std::min<int>(BLOCK, IJ.size() - b);

                    ct.resize(nij, alpha.size());
                    vt.resize(nij, V.cols());
            
                    {
                        std::vector<int> index;
                        for (int i = 0; i < nij; ++i) {
                            const auto &ij = IJ.at(i + b);
                            index.push_back(Jb_index[ij.J]);
                            vt.row(i) = V.col(ij.integral);
                            //ct.row(i) = c.col(Jb_index[ij.J]);
                        }

                        int N = 64;
                        Matrix t;
                        // pack C' matrix
                        for (int a = 0; a < alpha.size(); a += N) {
                            range ra(a, std::min<int>(a + N, alpha.size()));
                            t.resize(ra.size(), index.size());
                            for (int i = 0; i < index.size(); ++i) {
                                t.col(i) = c(ra, index[i]);
                            }
                            range ri(0, index.size());
                            ct(ri, ra) = t.transpose();
                        }
                    }

                    timer t;
                    size_t ops = 0;
                    {
                        Matrix s1, s2;
                        for (int a = 0; a < AA.list.size(); ++a) {
                            range ra = AA.ranges[a];
                            s1.resize(nij, ra.size());
                            s1.fill(0);
                            s2.resize(ra.size(), nij);
                            ops += sigma3(AA.list[a], vt, ct, s1, *ra.begin());
                            for (int i = 0; i < nij; ++i) {
                                const auto &ij = IJ.at(i + b);
                                s2.col(i) = ij.sgn *mpqc::Vector(s1.row(i).transpose());
                            }
#pragma omp critical
                            for (int i = 0; i < nij; ++i) {
                                const auto &ij = IJ.at(i + b);
                                s(ra, ij.I - *rb.begin()) += s2.col(i);
                            }
                        }
                    }

                    //#pragma omp master
                    if (omp::master()) {
                        time.kernel3 += t;
                        time.ops += ops;
                    }

                }

                //ci.comm.printf("S computing took %f\n", (double)t);

            }

            S(alpha, rb) << s;
            time.sigma3 += t;

        }

        S.sync();

        // std::cout << "sigma3 time: " << time.sigma3 << std::endl;
        // std::cout << "sigma3 kernel time: " << time.kernel3 << std::endl;
        // std::cout << "sigma3 flops: " << time.ops/time.kernel3 << std::endl;

    }

}
}

#ifndef MPQC_CI_INTEGRALS_HPP
#define MPQC_CI_INTEGRALS_HPP

#include "mpqc/math/matrix.hpp"
#include "mpqc/utility/foreach.hpp"

#include <chemistry/qc/basis/tbint.h>

namespace mpqc {
  namespace ci {

    template<class C, class U, class T>
    void transform(const range &p, const C &c, const U &u, T &t) {
      t += u.transpose() * c(range(c.rows()), p).transpose();
    }

    void integrals(sc::Ref<sc::GaussianBasisSet> basis, const Matrix &C,
                   sc::Ref<sc::TwoBodyInt> int2, Matrix &V) {
      size_t no = C.rows();
      range O(0, no);
      range shells(basis->nshell());

      V.resize((no * no + no) / 2, (no * no + no) / 2);
      V.fill(0);

      Matrix T1, T2, T3, T4;
      T4.resize(no * no * no, no);
      T4.fill(0);

      foreach (auto s, shells) {
        range S = basis->range(s);
        int ns = S.size();
        T3.resize(ns*no*no, no);
        T3.fill(0);
        foreach (auto r, shells) {
          range R = basis->range(r);
          int nr = R.size();
          T2.resize(nr*ns*no, no);
          T2.fill(0);
          foreach (auto q, shells) {
            range Q = basis->range(q);
            int nq = Q.size();
            T1.resize(nq*nr*ns, no);
            T1.fill(0);
            foreach (auto p, shells) {
              range P = basis->range(p);
              int np = P.size();
              int2->compute_shell(s, r, q, p);
              auto G = Matrix::Map(int2->buffer(), np, nq*nr*ns);
              transform(P, C, G, T1);
            }
            T1.resize(nq, nr*ns*no);
            transform(Q, C, T1, T2);
          }
          T2.resize(nr, ns*no*no);
          transform(R, C, T2, T3);
        }
        T3.resize(ns, no*no*no);
        transform(S, C, T3, T4);
      }

      for (size_t l = 0, kl = 0; l < no; ++l) {
        for (size_t k = 0; k <= l; ++k, ++kl) {
          for (size_t j = 0, ij = 0; j < no; ++j) {
            for (size_t i = 0; i <= j; ++i, ++ij) {
              MPQC_ASSERT(
                  fabs(
                      V(index(j, i), index(k, l)) - V(index(i, j), index(k, l))
                          < 1e-14));
              auto v = T4(i + j * no + k * no * no, l);
              if (fabs(v) < 1e-15)
                v = 0;
              V(ij, kl) = v;
              //printf("integral %i %i %e\n", ij, kl, v);
            }
          }
        }
      }

    }

    void integrals(const Matrix &C, const Matrix &h_ao, Vector &h) {
      size_t no = C.rows();
      Matrix T1 = C * h_ao;
      Matrix T2 = T1 * C.transpose();
      h.resize((no * no + no) / 2);
      h.fill(0);
      for (size_t l = 0, kl = 0; l < no; ++l) {
        for (size_t k = 0; k <= l; ++k, ++kl) {
          h(kl) = T2(k, l);
        }
      }
    }

    void integrals(sc::Ref<sc::GaussianBasisSet> basis, const Matrix &C,
                   sc::Ref<sc::OneBodyInt> int1, Vector &h) {
      range shells(basis->nshell());
      size_t n = C.cols();
      Matrix h_ao(n, n);
      h_ao.fill(0);
      foreach (auto s, shells) {
        range S = basis->range(s);
        int ns = S.size();
        foreach (auto r, shells) {
          range R = basis->range(r);
          int nr = R.size();
          int1->compute_shell(s, r);
          h_ao(R,S) = Matrix::Map(int1->buffer(), nr, ns);
        }
      }
      integrals(C, h_ao, h);
    }

  } // namespace ci
} // namespace mpqc

#endif /* MPQC_CI_INTEGRALS_HPP */

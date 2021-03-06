

#ifndef MPQC4_SRC_MPQC_MATH_TENSOR_CLR_DECOMPOSED_TENSOR_GEMM_HELPER_H_
#define MPQC4_SRC_MPQC_MATH_TENSOR_CLR_DECOMPOSED_TENSOR_GEMM_HELPER_H_

#include "mpqc/math/tensor/clr/decomposed_tensor.h"
#include "mpqc/math/tensor/clr/decomposed_tensor_addition.h"
#include "mpqc/math/tensor/clr/decomposed_tensor_algebra.h"

namespace mpqc {
namespace math {
namespace detail {

extern bool recompress;

template <typename U>
using Dtensor = DecomposedTensor<U>;

using GHelper = TA::math::GemmHelper;

template <unsigned long... Dims>
struct low_rank_gemm {
  template <typename T>
  Dtensor<T> operator()(Dtensor<T> const &a, Dtensor<T> const &b,
                        const T factor, GHelper const &gh) {
    std::cout << "No overload for gemm with dimensions " << gh.result_rank()
              << " = " << gh.left_rank() << " * " << gh.right_rank()
              << std::endl;
    assert(false);
    return Dtensor<T>{};
  }
};

// Eri3("X,a,b") * C("b,k") = W("X, a, k")
template <>
struct low_rank_gemm<3ul, 3ul, 2ul> {
  template <typename T>
  Dtensor<T> operator()(Dtensor<T> const &a, Dtensor<T> const &b, const T f,
                        GHelper const &gh) {
    // for now assume that b is always full rank.
    if (a.ndecomp() == 1) {  // Reg gemm
      return Dtensor<T>{a.cut(), a.tensor(0).gemm(b.tensor(0), f, gh)};
    } else if (a.ndecomp() == 2) {  // LR gemm
      auto Rp = a.tensor(1).gemm(b.tensor(0), f, gh);
      auto tensor = Dtensor<T>{a.cut(), a.tensor(0).clone(), std::move(Rp)};

      if (recompress) {
        return tensor;
      } else {
        return Dtensor<T>{a.cut(), combine(tensor)};
      }
    }
    assert(false);
    return Dtensor<T>{};
  }

  template <typename T>
  Dtensor<T> &operator()(Dtensor<T> &c, Dtensor<T> const &a,
                         Dtensor<T> const &b, const T f, GHelper const &gh) {
    // assume b is never decomposed.
    if (c.ndecomp() == 1) {
      if (a.ndecomp() == 1) {
        c.tensor(0).gemm(a.tensor(0), b.tensor(0), f, gh);
      } else {
        auto Rp = a.tensor(1).gemm(b.tensor(0), f, gh);
        auto NoT = gh.left_op();
        auto gh =
            TA::math::GemmHelper(NoT, NoT, c.tensor(0).range().rank(),
                                 a.tensor(0).range().rank(), Rp.range().rank());
        c.tensor(0).gemm(a.tensor(0), Rp, f, gh);
      }
      return c;
    } else {
      if (a.ndecomp() == 1) {
        auto ab_tensor = a.tensor(0).gemm(b.tensor(0), 1.0, gh);
        const auto NoT = gh.left_op();
        auto gh = TA::math::GemmHelper(NoT, NoT, 3, 2, 3);
        ab_tensor.gemm(c.tensor(0), c.tensor(1), 1.0, gh);
        c = DecomposedTensor<double>(c.cut(), std::move(ab_tensor));

        return c;
      }

      c = add(c, this->operator()(a, b, f, gh));
      auto const &c_left_extent = c.tensor(0).range().extent();
      auto const &c_right_extent = c.tensor(1).range().extent();
      const auto long_dim = c_right_extent[1] * c_right_extent[2];
      auto out_dim = c.rank();
      const auto full_rank = std::min(c_left_extent[0], long_dim);

      if (recompress) {
        recompress(c);
        out_dim = c.rank();
      }

      if (out_dim > full_rank / 2) {
        c = DecomposedTensor<T>(c.cut(), combine(c));
      }

      return c;
    }

    assert(false);
  }
};

// B("X,a,b") = V^{-1}("X,P") * X("P,a,b")
template <>
struct low_rank_gemm<3ul, 2ul, 3ul> {
  template <typename T>
  Dtensor<T> operator()(Dtensor<T> const &a, Dtensor<T> const &b, const T f,
                        GHelper const &gh) {
    Dtensor<T> out;
    if (a.ndecomp() == 1) {    // Full *
      if (b.ndecomp() == 1) {  // Full Full
        out = Dtensor<T>{a.cut(), a.tensor(0).gemm(b.tensor(0), f, gh)};
      } else {  // Full Low
        const auto gh_fl = GHelper(NoT, NoT, 2, 2, 2);
        auto sw = a.tensor(0).gemm(b.tensor(0), f, gh_fl);
        out = Dtensor<T>{a.cut(), std::move(sw), b.tensor(1).clone()};
      }
    } else {                   // Low *
      if (b.ndecomp() == 1) {  // Low Full
        auto Tw = a.tensor(1).gemm(b.tensor(0), f, gh);
        out = Dtensor<T>{a.cut(), a.tensor(0).clone(), std::move(Tw)};
      } else {  // Low Low
        const auto gh_m = GHelper(NoT, NoT, 2, 2, 2);
        auto M = a.tensor(1).gemm(b.tensor(0), 1.0, gh_m);

        // Want to contract over largest dim so output is small.
        const auto go_left = a.rank() >= b.rank();
        if (go_left) {
          const auto gh_left = GHelper(NoT, NoT, 2, 2, 2);
          out = Dtensor<T>{a.cut(), a.tensor(0).gemm(M, f, gh_left),
                           b.tensor(1).clone()};
        } else {
          const auto gh_right = GHelper(NoT, NoT, 3, 2, 3);
          out = Dtensor<T>{a.cut(), a.tensor(0).clone(),
                           M.gemm(b.tensor(1), f, gh_right)};
        }
      }
    }

    if (recompress && out.cut() != 0.0) {
      if (out.ndecomp() == 1) {
        auto test = two_way_decomposition(out);
        if (!test.empty()) {
          out = test;
        }
      } else {
        recompress(out);
      }
    }

    return out;
  }

  template <typename T>
  Dtensor<T> &operator()(Dtensor<T> &c, Dtensor<T> const &a,
                         Dtensor<T> const &b, const T f, GHelper const &gh) {
    if (c.ndecomp() == 1) {      // Full * *
      if (a.ndecomp() == 1) {    // Full Full *
        if (b.ndecomp() == 1) {  //  Full Full Full
          c.tensor(0).gemm(a.tensor(0), b.tensor(0), f, gh);
        } else {  // Full Full Low
          const auto gh_left = GHelper(NoT, NoT, 2, 2, 2);
          c.tensor(0).gemm(a.tensor(0).gemm(b.tensor(0), 1.0, gh_left),
                           b.tensor(1), f, gh);
        }
      } else {                   // Full Low *
        if (b.ndecomp() == 1) {  //  Full Low Full
          const auto gh_right = GHelper(NoT, NoT, 3, 2, 3);
          c.tensor(0).gemm(a.tensor(0),
                           a.tensor(1).gemm(b.tensor(0), 1.0, gh_right), f, gh);
        } else {  //  Full Low Low
          const auto gh_m = GHelper(NoT, NoT, 2, 2, 2);
          auto mid = a.tensor(1).gemm(b.tensor(0), 1.0, gh_m);

          auto const &sa_sizes = a.tensor(0).range().extent();
          auto const &tb_sizes = b.tensor(1).range().extent();
          auto const &m_sizes = mid.range().extent();
          const auto r1 = m_sizes[0];
          const auto r2 = m_sizes[1];
          const auto X = sa_sizes[0];
          const auto ab = tb_sizes[1] * tb_sizes[2];

          auto left_cost = X * r2 * (r1 + ab);
          auto right_cost = ab * r1 * (r2 + X);

          if (left_cost <= right_cost) {
            auto gh_left = GHelper(NoT, NoT, 2, 2, 2);
            auto left_tensor = a.tensor(0).gemm(mid, 1.0, gh_left);
            c.tensor(0).gemm(left_tensor, b.tensor(1), f, gh);
          } else {
            auto gh_right = GHelper(NoT, NoT, 3, 2, 3);
            auto right_tensor = mid.gemm(b.tensor(1), 1.0, gh_right);
            c.tensor(0).gemm(a.tensor(0), right_tensor, f, gh);
          }
        }
      }
    } else {                                       // Low * *
      if (a.ndecomp() == 1 && b.ndecomp() == 1) {  // Low Full Full
        auto temp = a.tensor(0).gemm(b.tensor(0), f, gh);
        temp.gemm(c.tensor(0), c.tensor(1), 1.0, gh);

        c = Dtensor<T>{c.cut(), std::move(temp)};

      } else {  // Other cases can be handled by the following LFL LLF and
                // LLL
        auto ab = this->operator()(a, b, f, gh);
        c = add(c, ab);
      }
    }

    if (recompress && c.cut() != 0.0) {
      if (c.ndecomp() == 1) {
        auto test = two_way_decomposition(c);
        if (!test.empty()) {
          c = test;
        }
      } else {
        recompress(c);

        auto out_dim = c.rank();
        auto left_dim = c.tensor(0).range().extent()[0];
        auto const &right_extent = c.tensor(1).range().extent();
        auto right_dim = right_extent[1] * right_extent[2];
        if (out_dim >= std::min(left_dim, right_dim) / 2) {
          c = Dtensor<T>(c.cut(), combine(c));
        }
      }
    }

    return c;
  }
};

// K("i,j") = W("X,k,i") * W("X,k,j")
template <>
struct low_rank_gemm<2ul, 3ul, 3ul> {
  // just use the 3 way functions
  template <typename T>
  Dtensor<T> operator()(Dtensor<T> const &a, Dtensor<T> const &b, const T f,
                        GHelper const &gh) {
    const auto left_dim = a.tensors().back().range().extent()[2];
    const auto right_dim = b.tensors().back().range().extent()[2];
    const auto out_range = TA::Range{left_dim, right_dim};
    const auto out_tensor = TA::Tensor<T>(std::move(out_range), T(0.0));
    auto out = Dtensor<T>(a.cut(), std::move(out_tensor));
    this->operator()(out, a, b, f, gh);
    return out;
  }

  // B^{T}_{a,Xk} * W_{Xk,b}
  template <typename T>
  Dtensor<T> &operator()(Dtensor<T> &c, Dtensor<T> const &a,
                         Dtensor<T> const &b, const T f, GHelper const &gh) {
    assert(1 == c.ndecomp());
    if (a.ndecomp() == 1) {    // Full *
      if (b.ndecomp() == 1) {  //  Full Full
        c.tensor(0).gemm(a.tensor(0), b.tensor(0), f, gh);
      } else {  // Full Low we are going to do this backwards.
        // S^{wT}_{r,X} * B_{X,ka} = M_{r,ka}
        const auto gh_right = GHelper(Tr, NoT, 3, 2, 3);
        auto M = b.tensor(0).gemm(a.tensor(0), 1.0, gh_right);

        // M^{T}_{a, rk} * T^{w}_{rk, b} = K_{a,b}
        const auto gh_last = GHelper(Tr, NoT, 2, 3, 3);
        c.tensor(0).gemm(M, b.tensor(1), f, gh_last);
      }
    } else {                   //  Low *
      if (b.ndecomp() == 1) {  //  Low Full
        // S^{bT}_{r,X} * W_{X,kb} = M_{r, kb}
        const auto gh_right = GHelper(Tr, NoT, 3, 2, 3);
        auto M = a.tensor(0).gemm(b.tensor(0), 1.0, gh_right);

        // T^{bT}_{a, rk} * M_{rk, b} = K_{a,b}
        const auto gh_last = GHelper(Tr, NoT, 2, 3, 3);
        c.tensor(0).gemm(a.tensor(1), M, f, gh_last);
      } else {  //  Low Low
        // S^{bT}_{r1,X} * S^{w}_{X,r2} = M_{r1,r2}
        const auto gh_mid = GHelper(Tr, NoT, 2, 2, 2);
        auto M = a.tensor(0).gemm(b.tensor(0), 1.0, gh_mid);

        // M_{r1,r2} * T^{w}_{r2,kb} = Mp_{r1,kb}
        const auto gh_right = GHelper(NoT, NoT, 3, 2, 3);
        auto Mp = M.gemm(b.tensor(1), 1.0, gh_right);

        // T^{bT}_{a,r1k} * Mp_{r1k, b} = K_{a,b}
        const auto gh_last = GHelper(Tr, NoT, 2, 3, 3);
        c.tensor(0).gemm(a.tensor(1), Mp, f, gh_last);
      }
    }
    return c;
  }
};

// Eri3("X,a,b") * D("a,b") = M("X")
template <>
struct low_rank_gemm<1ul, 3ul, 2ul> {
  template <typename T>
  Dtensor<T> operator()(Dtensor<T> const &a, Dtensor<T> const &b, const T f,
                        GHelper const &gh) {
    const auto X = a.tensor(0).range().extent()[0];  // X from above.
    auto range = TA::Range(X);
    auto out_tensor = Dtensor<T>(a.cut(), TA::Tensor<T>(std::move(range), 0.0));
    this->operator()(out_tensor, a, b, f, gh);
    return out_tensor;
  }

  template <typename T>
  Dtensor<T> &operator()(Dtensor<T> &c, Dtensor<T> const &a,
                         Dtensor<T> const &b, const T f, GHelper const &gh) {
    // assume b and c are never decomposed.
    if (c.ndecomp() == 1) {
      if (a.ndecomp() == 1) {
        c.tensor(0).gemm(a.tensor(0), b.tensor(0), f, gh);
      } else {
        auto gh_right = TA::math::GemmHelper(NoT, NoT, 1, 3, 2);
        auto Rp = a.tensor(1).gemm(b.tensor(0), 1.0, gh_right);
        auto gh_left = TA::math::GemmHelper(NoT, NoT, 1, 2, 1);
        c.tensor(0).gemm(a.tensor(0), Rp, f, gh_left);
      }
      return c;
    }
    assert(false);
  }
};

// ("X") = M("X,P") * ("P")
template <>
struct low_rank_gemm<1ul, 2ul, 1ul> {
  template <typename T>
  Dtensor<T> operator()(Dtensor<T> const &a, Dtensor<T> const &b, const T f,
                        GHelper const &gh) {
    const auto X = a.tensor(0).range().extent()[0];  // X from above.
    auto range = TA::Range(X);
    auto out_tensor = Dtensor<T>(a.cut(), TA::Tensor<T>(std::move(range), 0.0));
    this->operator()(out_tensor, a, b, f, gh);
    return out_tensor;
  }

  template <typename T>
  Dtensor<T> &operator()(Dtensor<T> &c, Dtensor<T> const &a,
                         Dtensor<T> const &b, const T f, GHelper const &gh) {
    // assume b and c are never decomposed.
    if (c.ndecomp() == 1) {
      if (a.ndecomp() == 1) {
        c.tensor(0).gemm(a.tensor(0), b.tensor(0), f, gh);
      } else {
        auto gh_right = TA::math::GemmHelper(NoT, NoT, 1, 2, 1);
        auto Rp = a.tensor(1).gemm(b.tensor(0), 1.0, gh_right);
        auto gh_left = TA::math::GemmHelper(NoT, NoT, 1, 2, 1);
        c.tensor(0).gemm(a.tensor(0), Rp, f, gh_left);
      }
      return c;
    }
    assert(false);
  }
};

// W("X,i,j") * M("X") = J("i,j");
template <>
struct low_rank_gemm<2ul, 3ul, 1ul> {
  template <typename T>
  Dtensor<T> operator()(Dtensor<T> const &a, Dtensor<T> const &b, const T f,
                        GHelper const &gh) {
    auto out_tensor = Dtensor<T>(a.cut());
    if (a.ndecomp() == 1) {
      const auto i = a.tensor(0).range().extent()[1];
      const auto j = a.tensor(0).range().extent()[2];
      auto range = TA::Range{i, j};
      out_tensor = Dtensor<T>(a.cut(), TA::Tensor<T>(std::move(range), 0.0));
    } else if (a.ndecomp() == 2) {
      const auto i = a.tensor(1).range().extent()[1];
      const auto j = a.tensor(1).range().extent()[2];
      auto range = TA::Range(i, j);
      out_tensor = Dtensor<T>(a.cut(), TA::Tensor<T>(std::move(range), 0.0));
    }
    this->operator()(out_tensor, a, b, f, gh);
    return out_tensor;
  }

  template <typename T>
  Dtensor<T> &operator()(Dtensor<T> &c, Dtensor<T> const &a,
                         Dtensor<T> const &b, const T f, GHelper const &gh) {
    // assume b and c are never decomposed.
    if (c.ndecomp() == 1) {
      if (a.ndecomp() == 1) {
        c.tensor(0).gemm(a.tensor(0), b.tensor(0), f, gh);
      } else {
        auto gh_right = TA::math::GemmHelper(Tr, NoT, 1, 2, 1);
        auto Rp = a.tensor(0).gemm(b.tensor(0), 1.0, gh_right);
        auto gh_left = TA::math::GemmHelper(Tr, NoT, 2, 3, 1);
        c.tensor(0).gemm(a.tensor(1), Rp, f, gh_left);
      }
      return c;
    }
    assert(false);
  }
};

// C("i,j") = A("i,k") * B("k,j") + C("i,j");
template <>
struct low_rank_gemm<2ul, 2ul, 2ul> {
  template <typename T>
  Dtensor<T> operator()(Dtensor<T> const &a, Dtensor<T> const &b, const T f,
                        GHelper const &gh) {
    if (a.ndecomp() == 1) {
      if (b.ndecomp() == 1) {
        return Dtensor<T>(a.cut(), a.tensor(0).gemm(b.tensor(0), f, gh));
      } else {
        auto left = a.tensor(0).gemm(b.tensor(0), f, gh);
        auto right = b.tensor(1).clone();
        return Dtensor<T>(a.cut(), std::move(left), std::move(right));
      }
    } else {
      if (b.ndecomp() == 1) {
        auto left = a.tensor(0).clone();
        auto right = a.tensor(1).gemm(b.tensor(0), f, gh);
        return Dtensor<T>(a.cut(), std::move(left), std::move(right));
      } else {
        auto mid = a.tensor(1).gemm(b.tensor(0), 1, gh);
        decltype(mid) left, right;
        if (a.rank() > b.rank()) {
          left = a.tensor(0).gemm(mid, f, gh);
          right = b.tensor(1).clone();
        } else {
          left = a.tensor(0).clone();
          right = mid.gemm(b.tensor(1), f, gh);
        }
        return Dtensor<T>(a.cut(), std::move(left), std::move(right));
      }
    }
    assert(false);
    return Dtensor<T>(a.cut());
  }

  template <typename T>
  Dtensor<T> &operator()(Dtensor<T> &c, Dtensor<T> const &a,
                         Dtensor<T> const &b, const T f, GHelper const &gh) {
    /* if (c.ndecomp() == 1) { */
    /*     if (a.ndecomp() == 1) { */
    /*         if(b.ndecomp() == 1){ */
    /*             c.tensor(0).gemm(a.tensor(0), b.tensor(0), f, gh); */
    /*             return c; */
    /*         } */
    /*     } */
    /* } */
    auto c_contrib = this->operator()(a, b, f, gh);
    add_to(c, c_contrib);
    return c;
  }
};

}  // namespace detail
}  // namespace math
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_MATH_TENSOR_CLR_DECOMPOSED_TENSOR_GEMM_HELPER_H_

#pragma once

#include "tensor_fwd.h"
#include "decomp_helper.h"
#include <vector>
#include <tuple>


namespace tcc {
namespace tensor {

class DecompTensor {
  public:
    using value_type = double;
    using numeric_type = double;

  public:
    // Compiler provided constructors
    DecompTensor() = default;
    DecompTensor(DecompTensor const &) = default;
    DecompTensor(DecompTensor &&) = default;
    DecompTensor &operator=(DecompTensor const &) = default;
    DecompTensor &operator=(DecompTensor &&) = default;

    DecompTensor(TARange const &r, double c) : range_(r), cut_(c) {}
    DecompTensor(TARange const &r, TATensor t, double c)
        : range_(r), cut_(c), tensors_{std::move(t)} {}

    DecompTensor(TARange const &r, TATensor t1, TATensor t2, double c)
        : range_(r), cut_(c), tensors_{std::move(t1), std::move(t2)} {}

    DecompTensor(TARange const &r, TATensor t1, TATensor t2, TATensor t3,
                 double c)
        : range_(r),
          cut_(c),
          tensors_{std::move(t1), std::move(t2), std::move(t3)} {}

    DecompTensor(TARange const &r, TATensor t1, TATensor t2, TATensor t3,
                 TATensor t4, double c)
        : range_(r),
          cut_(c),
          tensors_{std::move(t1), std::move(t2), std::move(t3), std::move(t4)} {
    }

    DecompTensor(TARange const &r, std::vector<TATensor> tvec, double c)
        : range_(r), cut_(c), tensors_{std::move(tvec)} {}

    DecompTensor(std::tuple<TARange, std::vector<TATensor>, double> &&tup)
        : range_(std::move(std::get<0>(tup))),
          cut_(std::move(std::get<2>(tup))),
          tensors_{std::move(std::get<1>(tup))} {}

    DecompTensor(DecompTensor const &other, Permutation const &perm) {}

    bool is_decomposable() const { return tensors_.size() >= 2; }
    const TARange &range() const { return range_; }
    unsigned long size() const { return range_.volume(); }
    bool empty() const { return tensors_.empty(); }
    unsigned long ndivisions() const { return tensors_.size(); }
    TATensor const &tensor(unsigned long n) const { return tensors_[n]; }
    std::vector<TATensor> const &tensors() const { return tensors_; }
    double cut() const { return cut_; }

    TATensor full_tensor() const { return detail::to_full(range(), tensors()); }

    template <typename Archive>
    typename madness::enable_if<madness::archive::is_output_archive<Archive>>::
        type
        serialize(Archive &ar) {}

    /// Serialize tensor data
    template <typename Archive>
    typename madness::enable_if<madness::archive::is_input_archive<Archive>>::
        type
        serialize(Archive &ar) {}


    void swap(DecompTensor &other) {
        using std::swap;
        swap(tensors_, other.tensors_);
        swap(range_, other.range_);
    }

    DecompTensor permute(Permutation const &perm) const {
        return DecompTensor(*this, perm);
    }

    template <typename Op>
    DecompTensor binary(DecompTensor const &other, Op const &op) const {
        return DecompTensor(*this, other, op);
    }

    template <typename Op>
    DecompTensor
    binary(DecompTensor &other, const Op &op, const Permutation &perm) const {
        return DecompTensor(*this, other, op, perm);
    }

    template <typename Op>
    DecompTensor unary(const Op &op) const {
        return DecompTensor(*this, op);
    }

    template <typename Op>
    DecompTensor unary(const Op &op, const Permutation &perm) const {
        return DecompTensor(*this, op, perm);
    }

    template <typename Op>
    DecompTensor &inplace_unary(const Op &op) {
        assert(false);
    }

    DecompTensor scale(double factor) const { assert(false); }

    DecompTensor scale(double factor, const Permutation &perm) const {
        assert(false);
    }

    DecompTensor &scale_to(double factor) { assert(false); }

    DecompTensor add(const DecompTensor &other) const {
        return this->add(other, 1.0);
    }

    DecompTensor add(const DecompTensor &other, const Permutation &perm) const {
        assert(false);
    }

    DecompTensor add(const DecompTensor &other, double factor) const {
        if (ndivisions() == 1) {
            if (other.ndivisions() == 1) { // Full Full
                return DecompTensor{
                    detail::add<1, 1>(range(), tensors(), cut(), other.range(),
                                      other.tensors(), other.cut(), factor)};
            } else if (other.ndivisions() == 2) { // Full Low
                return DecompTensor{
                    detail::add<1, 2>(range(), tensors(), cut(), other.range(),
                                      other.tensors(), other.cut(), factor)};
            } else {
                assert(false);
            }
        } else if (ndivisions() == 2) {
            if (other.ndivisions() == 1) { // Low Full
                return DecompTensor{
                    detail::add<2, 1>(range(), tensors(), cut(), other.range(),
                                      other.tensors(), other.cut(), factor)};
            } else if (other.ndivisions() == 2) { // Low Low
                return DecompTensor{
                    detail::add<2, 2>(range(), tensors(), cut(), other.range(),
                                      other.tensors(), other.cut(), factor)};
            } else {
                assert(false);
            }
        } else {
            assert(false);
        }
    }

    DecompTensor add(const DecompTensor &other, double factor,
                     const Permutation &perm) const {
        assert(false);
    }

    DecompTensor add(double value) const {}

    DecompTensor add(double value, const Permutation &perm) const {
        assert(false);
    }

    DecompTensor &add_to(const DecompTensor &other) {
        *this = this->add(other, 1.0);
        return *this;
    }

    template <typename U, typename AU>
    DecompTensor &add_to(const DecompTensor &other, double factor) {
        *this = this->add(other,factor);
        return *this;
    }

    DecompTensor &add_to(double value) { assert(false); }

    template <typename U, typename AU>
    DecompTensor subt(const DecompTensor &other) const {
        assert(false);
    }

    template <typename U, typename AU>
    DecompTensor
    subt(const DecompTensor &other, const Permutation &perm) const {
        assert(false);
    }

    DecompTensor subt(const DecompTensor &other, double factor) const {
        assert(false);
    }

    DecompTensor subt(const DecompTensor &other, double factor,
                      const Permutation &perm) const {
        assert(false);
    }

    DecompTensor subt(double value) const { return add(-value); }

    DecompTensor subt(double value, const Permutation &perm) const {
        return add(-value, perm);
    }

    DecompTensor &subt_to(const DecompTensor &other) { assert(false); }

    template <typename U, typename AU>
    DecompTensor &subt_to(const DecompTensor &other, double factor) {
        assert(false);
    }

    DecompTensor &subt_to(double value) { return add_to(-value); }

    DecompTensor mult(const DecompTensor &other) const { assert(false); }

    DecompTensor
    mult(const DecompTensor &other, const Permutation &perm) const {
        assert(false);
    }

    DecompTensor mult(const DecompTensor &other, double factor) const {
        assert(false);
    }

    template <typename U, typename AU>
    DecompTensor mult(const DecompTensor &other, double factor,
                      const Permutation &perm) const {
        assert(false);
    }

    DecompTensor &mult_to(const DecompTensor &other) { assert(false); }

    DecompTensor &mult_to(const DecompTensor &other, double factor) {
        assert(false);
    }

    DecompTensor neg() const { assert(false); }

    DecompTensor neg(const Permutation &perm) const { assert(false); }

    DecompTensor &neg_to() { assert(false); }

    DecompTensor gemm(const DecompTensor &other, double factor,
                      const math::GemmHelper &gemm_helper) const {
        auto out_range
            = gemm_helper.make_result_range<TARange>(range(), other.range());
        if (ndivisions() == 1) {
            if (other.ndivisions() == 1) { // Full Full
                return DecompTensor(detail::gemm<1, 1>(
                    out_range, range(), tensors(), cut(), other.range(),
                    other.tensors(), other.cut(), factor, gemm_helper));
            } else if (other.ndivisions() == 2) { // Full Low
                return DecompTensor(detail::gemm<1, 2>(
                    out_range, range(), tensors(), cut(), other.range(),
                    other.tensors(), other.cut(), factor, gemm_helper));
            }
        } else if (ndivisions() == 2) {           // Low Low
            if (other.ndivisions() == 1) {        // Full Low
                return DecompTensor(detail::gemm<2, 1>(
                    out_range, range(), tensors(), cut(), other.range(),
                    other.tensors(), other.cut(), factor, gemm_helper));
            } else if (other.ndivisions() == 2) { // Low Low
                return DecompTensor(detail::gemm<2, 2>(
                    out_range, range(), tensors(), cut(), other.range(),
                    other.tensors(), other.cut(), factor, gemm_helper));
            }
        }
        assert(false);
    }

    DecompTensor &gemm(const DecompTensor &left, const DecompTensor &right,
                       double factor, const math::GemmHelper &gemm_helper) {
        return add_to(left.gemm(right, factor, gemm_helper));
    }

    double trace() const { assert(false); }

  public:
    template <typename Op>
    double reduce(double init_value, const Op &op) const {
        assert(false);
    }

    template <typename Op>
    double
    reduce(const DecompTensor &other, double init_value, const Op &op) const {
        assert(false);
    }

    double sum() const { assert(false); }

    double product() const { assert(false); }

    double squared_norm() const { assert(false); }

    double norm() const { assert(false); }

    double min() const { assert(false); }

    double max() const { assert(false); }

    double abs_min() const { assert(false); }

    double abs_max() const { assert(false); }

    double dot(const DecompTensor &other) const { assert(false); }


  private:
    TARange range_;
    double cut_;
    std::vector<TATensor> tensors_;
};

inline DecompTensor
operator+(const DecompTensor &left, const DecompTensor &right) {
    assert(false);
}

inline DecompTensor
operator-(const DecompTensor &left, const DecompTensor &right) {
    assert(false);
}

inline DecompTensor
operator*(const DecompTensor &left, const DecompTensor &right) {
    assert(false);
}

inline DecompTensor operator^(const Permutation &perm,
                              const DecompTensor &tensor) { assert(false); }

inline std::ostream &
operator<<(std::ostream &os, const DecompTensor &t) {
    assert(false);
}

} // namespace tensor
} // namespace tcc

#pragma once

#include "tensor_fwd.h"
#include "decomp_helper.h"
#include <vector>
#include <algorithm>
#include <numeric>
#include <tuple>


namespace tcc {
namespace tensor {

/**
 * Class which specializes in holding an Exchange Tensor.
 **/
class CoulombTensor {
  private:
    TARange range_;
    double cut_;
    std::vector<TATensor> tensors_;

  public:
    using value_type = double;
    using numeric_type = double;

  public:
    // Compiler provided constructors
    CoulombTensor() = default;
    CoulombTensor(CoulombTensor const &) = default;
    CoulombTensor(CoulombTensor &&) = default;
    CoulombTensor &operator=(CoulombTensor const &) = default;
    CoulombTensor &operator=(CoulombTensor &&) = default;

    CoulombTensor(TARange const &r, double c) : range_(r), cut_(c) {}
    CoulombTensor(TARange const &r, TATensor t, double c)
        : range_(r), cut_(c), tensors_{std::move(t)} {}

    CoulombTensor(TARange const &r, TATensor t1, TATensor t2, TATensor t3,
                   double c)
        : range_(r),
          cut_(c),
          tensors_{std::move(t1), std::move(t2), std::move(t3)} {}

    CoulombTensor(TARange const &r, std::vector<TATensor> tvec, double c)
        : range_(r), cut_(c), tensors_{std::move(tvec)} {}

    CoulombTensor(std::tuple<TARange, std::vector<TATensor>> const &tup,
                   double c)
        : range_(std::get<0>(tup)), cut_(c), tensors_(std::get<1>(tup)) {}

    CoulombTensor(CoulombTensor const &other, Permutation const &perm) {}

    // Normal Tiledarary tensor functions
    const TARange &range() const { return range_; }
    unsigned long size() const { return range_.volume(); }
    bool empty() const { return tensors_.empty(); }

    // Functions for the low rank tensor.
    bool is_decomposed() const { return tensors_.size() >= 2; }

    std::vector<unsigned int> ranks() const {
        assert(!tensors_.empty());
        std::vector<unsigned int> ranks;

        if (tensors_.size() == 1) {
            ranks.push_back(range_.size()[0] * range_.size()[1]);
        } else {
            assert(tensors_.size() == 2);
            ranks.push_back(tensors_[0].range().size()[1]);
        }

        return ranks;
    }


    unsigned long compressed_size() const {
        return std::accumulate(std::begin(tensors_), std::end(tensors_), 0,
                               [](unsigned long int size, TATensor const &t) {
            return size + t.range().volume();
        });
    }

    double cut() const { return cut_; }

    template <typename Archive>
    typename mad::enable_if<mad::archive::is_output_archive<Archive>>::type
    serialize(Archive &ar) {}

    /// Serialize tensor data
    template <typename Archive>
    typename mad::enable_if<mad::archive::is_input_archive<Archive>>::type
    serialize(Archive &ar) {}

    void swap(CoulombTensor &other) {
        using std::swap;
        swap(tensors_, other.tensors_);
        swap(range_, other.range_);
    }

    CoulombTensor permute(Permutation const &perm) const {
        return CoulombTensor(*this, perm);
    }

    template <typename Op>
    CoulombTensor binary(CoulombTensor const &other, Op const &op) const {
        return CoulombTensor(*this, other, op);
    }

    template <typename Op>
    CoulombTensor
    binary(CoulombTensor &other, const Op &op, const Permutation &perm) const {
        return CoulombTensor(*this, other, op, perm);
    }

    template <typename Op>
    CoulombTensor unary(const Op &op) const {
        return CoulombTensor(*this, op);
    }

    template <typename Op>
    CoulombTensor unary(const Op &op, const Permutation &perm) const {
        return CoulombTensor(*this, op, perm);
    }

    template <typename Op>
    CoulombTensor &inplace_unary(const Op &op) {
        assert(false);
    }

    CoulombTensor scale(double factor) const { assert(false); }

    CoulombTensor scale(double factor, const Permutation &perm) const {
        assert(false);
    }

    CoulombTensor &scale_to(double factor) { assert(false); }

    CoulombTensor add(const CoulombTensor &other) const {
        return this->add(other, 1.0);
    }

    CoulombTensor
    add(const CoulombTensor &other, const Permutation &perm) const {
        assert(false);
    }

    CoulombTensor add(const CoulombTensor &other, double factor) const {
        assert(false);
    }

    CoulombTensor add(const CoulombTensor &other, double factor,
                       const Permutation &perm) const {
        assert(false);
    }

    CoulombTensor add(double value) const { assert(false); }

    CoulombTensor add(double value, const Permutation &perm) const {
        assert(false);
    }

    CoulombTensor &add_to(const CoulombTensor &other) { assert(false); }

    template <typename U, typename AU>
    CoulombTensor &add_to(const CoulombTensor &other, double factor) {
        assert(false);
    }

    CoulombTensor &add_to(double value) { assert(false); }

    template <typename U, typename AU>
    CoulombTensor subt(const CoulombTensor &other) const {
        assert(false);
    }

    template <typename U, typename AU>
    CoulombTensor
    subt(const CoulombTensor &other, const Permutation &perm) const {
        assert(false);
    }

    CoulombTensor subt(const CoulombTensor &other, double factor) const {
        assert(false);
    }

    CoulombTensor subt(const CoulombTensor &other, double factor,
                        const Permutation &perm) const {
        assert(false);
    }

    CoulombTensor subt(double value) const { assert(false); }

    CoulombTensor subt(double value, const Permutation &perm) const {
        assert(false);
    }

    CoulombTensor &subt_to(const CoulombTensor &other) { assert(false); }

    template <typename U, typename AU>
    CoulombTensor &subt_to(const CoulombTensor &other, double factor) {
        assert(false);
    }

    CoulombTensor &subt_to(double value) { assert(false); }

    CoulombTensor mult(const CoulombTensor &other) const { assert(false); }

    CoulombTensor
    mult(const CoulombTensor &other, const Permutation &perm) const {
        assert(false);
    }

    CoulombTensor mult(const CoulombTensor &other, double factor) const {
        assert(false);
    }

    template <typename U, typename AU>
    CoulombTensor mult(const CoulombTensor &other, double factor,
                        const Permutation &perm) const {
        assert(false);
    }

    CoulombTensor &mult_to(const CoulombTensor &other) { assert(false); }

    CoulombTensor &mult_to(const CoulombTensor &other, double factor) {
        assert(false);
    }

    CoulombTensor neg() const { assert(false); }

    CoulombTensor neg(const Permutation &perm) const { assert(false); }

    CoulombTensor &neg_to() { assert(false); }

    CoulombTensor gemm(const CoulombTensor &other, double factor,
                        const math::GemmHelper &gemm_helper) const {
        assert(false);
    }

    CoulombTensor &gemm(const CoulombTensor &left,
                         const CoulombTensor &right, double factor,
                         const math::GemmHelper &gemm_helper) {
        assert(false);
    }

    double trace() const { assert(false); }

  public:
    template <typename Op>
    double reduce(double init_value, const Op &op) const {
        assert(false);
    }

    template <typename Op>
    double
    reduce(const CoulombTensor &other, double init_value, const Op &op) const {
        assert(false);
    }

    double sum() const { assert(false); }

    double product() const { assert(false); }

    double squared_norm() const { assert(false); }

    double norm() const {
        if (tensors_.empty()) {
            return 0.0;
        } else {
            auto sum = 0.0;
            for (auto const &t : tensors_) {
                if(!t.empty()){
                    sum += t.norm();
                }
            }
            return sum;
        }
    }

    double min() const { assert(false); }

    double max() const { assert(false); }

    double abs_min() const { assert(false); }

    double abs_max() const { assert(false); }

    double dot(const CoulombTensor &other) const { assert(false); }
};

inline CoulombTensor
operator+(const CoulombTensor &left, const CoulombTensor &right) {
    assert(false);
}

inline CoulombTensor
operator-(const CoulombTensor &left, const CoulombTensor &right) {
    assert(false);
}

inline CoulombTensor
operator*(const CoulombTensor &left, const CoulombTensor &right) {
    assert(false);
}

inline CoulombTensor operator^(const Permutation &perm,
                                const CoulombTensor &tensor) { assert(false); }

inline std::ostream &
operator<<(std::ostream &os, const CoulombTensor &t) {
    assert(false);
}

} // namespace tensor
} // namespace tcc

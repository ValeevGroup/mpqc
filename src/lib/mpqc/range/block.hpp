#ifndef MPQC_RANGE_BLOCK_HPP
#define MPQC_RANGE_BLOCK_HPP

#include "mpqc/range.hpp"

namespace mpqc {

    struct range::block_list {

        block_list(const std::vector<range> &dims,
                         const std::vector<size_t> &block) {
            initialize(dims, block);
        }

        block_list(const std::vector<range> &dims, size_t N) {
            MPQC_ASSERT(N > 0);
            std::vector<size_t> block(dims.size(), 1);
            for (size_t i = 0; i < dims.size(); ++i) {
                block[i] = std::min<size_t>(dims[i].size(), N);
                N /= dims[i].size();
                if (!N) break;
            }
            initialize(dims, block);
        }

        size_t size() const {
            return size_;
        }

        std::vector<range> operator[](size_t index) const {
            MPQC_ASSERT(index < size_);
            size_t N = dims_.size();
            std::vector<size_t> t(N);
            for (size_t i = N; i > 0; ) {
                --i;
                t[i] = index/strides_[i];
                index = index%strides_[i];
            }
            std::vector<range> r;
            for (size_t i = 0; i < N; ++i) {
                int begin = t[i]*block_[i];
                int end = std::min<int>(begin + block_[i], *dims_[i].end());
                r.push_back(range(begin, end));
            }
            return r;
        }

    private:
        std::vector<range> dims_;
        std::vector<size_t> block_;
        std::vector<size_t> strides_;
        size_t size_;

        void initialize(const std::vector<range> &dims,
                   const std::vector<size_t> &block) {
            dims_ = dims;
            block_ = block;
            size_t N = dims_.size();
            size_ = 1;
            for (size_t i = 0; i < N; ++i) {
                strides_.push_back(size_);
                size_ *= (dims_[i].size() + block_[i] + 1)/block_[i];
            }
        }
    };

} // namespace mpqc

#endif // MPQC_RANGE_BLOCK_HPP

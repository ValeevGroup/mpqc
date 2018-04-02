//
// Created by Chong Peng on 10/6/17.
//

#ifndef SRC_MPQC_MATH_EXTERNAL_TILEDARRAY_REDUCTION_H_
#define SRC_MPQC_MATH_EXTERNAL_TILEDARRAY_REDUCTION_H_

#include <algorithm>
#include <functional>
#include <numeric>

namespace mpqc {
namespace detail {

template <typename Element>
struct AbsMaxCompare {
  bool operator()(const Element& i, const Element& j) {
    return std::abs(i) > std::abs(j);
  }
};

template <typename Element, typename Index>
struct MaxPairCompare {
  bool operator()(const std::pair<Element, Index>& i,
                  const std::pair<Element, Index>& j) {
    return i.first > j.first;
  }
};

template <typename Element, typename Index>
struct AbsMaxPairCompare {
  bool operator()(const std::pair<Element, Index>& i,
                  const std::pair<Element, Index>& j) {
    return std::abs(i.first) > std::abs(j.first);
  }
};

template <typename Element, typename Index>
struct PairApproxUnique {
  bool operator()(const std::pair<Element, Index>& i,
                  const std::pair<Element, Index>& j) {
    return std::abs(i.first - j.first) < 1.0e-14;
  }
};

template <typename Tile>
class MaxNReduction {
 public:
  // typedefs
  using element_type = typename Tile::numeric_type;
  using result_type = std::vector<element_type>;
  using argument_type = Tile;

  // constructor
  MaxNReduction(std::size_t n)
      : n_(n), default_(std::numeric_limits<element_type>::min()) {}
  ~MaxNReduction() = default;

  // Reduction functions

  // Make an empty result object
  result_type operator()() const { return result_type(n_, default_); }

  // Post process the result
  const result_type& operator()(const result_type& result) const {
    return result;
  }

  // Reduce two result objects
  void operator()(result_type& result, const result_type& arg) const {
    result.insert(result.end(), arg.begin(), arg.end());
    if (result.size() > n_) {
      std::partial_sort(result.begin(), result.begin() + n_, result.end(),
                        std::greater<element_type>());
      result.erase(result.begin() + n_, result.end());
    }
  }

  // Reduce an argument
  void operator()(result_type& result, const argument_type& arg) const {
    // make sure n does not exceed number of elements in arg
    std::size_t n = std::min(n_, arg.size());
    result_type arg_result(n);
    std::partial_sort_copy(arg.begin(), arg.end(), arg_result.begin(),
                           arg_result.end(), std::greater<element_type>());
    operator()(result, arg_result);
  }

 private:
  std::size_t n_;         // number of max number to return
  element_type default_;  // default initialized value

};  // class MaxNReduction

template <typename Tile>
class AbsMaxNReduction {
 public:
  // typedefs
  using element_type = typename Tile::numeric_type;
  using result_type = std::vector<element_type>;
  using argument_type = Tile;

  // constructor
  AbsMaxNReduction(std::size_t n) : n_(n), default_(0) {}
  ~AbsMaxNReduction() = default;

  // Reduction functions

  // Make an empty result object
  result_type operator()() const { return result_type(n_, default_); }

  // Post process the result
  const result_type& operator()(const result_type& result) const {
    return result;
  }

  // Reduce two result objects
  void operator()(result_type& result, const result_type& arg) const {
    result.insert(result.end(), arg.begin(), arg.end());
    if (result.size() > n_) {
      std::partial_sort(result.begin(), result.begin() + n_, result.end(),
                        AbsMaxCompare<element_type>());
      result.erase(result.begin() + n_, result.end());
    }
  }

  // Reduce an argument
  void operator()(result_type& result, const argument_type& arg) const {
    // make sure n does not exceed number of elements in arg
    std::size_t n = std::min(n_, arg.size());
    result_type arg_result(n);
    std::partial_sort_copy(arg.begin(), arg.end(), arg_result.begin(),
                           arg_result.end(), AbsMaxCompare<element_type>());
    operator()(result, arg_result);
  }

 private:
  std::size_t n_;         // number of max number to return
  element_type default_;  // default initialized value

};  // class AbsMaxNReduction

template <typename Tile>
class MaxNIndexReduction {
 public:
  // typedefs
  using element_type = typename Tile::numeric_type;
  using result_type =
      std::vector<std::pair<element_type, std::vector<std::size_t>>>;
  using argument_type = Tile;

  // constructor
  MaxNIndexReduction(std::size_t n)
      : n_(n), default_(std::numeric_limits<element_type>::min()) {}
  ~MaxNIndexReduction() = default;

  // Reduction functions

  // Make an empty result object
  result_type operator()() const {
    return result_type(n_,
                       std::make_pair(default_, std::vector<std::size_t>()));
  }

  // Post process the result
  const result_type& operator()(const result_type& result) const {
    return result;
  }

  // Reduce two result objects
  void operator()(result_type& result, const result_type& arg) const {
    result.insert(result.end(), arg.begin(), arg.end());
    if (result.size() > n_) {
      std::partial_sort(
          result.begin(), result.begin() + n_, result.end(),
          MaxPairCompare<element_type, std::vector<std::size_t>>());
      result.erase(result.begin() + n_, result.end());
    }
  }

  // Reduce an argument
  void operator()(result_type& result, const argument_type& arg) const {
    // make a vector of pairs

    const TA::Range& range = arg.range();

    std::size_t total = arg.size();

    std::size_t start = 0;

    std::vector<std::pair<element_type, std::size_t>> arg_result(total);

    arg_result.shrink_to_fit();

    for (const auto& element : arg) {
      arg_result[start] = std::move(std::make_pair(element, start));
      start++;
    }

    // make sure n does not exceed number of elements in arg
    std::size_t n = std::min(n_, total);

    std::partial_sort(arg_result.begin(), arg_result.begin() + n,
                      arg_result.end(),
                      MaxPairCompare<element_type, std::size_t>());

    arg_result.erase(arg_result.begin() + n, arg_result.end());

    result_type arg_index_result(n);

    for (std::size_t i = 0; i < n; i++) {
      arg_index_result[i] =
          std::make_pair(arg_result[i].first, range.idx(arg_result[i].second));
    }

    operator()(result, arg_index_result);
  }

 private:
  std::size_t n_;         // number of max number to return
  element_type default_;  // default initialized value

};  // class MaxNIndexReduction

template <typename Tile>
class AbsMaxNIndexReduction {
 public:
  // typedefs
  using element_type = typename Tile::numeric_type;
  using result_type =
      std::vector<std::pair<element_type, std::vector<std::size_t>>>;
  using argument_type = Tile;

  // constructor
  AbsMaxNIndexReduction(std::size_t n) : n_(n), default_(0.0) {}
  ~AbsMaxNIndexReduction() = default;

  // Reduction functions

  // Make an empty result object
  result_type operator()() const {
    return result_type(n_,
                       std::make_pair(default_, std::vector<std::size_t>()));
  }

  // Post process the result
  const result_type& operator()(const result_type& result) const {
    return result;
  }

  // Reduce two result objects
  void operator()(result_type& result, const result_type& arg) const {
    result.insert(result.end(), arg.begin(), arg.end());
    if (result.size() > n_) {
      std::partial_sort(
          result.begin(), result.begin() + n_, result.end(),
          AbsMaxPairCompare<element_type, std::vector<std::size_t>>());
      result.erase(result.begin() + n_, result.end());
    }
  }

  // Reduce an argument
  void operator()(result_type& result, const argument_type& arg) const {
    // make a vector of pairs

    const TA::Range& range = arg.range();

    std::size_t total = arg.size();

    std::size_t start = 0;

    std::vector<std::pair<element_type, std::size_t>> arg_result(total);

    arg_result.shrink_to_fit();

    for (const auto& element : arg) {
      arg_result[start] = std::move(std::make_pair(element, start));
      start++;
    }

    // make sure n does not exceed number of elements in arg
    std::size_t n = std::min(n_, total);

    std::partial_sort(arg_result.begin(), arg_result.begin() + n,
                      arg_result.end(),
                      AbsMaxPairCompare<element_type, std::size_t>());

    arg_result.erase(arg_result.begin() + n, arg_result.end());

    result_type arg_index_result(n);

    for (std::size_t i = 0; i < n; i++) {
      arg_index_result[i] =
          std::make_pair(arg_result[i].first, range.idx(arg_result[i].second));
    }

    operator()(result, arg_index_result);
  }

 private:
  std::size_t n_;         // number of max number to return
  element_type default_;  // default initialized value

};  // class MaxNIndexReduction

}  // namespace detail

}  // namespace mpqc

#endif  // SRC_MPQC_MATH_EXTERNAL_TILEDARRAY_REDUCTION_H_

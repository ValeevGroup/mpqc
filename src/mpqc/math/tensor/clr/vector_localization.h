#ifndef MPQC4_SRC_MPQC_MATH_TENSOR_CLR_VECTOR_LOCALIZATION_H_
#define MPQC4_SRC_MPQC_MATH_TENSOR_CLR_VECTOR_LOCALIZATION_H_

#include <random>
#include <utility>
#include <vector>

#include <tiledarray.h>

#include "mpqc/math/external/eigen/eigen.h"

namespace mpqc {
namespace tensor {
template <typename T>
using Matrix =
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

template <typename T>
using Vector = Eigen::Matrix<T, Eigen::Dynamic, 1>;

template <typename T>
class VectorCluster {
 private:
  using VecIdx = std::pair<Vector<T>, long>;
  std::vector<VecIdx> elems_;
  Vector<double> center_;

  void update_center() {
    Vector<double> new_center = Vector<double>::Zero(center_.size());
    for (auto const &elem : elems_) {
      new_center += elem.first;
    }
    if (!elems_.empty()) {
      new_center /= elems_.size();
      center_ = new_center;
    }
  }

  void clear_elems() { elems_.clear(); }

 public:
  VectorCluster() = default;
  VectorCluster(VectorCluster const &) = default;
  VectorCluster(VectorCluster &&) = default;
  VectorCluster &operator=(VectorCluster const &) = default;
  VectorCluster &operator=(VectorCluster &&) = default;

  VectorCluster(Vector<T> const &center) : elems_(0), center_(center) {}

  Vector<double> const &center() const { return center_; }

  void add_elem(Vector<T> const &vec, long idx) {
    elems_.push_back(std::make_pair(vec, idx));
  }

  T sum_of_norms() {
    T sum = 0;
    for (auto const &elem : elems_) {
      sum += (elem.first - center_).norm();
    }
    return sum;
  }

  std::vector<VecIdx> const &elems() const { return elems_; }

  bool empty() const { return elems_.empty(); }

  void update() {
    update_center();
    clear_elems();
  }
};

template <typename T>
void attach_to_closest(Matrix<T> const &D,
                       std::vector<VectorCluster<T>> &clusters,
                       bool init = false) {
  if (!init) {
    for (auto &cluster : clusters) {
      cluster.update();
    }
  }

  for (auto i = 0; i < D.rows(); ++i) {
    Vector<T> current = D.row(i);
    auto &closest = *std::min_element(
        clusters.begin(), clusters.end(),
        [&current](VectorCluster<T> const &a, VectorCluster<T> const &b) {
          return (a.center() - current).squaredNorm() <
                 (b.center() - current).squaredNorm();
        });

    closest.add_elem(current, i);
  }
}

template <typename T>
std::vector<VectorCluster<T>> init_rows(Matrix<T> const &D,
                                        unsigned long num_clusters) {
  std::vector<double> weights(D.rows(), 1.0);
  std::minstd_rand engine(42);

  std::vector<VectorCluster<T>> clusters(num_clusters);

  for (auto i = 0ul; i < num_clusters; ++i) {
    std::discrete_distribution<unsigned int> random_index(weights.begin(),
                                                          weights.end());

    auto idx = random_index(engine);
    clusters[i] = VectorCluster<T>(D.row(idx));

    weights[idx] = 0;

    // for (auto j = 0; j < D.rows(); ++j) {
    //     Vector<T> current = D.row(j);

    //     // Need to search for this
    //     auto closest_cluster = 0;
    //     for(

    //     auto smallest_norm
    //           = (clusters[closest_cluster].center() - current).norm();

    //     for (auto k = 1ul; k <= i; ++k) {
    //         auto diff_norm = (current - clusters[k].center()).norm();
    //         if (diff_norm < smallest_norm) {
    //             smallest_norm = diff_norm;
    //         }
    //     }

    //     weights[j] = smallest_norm * smallest_norm;
    // }
  }
  bool initial_clustering = true;
  attach_to_closest(D, clusters, initial_clustering);
  return clusters;
}

template <typename T>
Vector<unsigned long> get_pivots(
    std::vector<VectorCluster<T>> const &clusters) {
  std::vector<unsigned long> pivs;
  for (auto const &cluster : clusters) {
    for (auto const &elem : cluster.elems()) {
      pivs.push_back(elem.second);
    }
  }

  Vector<unsigned long> J =
      Eigen::Map<Vector<unsigned long>>(pivs.data(), pivs.size(), 1);
  return J;
}

template <typename T>
auto all_have_elems(std::vector<T> &vs) {
  auto end = vs.end();
  for (auto it = vs.begin(); it != end; ++it) {
    if (it->elems().size() == 0) {
      return it;
    }
  }
  return end;
}

template <typename T>
TA::TiledRange1 localize_vectors_with_kmeans(Matrix<T> const &xyz, Matrix<T> &D,
                                             unsigned long num_clusters) {
  auto clusters = init_rows(xyz, num_clusters);
  for (auto i = 0; i < 5; ++i) {
    attach_to_closest(xyz, clusters);

    for (auto j = 0; j < clusters.size(); ++j) {
      if (clusters[j].elems().empty()) {
        clusters[j] = VectorCluster<T>(xyz.row(j));
        attach_to_closest(xyz, clusters);
      }
    }
  }

  for (auto i = 0; i < clusters.size(); ++i) {
    if (clusters[i].elems().empty()) {
      std::cout << "Cluster " << i << " was emepty with center "
                << clusters[i].center().transpose() << std::endl;
      std::cout << "Other clusters had\n";
      for (auto j = 0; j < clusters.size(); ++j) {
        if (i != j) {
          std::cout << "\t" << j << ": " << clusters[j].elems().size()
                    << std::endl;
        }
      }
      throw std::runtime_error("Empty cluster");
    }
  }

  Vector<unsigned long> J = get_pivots(clusters);
  Eigen::PermutationWrapper<Vector<unsigned long>> P(J);
  P.applyThisOnTheRight(D);

  std::vector<unsigned long> blocks(1, 0);
  blocks.reserve(num_clusters);

  for (auto const &cluster : clusters) {
    assert(!cluster.empty());
    blocks.push_back(cluster.elems().size() + blocks.back());
  }

  return TA::TiledRange1(blocks.begin(), blocks.end());
}

}  // namespace tensor
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_MATH_TENSOR_CLR_VECTOR_LOCALIZATION_H_

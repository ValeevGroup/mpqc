/* Cluster concept
 *
 * 2016 Drew Lewis
 */
#pragma once
#ifndef MPQC_CLUSTERING_CLUSTER_CONCEPT_H
#define MPQC_CLUSTERING_CLUSTER_CONCEPT_H

#include <vector>
#include <memory>
#include <type_traits>

#include "../include/eigen.h"

#include "collapse_to_base.h"

namespace mpqc {
namespace clustering {

template <typename Base>
class ClusterConcept {
 public:
  virtual ~ClusterConcept() = default;

  virtual Eigen::Vector3d center_() const = 0;
  virtual std::vector<Base> collapse_to_base_() const = 0;
};

template <typename T, typename Base>
class ClusterModel : public ClusterConcept<Base> {
 private:
  T element_;

 public:
  ClusterModel(T t) : element_(std::move(t)) {}
  ClusterModel(ClusterModel const &c) = default;
  ClusterModel &operator=(ClusterModel c) {
    element_ = std::move(c.element_);
    return *this;
  }

  ClusterModel(ClusterModel &&c) = default;
  ClusterModel &operator=(ClusterModel &&c) = default;

  Eigen::Vector3d center_() const override final { return center(element_); }

  std::vector<Base> collapse_to_base_() const override final {
    return collapse_to_base<Base>{}.template
    operator()<typename std::remove_const<decltype(element_)>::type>(element_);
  }
};

template <typename Base>
class Clusterable {
 private:
  std::shared_ptr<const ClusterConcept<Base>> element_impl_;

 public:
  template <typename T>
  Clusterable(T t)
      : element_impl_(std::make_shared<ClusterModel<T, Base>>(std::move(t))) {}

  Clusterable(Clusterable const &c) = default;
  Clusterable(Clusterable &&c) = default;
  Clusterable &operator=(Clusterable const &c) = default;
  Clusterable &operator=(Clusterable &&c) = default;

  Eigen::Vector3d center() const { return element_impl_->center_(); }

  std::vector<Base> collapse_to_base() const {
    return element_impl_->collapse_to_base_();
  }
};

template<typename T>
Eigen::Vector3d center(Clusterable<T> const &c){
    return c.center();
}

template <typename Base>
class Cluster {
 private:
  std::vector<Clusterable<Base>> cluster_;

 public:
  Cluster(std::vector<Clusterable<Base>> cluster)
      : cluster_(std::move(cluster)) {}

  std::vector<Base> flatten() const {
    std::vector<Base> bases;
    for (auto const &c : cluster_) {
      auto temp = c.collapse_to_base();
      bases.insert(bases.end(), temp.begin(), temp.end());
    }

    return bases;
  }

  auto begin() const -> decltype(cluster_.cbegin()) {
    return cluster_.cbegin();
  }

  auto end() const -> decltype(cluster_.cend()) { return cluster_.cend(); }

  auto size() const -> decltype(cluster_.size()) { return cluster_.size(); }

};

template<typename T>
Eigen::Vector3d center(Cluster<T> const &c){
    Eigen::Vector3d center;
    center.setZero();
    for(auto const &e : c){
        auto e_center = e.center();
        center[0] += e_center[0];
        center[1] += e_center[1];
        center[2] += e_center[2];
    }

    center /= double(c.size());
    return center;
}

template <typename Base>
class Clustered {
 private:
  std::vector<Cluster<Base>> clusters_;

 public:
  Clustered(std::vector<Cluster<Base>> clusters)
      : clusters_(std::move(clusters)) {}

  std::vector<std::vector<Base>> flatten_clusters() const {
    std::vector<Base> bases;
    for (auto const &c : clusters_) {
      auto temp = c.flatten();
      bases.insert(bases.end(), temp.begin(), temp.end());
    }

    return bases;
  }
};

}  // namespace clustering
}  // namespace mpqc

#endif  // MPQC_CLUSTERING_CLUSTER_CONCEPT_H

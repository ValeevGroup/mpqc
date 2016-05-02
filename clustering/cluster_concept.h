/* Cluster concept
 *
 * 2016 Drew Lewis
 */
#pragma once
#ifndef MPQC_CLUSTERING_CLUSTER_CONCEPT_H
#define MPQC_CLUSTERING_CLUSTER_CONCEPT_H

#include <vector>
#include <memory>

#include "../include/eigen.h"

#include "mpqc/util/misc/type_traits.h"
#include "collapse_to_base.h"

namespace mpqc {
namespace clustering {

namespace detail {

template <typename, class = void_t<>>
struct can_compute_mass_ : std::false_type {};

template <typename T>
struct can_compute_mass_<T, void_t<decltype(mass(std::declval<T>()))>>
    : std::true_type {};

template <typename, class = void_t<>>
struct can_compute_com_ : std::false_type {};

template <typename T>
struct can_compute_com_<T, void_t<decltype(center_of_mass(std::declval<T>()))>>
    : can_compute_mass_<T> {};

template <typename Base>
struct ClusterConceptAbstract {
  virtual ~ClusterConceptAbstract() = default;
  virtual Eigen::Vector3d center_() const = 0;
  virtual std::vector<Base> collapse_to_base_() const = 0;
};

template <typename Base>
struct ContainsCom : ClusterConceptAbstract<Base> {
  virtual Eigen::Vector3d com_() const = 0;
  virtual double mass_() const = 0;
};

}  // namespace detail

template <typename Base, typename = void_t<>>
class ClusterConcept {};

template <typename Base>
class ClusterConcept<Base, enable_if_t<detail::can_compute_com_<Base>::value>>
    : public detail::ContainsCom<Base> {};

template <typename Base>
class ClusterConcept<Base, enable_if_t<!detail::can_compute_com_<Base>::value>>
    : public detail::ClusterConceptAbstract<Base> {};

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

  Eigen::Vector3d center_() const { return center(element_); }

  Eigen::Vector3d com_() const {
    static_assert(
        detail::can_compute_com_<T>::value,
        "ClusterModel T class must define center_of_mass external function.");
    return center_of_mass(element_);
  }

  double mass_() const {
    static_assert(
        detail::can_compute_mass_<T>::value,
        "ClusterModel T class must define center_of_mass external function.");
    return mass(element_);
  }

  std::vector<Base> collapse_to_base_() const {
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

  Eigen::Vector3d com() const {
    static_assert(
        detail::can_compute_com_<Base>::value,
        "Clusterable Base class must define center_of_mass external function.");
    return element_impl_->com_();
  }

  double mass() const {
    static_assert(detail::can_compute_mass_<Base>::value,
                  "Clusterable Base class must define mass external function.");
    return element_impl_->mass_();
  }

  std::vector<Base> collapse_to_base() const {
    return element_impl_->collapse_to_base_();
  }
};

template <typename T>
Eigen::Vector3d center(Clusterable<T> const &c) {
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

template <typename T>
Eigen::Vector3d center(Cluster<T> const &c) {
  Eigen::Vector3d center;
  center.setZero();
  for (auto const &e : c) {
    auto e_center = e.center();
    center[0] += e_center[0];
    center[1] += e_center[1];
    center[2] += e_center[2];
  }

  center /= double(c.size());
  return center;
}

template <typename T,
          enable_if_t<detail::can_compute_com_<T>::value> * = nullptr>
Eigen::Vector3d center_of_mass(Cluster<T> const &c) {
  Eigen::Vector3d com;
  com.setZero();
  auto total_mass = 0.0;

  for (auto const &e : c) {
    auto e_center = e.center();
    auto e_mass = e.mass();
    com[0] += e_mass * e_center[0];
    com[1] += e_mass * e_center[1];
    com[2] += e_mass * e_center[2];
    total_mass += e_mass;
  }

  com /= total_mass;
  return com;
}

template <typename T,
          enable_if_t<detail::can_compute_com_<T>::value> * = nullptr>
double mass(Cluster<T> const &c) {
  auto total_mass = 0.0;
  for (auto const &e : c) {
    total_mass += e.mass();
  }

  return total_mass;
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

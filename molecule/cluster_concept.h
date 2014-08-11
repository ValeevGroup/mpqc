#ifndef CLUSTER_CONCEPT_H
#define CLUSTER_CONCEPT_H

#include <memory>
#include <vector>
#include "../include/eigen.h"
#include "molecule_fwd.h"
#include "cluster_collapse.h"

/**
 * @brief ClusterConcept is the base class which defines the operations which
 * different clusterable types must have.
 */
class ClusterConcept {
public:
  virtual ~ClusterConcept() = default;

  virtual ClusterConcept *copy() const = 0;
  virtual Eigen::Vector3d center() const = 0;
  virtual double mass() const = 0;
  virtual double charge() const = 0;
  virtual std::vector<Atom> atoms() const = 0;
};

/**
 * @brief ClusterModel is a class which is basically used for type erasure
 */
template <typename T> class ClusterModel : public ClusterConcept {
public:
  ClusterModel(T t) : element_(std::move(t)) {}
  ClusterModel(const ClusterModel &c) = default;

  ClusterModel &operator=(ClusterModel c) {
    element_ = std::move(c.element_);
    return *this;
  }

  ClusterModel(ClusterModel &&c) = default;
  ClusterModel &operator=(ClusterModel &&c) = default;

  ClusterConcept *copy() const { return new ClusterModel(*this); }

  Eigen::Vector3d center() const { return element_.center(); }
  double mass() const { return element_.mass(); }
  double charge() const { return element_.charge(); }
  std::vector<Atom> atoms() const { return collapse_to_atoms(element_); }

private:
  T element_;
};

/**
 * @brief The Clusterable is a class that holds any clusterable type.
 * The requirements on this type are set by ClusterConcept.
 */
class Clusterable {
public:
  using position_t = Eigen::Vector3d;

  template <typename T>
  Clusterable(T t)
      : element_impl_(new ClusterModel<T>(std::move(t))) {}

  Clusterable(const Clusterable &c) : element_impl_(std::move(c.element_impl_->copy())) {}

  // for operator make a copy and then move that copy into this.
  Clusterable &operator=(const Clusterable &c) {
    *this = std::move(Clusterable(c));
    return *this;
  }

  Clusterable(Clusterable &&c) = default;
  Clusterable &operator=(Clusterable &&c) = default;

  position_t center() const { return element_impl_->center(); }
  double mass() const { return element_impl_->mass(); }
  double charge() const { return element_impl_->charge(); }
  std::vector<Atom> atoms() const { return element_impl_->atoms(); }

private:
  std::unique_ptr<ClusterConcept> element_impl_;
};

#endif // CLUSTER_CONCEPT_H

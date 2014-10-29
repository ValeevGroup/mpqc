#pragma once
#ifndef TCC_MOLECULE_CLUSTERCONCEPT_H
#define TCC_MOLECULE_CLUSTERCONCEPT_H

#include "molecule_fwd.h"
#include "cluster_collapse.h"

#include <memory>
#include <vector>

namespace tcc {
namespace molecule {

/**
 * @brief ClusterConcept is the base class which defines the operations which
 * different clusterable types must have.
 */
class ClusterConcept {
  public:
    virtual ~ClusterConcept() = default;

    virtual ClusterConcept *clone() const = 0;
    virtual Eigen::Vector3d center() const = 0;
    virtual double mass() const = 0;
    virtual double charge() const = 0;
    virtual std::vector<Atom> atoms() const = 0;
};

/**
 * @brief ClusterModel is a class which is basically used for type erasure
 */
template <typename T>
class ClusterModel : public ClusterConcept {
  public:
    ClusterModel(T t) : element_(std::move(t)) {}
    ClusterModel(const ClusterModel &c) = default;
    ClusterModel &operator=(ClusterModel c) {
        element_ = std::move(c.element_);
        return *this;
    }
    ClusterModel(ClusterModel &&c) = default;
    ClusterModel &operator=(ClusterModel &&c) = default;

    ClusterConcept *clone() const override final {
        return new ClusterModel(*this);
    }

    Eigen::Vector3d center() const override final { return element_.center(); }
    double mass() const override final { return element_.mass(); }
    double charge() const override final { return element_.charge(); }
    std::vector<Atom> atoms() const override final {
        return collapse_to_atoms(element_);
    }

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
        : element_impl_(std::make_shared<ClusterModel<T>>(std::move(t))) {}
    Clusterable(Clusterable const &c) = default;
    Clusterable &operator=(Clusterable const &c) = default;
    Clusterable(Clusterable &&c) = default;
    Clusterable &operator=(Clusterable &&c) = default;

    position_t center() const { return element_impl_->center(); }
    double mass() const { return element_impl_->mass(); }
    double charge() const { return element_impl_->charge(); }
    std::vector<Atom> atoms() const { return element_impl_->atoms(); }

  private:
    std::shared_ptr<const ClusterConcept> element_impl_;
};

} // namespace molecule
} // namespace tcc

#endif // TCC_MOLECULE_CLUSTERCONCEPT_H

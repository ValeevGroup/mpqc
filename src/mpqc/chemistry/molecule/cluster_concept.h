#pragma once
#ifndef MPQC_MOLECULE_CLUSTERCONCEPT_H
#define MPQC_MOLECULE_CLUSTERCONCEPT_H

#include <mpqc/chemistry/molecule/molecule_fwd.h>

#include <memory>
#include <iostream>

namespace mpqc {

/*!
 * \addtogroup Molecule
 *
 * @{
 */

/*
 * ClusterConcept is the base class which defines the operations which
 * different clusterable types must have.
 */
class ClusterConcept {
  public:
    virtual ~ClusterConcept() = default;

    virtual ClusterConcept *clone_() const = 0;
    virtual Vec3D const & center_() const = 0;
    virtual std::ostream &print_(std::ostream &) const = 0;
};

/*
 * ClusterModel is a class which is basically used for type erasure
 */
template <typename T>
class ClusterModel : public ClusterConcept {
  private:
    T element_;

  public:
    ClusterModel(T t) : element_(std::move(t)) {}
    ClusterModel(const ClusterModel &c) = default;
    ClusterModel &operator=(ClusterModel c) {
        element_ = std::move(c.element_);
        return *this;
    }

    ClusterModel(ClusterModel &&c) = default;
    ClusterModel &operator=(ClusterModel &&c) = default;

    ClusterConcept *clone_() const override final {
        return new ClusterModel(*this);
    }

    Vec3D const &center_() const override final { return center(element_); }
    std::ostream &print_(std::ostream &os) const override final {
        os << element_;
        return os;
    }
};

/*!
 * \brief The Clusterable is a class that holds any clusterable type.
 *
 * The any type which fits the Clusterable concept must provide the following
 * functions
 */
class Clusterable {
  private:
    std::shared_ptr<const ClusterConcept> element_impl_;

  public:
    template <typename C>
    Clusterable(C c)
            : element_impl_(std::make_shared<ClusterModel<C>>(std::move(c))) {}
    Clusterable(Clusterable const &c) = default;
    Clusterable &operator=(Clusterable const &c) = default;
    Clusterable(Clusterable &&c) = default;
    Clusterable &operator=(Clusterable &&c) = default;

    Vec3D const &center() const { return element_impl_->center_(); }

    std::ostream &print(std::ostream &os) const {
        return element_impl_->print_(os);
    }
};

/*! @} */

} // namespace mpqc

#endif // MPQC_MOLECULE_CLUSTERCONCEPT_H

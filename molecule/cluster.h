#ifndef CLUSTER_H
#define CLUSTER_H

#include <vector>
#include <numeric>
#include "../include/eigen.h"
#include "cluster_concept.h"

//TODO move code to cpp file
class Cluster {
public:
  using position_t = Clusterable::position_t;
  Cluster() = default;

  Cluster(const Cluster &c) = default;
  Cluster& operator=(const Cluster &c) = default;

  // Eigen doesn't have move constructors so use their swap.
  Cluster(Cluster &&c) noexcept : elements_(std::move(c.elements_)), charge_(std::move(c.charge_)), mass_(std::move(c.mass_)){
    c.center_.swap(center_);
  }
  Cluster& operator=(Cluster &&c) {
    elements_ = std::move(c.elements_);
    charge_ = std::move(c.charge_);
    mass_ = std::move(c.mass_);
    c.center_.swap(center_);
    return *this;
  }

  template <typename T> void add_clusterable(T t) {
    mass_ += t.mass();
    charge_ += t.charge();
    elements_.emplace_back(t);
  }

  unsigned long nelements() const {return elements_.size();}

  void clear() {
    charge_ = 0.0;
    mass_ = 0.0;
    elements_.clear();
  }

  void init_center(position_t point) { center_ = point; }

  void guess_center() {
    center_ = position_t::Zero();
    for(const auto element : elements_){
      center_ += element.center() * element.mass();
    }
    center_ /= mass_;
  }

  double sum_distances_from_center() const {
    return std::accumulate(elements_.begin(), elements_.end(),0.0,[&](
                           double d, const Clusterable &c){
      return d + (c.center() - center_).norm();
    });
  }

  position_t center() const { return center_; }
  double mass() const { return mass_; }
  double charge() const { return charge_; }

  std::vector<const Clusterable>::iterator begin() const {
    return elements_.begin();
  }

  std::vector<const Clusterable>::iterator end() const {
    return elements_.end();
  }

private:
  std::vector<Clusterable> elements_;
  position_t center_ = {0,0,0};
  double charge_ = 0.0;
  double mass_ = 0.0;
};

#endif // CLUSTER_H

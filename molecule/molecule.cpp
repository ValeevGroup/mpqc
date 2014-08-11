
#include "molecule.h"
#include "cluster_concept.h"
#include "cluster.h"

namespace molecule_detail {
inline double calculate_mass(const std::vector<Clusterable> &c) {
  return std::accumulate(
      c.begin(), c.end(), 0.0,
      [](double x, const Clusterable &c) { return x + c.mass(); });
}

inline double calculate_charge(const std::vector<Clusterable> &c) {
  return std::accumulate(
      c.begin(), c.end(), 0.0,
      [](double x, const Clusterable &c) { return x + c.charge(); });
}

inline Clusterable::position_t
calculate_center(const std::vector<Clusterable> &c, double mass) {
  Clusterable::position_t center = Clusterable::position_t::Zero();
  for (const auto &elem : c) {
    center += elem.center() * elem.mass();
  }
  return center /= mass;
}

class sort_by_distance_from_point {
public:
  sort_by_distance_from_point(const Cluster::position_t &point)
      : point_(point) {}

  bool operator()(const Clusterable &a, const Clusterable &b) {
    Atom::position_t a_dist = a.center() - point_;
    Atom::position_t b_dist = b.center() - point_;
    if (!(a_dist.norm() == b_dist.norm())) {
      return a_dist.norm() < b_dist.norm();
    } else if (a_dist[0] == b_dist[0]) {
      return a_dist[0] < b_dist[0];
    } else if (a_dist[1] == b_dist[1]) {
      return a_dist[1] < b_dist[1];
    } else
      return a_dist[2] < b_dist[2];
  }

private:
  Cluster::position_t point_;
};

} // namespace moleucle detail

Molecule::Molecule(std::vector<Clusterable> c) : elements_(std::move(c)) {
  mass_ = molecule_detail::calculate_mass(elements_);
  charge_ = molecule_detail::calculate_charge(elements_);
  center_ = molecule_detail::calculate_center(elements_, mass_);
  std::sort(elements_.begin(), elements_.end(),
            molecule_detail::sort_by_distance_from_point(center_));
}

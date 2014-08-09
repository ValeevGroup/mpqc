#ifndef ATOM_H
#define ATOM_H

#include "../include/eigen.h"

class Atom {
public:
  using position_t = Eigen::Vector3d;

  Atom() = default;
  Atom(const Atom &atom) = default;
  Atom(Atom &&atom) = default;
  Atom &operator=(const Atom &atom) = default;
  Atom &operator=(Atom &&atom) = default;
  Atom(position_t center, double mass, double charge)
      : center_(std::move(center)),
        charge_(std::move(charge)),
        mass_(std::move(mass)) {}

  position_t center() const { return center_; }
  double charge() const {return charge_;}
  double mass() const {return mass_;}

private:
  position_t center_ = position_t::Zero();
  double charge_ = 0;
  double mass_ = 0;
};

#endif // ATOM_H

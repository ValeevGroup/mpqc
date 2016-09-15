#include "cluster_concept.h"

#include <vector>
#include <iostream>

struct Atom {
  Eigen::Vector3d r;
  double mass;
  Atom(double m, Eigen::Vector3d xyz) : r(xyz), mass(m) {}
};

struct Shell {
  Eigen::Vector3d r;
  Shell(Eigen::Vector3d xyz) : r(xyz) {}
};

Eigen::Vector3d center(Shell const &s) { return s.r; }
Eigen::Vector3d center_of_mass(Shell const &s) { return s.r; }

Eigen::Vector3d center(Atom const &a) { return a.r; }

Eigen::Vector3d center_of_mass(Atom const &a) { return a.r; }

double mass(Atom const &a) { return a.mass; }

struct Molecule {
  std::vector<Atom> atoms;
  auto begin() const -> decltype(atoms.cbegin()) { return atoms.cbegin(); }
  auto end() const -> decltype(atoms.cend()) { return atoms.cend(); }
};

Eigen::Vector3d center(Molecule const &m) {
  Eigen::Vector3d c;
  c.setZero();

  for (auto const &a : m.atoms) {
    c[0] += a.r[0];
    c[1] += a.r[1];
    c[2] += a.r[2];
  }

  c /= double(m.atoms.size());

  return c;
}

double mass(Molecule const &m) {
  auto total_mass = 0.0;
  for (auto const &a : m.atoms) {
    total_mass += a.mass;
  }

  return total_mass;
}

Eigen::Vector3d center_of_mass(Molecule const &m) {
  Eigen::Vector3d c;
  c.setZero();

  auto total_mass = 0.0;
  for (auto const &a : m.atoms) {
    c[0] += a.mass * a.r[0];
    c[1] += a.mass * a.r[1];
    c[2] += a.mass * a.r[2];
    total_mass += a.mass;
  }

  c /= total_mass;

  return c;
}

int main() {
  Shell sh({0, 0, 0});

  std::vector<Atom> atoms;
  for (auto i = 0; i < 10; ++i) {
    atoms.push_back(Atom(i, {i, 0, 2}));
  }

  Molecule mol;
  mol.atoms = atoms;

  // Clusterables of atom and molecule
  mpqc::clustering::Clusterable<Shell> s1(sh);
  s1.mass();
  mpqc::clustering::Clusterable<Atom> a1(atoms[0]);
  mpqc::clustering::Clusterable<Atom> m1(mol);
  std::cout << "Center of mass of the molecule is: " << m1.com().transpose()
            << std::endl;

  // cluster of 2 clusterables
  mpqc::clustering::Cluster<Atom> c1({a1, m1});

  // cluster of 2 clusterables one of which is a cluster
  mpqc::clustering::Cluster<Atom> c2({a1, c1});

  // clusterable of clusterable
  mpqc::clustering::Clusterable<Atom> cluster_ception(m1);
  std::cout << "Center of mass of the clusterable molecule is: "
            << cluster_ception.com().transpose() << std::endl;

  mpqc::clustering::Clusterable<Atom> cluster_super_ception(c2);
  mpqc::clustering::Cluster<Atom> c3({cluster_super_ception, m1, a1, c1});
  std::cout << "Center of mass of everything: "
            << center_of_mass(c3).transpose() << std::endl;

  auto all_atoms = c3.flatten();

  for (auto const &a : all_atoms) {
    std::cout << a.r.transpose() << std::endl;
  }

  return 0;
}

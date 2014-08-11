#ifndef CLUSTER_COLLAPSE_H
#define CLUSTER_COLLAPSE_H

#include "atom.h"

/**
 * @brief collapse_to_atoms ends the recursive loop of the templated version of
 * the function
 */
inline std::vector<Atom> collapse_to_atoms(const Atom &a) {
  return {a};
}

/**
 * collapse_to_atoms must have a begin and end iterator which reference
 * a vector of Clusterables
 */
template <typename T> std::vector<Atom> collapse_to_atoms(const T &t) {
  std::vector<Atom> atoms;
  for (const auto &element : t) {
    std::vector<Atom> temp_atoms = element.atoms();
    atoms.insert(atoms.end(), temp_atoms.begin(), temp_atoms.end());
  }
  return atoms;
}

#endif // CLUSTER_COLLAPSE_H

#pragma once
#ifndef TCC_MOLECULE_MOLECULE_H
#define TCC_MOLECULE_MOLECULE_H

#include "molecule_fwd.h"
#include <functional>
#include <vector>

namespace tcc {
namespace molecule {

class Molecule {
  public:
    using cluster_fn_t = std::function<std::vector<Cluster>(
        std::vector<Clusterable>, unsigned long)>;

    Molecule(std::vector<Clusterable> c);

    position_t center() const;
    double charge() const;
    double mass() const;

    std::vector<Clusterable>::const_iterator begin() const;
    std::vector<Clusterable>::const_iterator end() const;

    unsigned long nelements() const;

    std::vector<Cluster>
    cluster_molecule(cluster_fn_t fn, unsigned long nclusters) const;

  private:
    std::vector<Clusterable> elements_;

    position_t center_ = {0, 0, 0};
    double mass_ = 0.0;
    double charge_ = 0.0;
};

} // namespace molecule
} // namespace tcc

#endif // TCC_MOLECULE_MOLECULE_H

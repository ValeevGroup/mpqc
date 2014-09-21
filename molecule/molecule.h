#ifndef CLUSTERED_MOLECULE_H
#define CLUSTERED_MOLECULE_H

#include "molecule_fwd.h"
#include <functional>
#include <vector>

class Molecule {
  public:
    using cluster_fn_t = std::function
        <std::vector<Cluster>(std::vector<Clusterable>, unsigned long)>;

    Molecule(std::vector<Clusterable> c);

    position_t center() const { return center_; }
    double charge() const { return charge_; }
    double mass() const { return mass_; }

    std::vector<const Clusterable>::iterator begin() const {
        return elements_.begin();
    }
    std::vector<const Clusterable>::iterator end() const {
        return elements_.end();
    }

    unsigned long nelements() const { return elements_.size(); }

    std::vector<Cluster>
    cluster_molecule(cluster_fn_t fn, unsigned long nclusters) const {
        return fn(elements_, nclusters);
    }

  private:
    std::vector<Clusterable> elements_;

    position_t center_ = {0, 0, 0};
    double mass_ = 0.0;
    double charge_ = 0.0;
};

#endif // CLUSTERED_MOLECULE_H

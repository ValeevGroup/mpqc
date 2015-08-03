#pragma once
#ifndef TCC_MOLECULE_MOLECULE_H
#define TCC_MOLECULE_MOLECULE_H

#include "molecule_fwd.h"
#include <functional>
#include <vector>

namespace tcc {
namespace molecule {

/*!
 * \defgroup Molecule Molecule
 *
 * \brief The molecule module contains information about how to make and cluster
 *molecules
 *
 *
 *  The molecule module contains all of the classes which are needed for 
 *  clustering.  Ulitimately everything in the molecule module should support 
 *  the Clusterable interface.
 *  
 * @{
 */

/*! \brief Molecule is a class which contains a vector of clusterables
 *
 * At its core molecule is a collection of things which can be collapsed
 * to atoms.  Its main job is allow for clustering of its clusterables.
 *
 */
class Molecule {
  public:
    using cluster_fn_t
          = std::function<std::vector<Cluster>(std::vector<Clusterable>,
                                               unsigned long)>;

    Molecule(std::vector<Clusterable> c);

    position_t center() const;
    int charge() const;
    double mass() const;

    std::vector<Clusterable>::const_iterator begin() const;
    std::vector<Clusterable>::const_iterator end() const;

    unsigned long nelements() const;
    unsigned long occupation(unsigned long total_charge) const;
    // auto compute core electrons
    unsigned long core_electrons() const;
    double nuclear_repulsion() const;

    std::vector<Cluster>
    cluster_molecule(cluster_fn_t fn, unsigned long nclusters) const;

    // Will attach any hydrogens to their closest heavy atom.
    std::vector<Cluster> attach_hydrogens() const;

    /// Will attach hydrogens by calling attach_hydrogens and will then do a
    /// search for the best k-means clustering by attempting to cluster multiple
    /// with different seeds.  The criteria for best is determined by taking the
    /// minimimum sum of the square distances from each cluster. Any groupings
    /// that include a cluster with zero memebers are thrown away.
    std::vector<Cluster>
    attach_H_and_kmeans(unsigned long nclusters,
                        unsigned long init_seed = 42) const;

    /// Kmeans without h-attach
    std::vector<Cluster>
    kmeans(unsigned long nclusters,
                        unsigned long init_seed = 42) const;

  private:
    std::vector<Clusterable> elements_;

    position_t center_ = {0, 0, 0};
    double mass_ = 0.0;
    int charge_ = 0;
};

Molecule read_xyz(std::string const &);

/*!
 * @}
 */

} // namespace molecule
} // namespace tcc

#endif // TCC_MOLECULE_MOLECULE_H

//
// kmeans.hpp
//
// Copyright (C) 2013 Drew Lewis
//
// Authors: Drew Lewis
// Maintainer: Drew Lewis and Edward Valeev
//
// This file is part of the MPQC Toolkit.
//
// The MPQC Toolkit is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The MPQC Toolkit is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the MPQC Toolkit; see the file COPYING.LIB.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//
// The U.S. Government is granted a limited license as per AL 91-7.
//

#ifndef MPQC_BASIS_KMEANS_HPP
#define MPQC_BASIS_KMEANS_HPP

#include<vector>
#include<string>
#include<cassert>

#include<Eigen/Dense>
#include<chemistry/moleucle/atom.h>
#include<chemistry/molecule/molecule.h>
#include<chemistry/qc/basis/basis.h>

#include <mpqc/utility/foreach.hpp>
#include "kcluster.hpp"

namespace mpqc{
namespace basis{

    /**
     * Determines the clustering of shells based on k-means clustering.
     */
    class ShellOrder {
        using namespace TA = TiledArray;
        using Atom = cluster::Atom;
        using Shell = sc::Shell;

    public:

    /**
     * Initializes ShellOrder with the atoms from the molecule and the shells
     * from the basis.
     */
    ShellOrder(const sc::Ref<sc::GausianBasisSet> &basis) :
        clusters_(),
        atoms_(),
        basis_(basis)
    {
        // Get the molecule.
        sc::Ref<sc::Molecule> mol = basis->molecule();

        // Get the Atoms and their index out of the molecule for the clusters.
        for(auto i = 0; i < mol->natom(); ++i){
           atoms_.push_back(Atom(mol->atom(i), i));
        }
    }

    /**
     * Returns a list of shells ordered by which cluster they belong to.
     */
    std::vector<Shell> ordered_shells(std::size_t nclusters){
        // Number of clusters desired.
        nclusters_ = nclusters;

        // Compute clusters using Lloyd's Algorithm
        compute_clusters(nclusters_);

        std::vector<Shell> shells = cluster_shells();


        return shells;
    }

    /**
     * Returns a TiledArray::TiledRange1 based on the number of clusters  that
     * were computed.
     */
    /*
    TA::TiledRange1 trange1(){
        assert(nclusters_ != 0);

    }
    */



    /**
     * Find clusters using Lloyd's algorithm.
     */
    void compute_clusters(std::size_t nclusters){
        // Make initial guess at clusters
        init_clusters();
        // Search for local minimium in terms of clustering.
        k_means_search();
    }

    /**
     * Determines initial guess for clusters.
     */
    void init_clusters(){
        // Sort atoms in molecule by mass
        std::sort(atoms_.begin(), atoms_.end(),
            [](const Atom &a, const Atom &b){return a.mass() > b.mass();});

        // Initialize the kcluster guess at the position of the heaviest atoms.
        for(auto i = 0; i < nclusters_; ++i){
            kclusters_.push_back(
                Vector3(atoms_[i].xyz(0), atoms_[i].xyz(1), atoms_[i],xyz(2))
            );
        }

        // Attach atoms to the cluster which they are closest too.
        attach_to_closest_cluster();
    }

    /**
     * Attaches each atom to its closest cluster.
     */
    void attach_to_closest_cluster(){
        // Loop over all the atoms.
        foreach(const auto atom, atoms_){
            // Guess that first cluster is closest
            double smallest = kclusters_[0].distance(atom);

            // To which cluster the atom belongs.
            std::size_t kindex = 0;

            // Loop over kclusters
            for(auto i = 1; i < nclusters_; ++i){
                // Compute distance from atom to next kcluster
                double dist = kclusters_[i].distance(atom);

                // if closer update index info
                if(dist < smallest){
                    kindex = i;
                    smallest = dist;
                }
            }

            // Add atom to the closest kcluster
            kclusters_[kindex].add_member(atom);
        }
    }

    /// Computes Lloyd's algorith to find a local minimium for the clusters.
    void k_means_search(std::size_t niter = 100){

        // Lloyd's algorithm iterations.
        for(auto i = 0; i < niter; ++i){

           // Recompute the center of the cluster using the centroid
           // the atoms.  Will lose information about which atoms
           // go with which center.
           foreach(auto &cluster, kclusters_){ cluster.guess_center(); }

           attach_to_closest_cluster();
       }
    }

    std::vector<Shell> cluster_shells(){

        std::vector<Shell> shells;
        // Loop over clusters
        foreach(const auto& cluster, clusters_){
            // Loop over atoms in cluster
            foreach(const auto& atom, cluster.atoms()){
                // Figure out where in the molecule the atom is located.
                std::size_t atom_index = atom.mol_index();
                // Figure out how many shells are on the atom.
                std::size_t nshells_on_atom =
                                basis_->nshell_on_center(atom_index);
                // Loop over the shells on the atom and pack them into
                // shells.
                for(auto i = 0; i < nshells_on_atom; ++i){
                    shells.push_back(basis_(atom_index, i));
                }
            }
        }

        return shells;
    }


    private:
        std::size_t nclusters_ = 0;
        std::vector<Atom> atoms_;
        std::vector<KCluster> clusters_;
        sc::Ref<sc::GaussianBasisSet> basis_;
    };

}
}





#endif /* MPQC_BASIS_KMEANS_HPP */

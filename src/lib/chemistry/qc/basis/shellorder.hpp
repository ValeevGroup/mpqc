//
// shellorder.hpp
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

#ifndef CHEMISTRY_QC_BASIS_SHELLORDER_HPP
#define CHEMISTRY_QC_BASIS_SHELLORDER_HPP

#include<vector>
#include<string>
#include<cassert>

#include <Eigen/Dense>
#include <chemistry/molecule/molecule.h>
#include <chemistry/qc/basis/basis.h>

#include <mpqc/utility/foreach.hpp>
#include "kcluster.hpp"

namespace mpqc{
namespace TA{

    /**
     * Determines the clustering of shells based on k-means clustering.
     */
    class ShellOrder {

    public:
        using Shell = sc::GaussianBasisSet::Shell;
        using Atom = KCluster::Atom;
        /// Each element represents the shell a new tile starts on.  So if the
        /// vector looks like | 0, 5, 6 ) then tile 0 is from 0-4, tile 1 has 5,
        /// tile 3 has 6.
        using ShellRange = std::vector<std::size_t>;
        using Vector3 = KCluster::Vector3;

        /**
         * Initializes ShellOrder with the atoms from the molecule and the shells
         * from the basis.
         */
        ShellOrder(const sc::Ref<sc::GaussianBasisSet> &basis) :
            clusters_(),
            atoms_(),
            basis_(basis),
            com_(Vector3::Zero(3))
        {
            // Get the molecule.
            sc::Ref<sc::Molecule> mol = basis->molecule();
            com_[0] = mol->center_of_mass()[0];
            com_[1] = mol->center_of_mass()[1];
            com_[2] = mol->center_of_mass()[2];

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
         * Returns a a ShellRange which specifies what shell each tile starts on.
         */
        ShellRange shell_ranges() const {
            return compute_shell_ranges();
        }

    private:
        /*
         * Find clusters using Lloyd's algorithm.
         */
        void compute_clusters(std::size_t nclusters){
            // Make initial guess at clusters
            init_clusters();
            // Search for local minimium in terms of clustering.
            k_means_search();
        }

        /*
         * Determines initial guess for clusters.
         */
        void init_clusters(){
            // Sort atoms in molecule by distance from COM such that 1. Closest
            // . . . N. Farthest.  If two atoms are equidistant then sort based
            // on x, then y, and finally z.
            auto ordering = [&](const Atom &a, const Atom &b){
                // Get position vector for each atom.
                const Vector3 va = Eigen::Map<const Vector3>(a.r(),3);
                const Vector3 vb = Eigen::Map<const Vector3>(b.r(),3);

                // Find distance from center of mass for each atom
                double da = Vector3(va - com_).norm();
                double db = Vector3(vb - com_).norm();

                // Begin to check coordinates
                if(da != db) { // Check not equidistant
                    return ((da < db) ? true : false);
                }
                else if(va[0] != vb[0]){ // Check x isn't equidistant
                    return ((va[0] < vb[0]) ? true : false);
                }
                else if(va[0] != vb[0]){ // Check y isn't equidistant
                    return ((va[0] < vb[0]) ? true : false);
                }
                else if(va[0] != vb[0]){ // Check z isn't equidistant
                    return ((va[0] < vb[0]) ? true : false);
                }
                else {
                    return false;
                }
            };

            std::sort(atoms_.begin(), atoms_.end(), ordering);

            // Initialize the kcluster guess at the position of atoms closest to
            // COM.
            for(auto i = 0; i < nclusters_; ++i){
                clusters_.push_back(
                    Vector3(atoms_[i].r(0), atoms_[i].r(1), atoms_[i].r(2))
                );
            }

            // Attach atoms to the cluster which they are closest too.
            attach_to_closest_cluster();
        }

        /*
         * Attaches each atom to its closest cluster.
         */
        void attach_to_closest_cluster(){
            // Loop over all the atoms.
            foreach(const auto atom, atoms_){
                // Guess that first cluster is closest
                double smallest = clusters_[0].distance(atom);

                // To which cluster the atom belongs.
                std::size_t kindex = 0;

                // Loop over kclusters using i = 1 because we already computed 0
                for(auto i = 1; i < nclusters_; ++i){
                    // Compute distance from atom to next kcluster
                    double dist = clusters_[i].distance(atom);

                    // if closer update index info
                    if(dist < smallest){
                        kindex = i;
                        smallest = dist;
                    }
                }

                // Add atom to the closest kcluster
                clusters_[kindex].add_atom(atom);
            }
            // Put atoms in correct order inside clusters so that the
            // Shell ordering becomes deterministic.
            foreach(auto& cluster, clusters_){
                cluster.sort_atoms();
            }
        }

        // Computes Lloyd's algorith to find a local minimium for the clusters.
        void k_means_search(std::size_t niter = 100){

            // Lloyd's algorithm iterations.
            for(auto i = 0; i < niter; ++i){

               // Recompute the center of the cluster using the centroid of
               // the atoms.  Will lose information about which atoms
               // go with which center.
               foreach(auto &cluster, clusters_){ cluster.guess_center(); }
               sort_clusters();

               attach_to_closest_cluster();
           }
        }

        /*
         * Returns a vector of shells in the order they appear in the clusters.
         */
        std::vector<Shell> cluster_shells() const {

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
                    // the shell contatiner.
                    for(auto i = 0; i < nshells_on_atom; ++i){
                        shells.push_back(basis_->operator()(atom_index, i));
                    }
                }
            }

            return shells;
        }

        /*
         * Returns a ShellRange that contains the shells included on each cluster.
         */
        ShellRange compute_shell_ranges() const {
            ShellRange range;
            range.reserve(nclusters_);

            // First range is easy
            range.push_back(0);

            // Loop over clusters
            for(auto i = 0; i < nclusters_; ++i){
                // Holds the number of shells on the cluster.
                std::size_t shells_in_cluster = 0;
                // Loop over atoms
                foreach(const auto atom, clusters_[i].atoms()){
                    // position of atom in molecule
                    std::size_t atom_index = atom.mol_index();
                    // Get number of shells on the atom.
                    shells_in_cluster += basis_->nshell_on_center(atom_index);
                }
                // Compute the Starting Shell of the next tile.
                range.push_back(range.at(i) + shells_in_cluster);
            }

            return  range;
        }

        /*
         * Sort clusters on distance from com.
         */
        void sort_clusters(){
            std::sort(clusters_.begin(), clusters_.end(),
                  [&](const KCluster &a, const KCluster &b){
                   return ((a.center()-com_).norm() < (b.center()-com_).norm());}
            );
        }

    private:
        std::size_t nclusters_ = 0;
        std::vector<Atom> atoms_;
        Vector3 com_;
        std::vector<KCluster> clusters_; // Note that KCluster contains fixed
                                         // Eigen objects, but these are not vectorizable
                                         // Thus we don't have to worry about Alignment issues.
        sc::Ref<sc::GaussianBasisSet> basis_;
    };

} // namespace TA
} // namespace mpqc

#endif /* CHEMISTRY_QC_BASIS_SHELLORDER_HPP */

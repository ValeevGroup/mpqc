//
// kcluster.hpp
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

#ifndef MPQC_BASIS_KCLUSTER_HPP
#define MPQC_BASIS_KCLUSTER_HPP

#include <chemistry/molecule/atom.h>
#include <chemistry/qc/basis/gaussbas.h>
#include <Eigen/Dense>
#include <vector>

namespace mpqc {
namespace basis {
    namespace cluster {
        // A wrapper around sc::Atom which also knows the atoms index in the
        // molecule.
        class ClusterAtom : public sc::Atom {
        public:

            // Takes the atom we want along with its molecular index.
            ClusterAtom(const sc::Atom &atom, std::size_t index) :
                sc::Atom(atom),
                mol_index_(index)
            {}

            // return the index of the atom in the molecule.
            std::size_t mol_index() const {return mol_index_;}

        private:
            // Don't use these
            ClusterAtom();

            std::size_t mol_index_;
        };
    }

    /**
     * class holds the information about the differnt clusters in k-means
     * tiling.
     */
    class KCluster {
    public:
        using Atom = cluster::ClusterAtom;
        using Vector3 = Eigen::Vector3d;

        /**
         * Constructor takes an Eigen::Vector3d which designates the center
         * of the cluster.
         */
        KCluster(const Vector3 &pos = Vector3(0,0,0)) : center_(pos)
        {}

        KCluster& operator=(const KCluster &rhs){
            center_ = rhs.center_;
            return *this;
        }

        /// Adds an atom to the cluster. Must know its index as well
        void add_atom(const sc::Atom &atom, std::size_t index){
            atoms_.push_back(Atom(atom, index));
        }

        /// Adds an atom to the cluster. Must know its index as well
        void add_atom(const Atom &atom){
            atoms_.push_back(atom);
        }

        /// Finds distance to any atom to the center of the cluster.
        double distance(const Atom &atom){
            Vector3 atom_vec(atom.xyz(0), atom.xyz(1), atom.xyz(2));
            // Computes the length of the vector to the atom from the center.
            return (atom_vec - center_).norm();
        }
        /// Returns the position vector of the center of the cluster.
        const Vector3& center() const { return center_; }


        /// Returns the centorid of the cluster.
        Vector3 centroid(){

            std::size_t n_atoms = natoms();
            Vector3 centroid(0,0,0);

            // Loop over all of the members of the cluster and total their
            // positions in each diminsion.
            for(auto i = 0; i < natoms(); ++i){
                centroid[0] += atoms_[i].xyz(0);
                centroid[1] += atoms_[i].xyz(1);
                centroid[2] += atoms_[i].xyz(2);
            }

            // Get the average position in each dimension.
            centroid[0] = centroid[0]/natoms();
            centroid[1] = centroid[1]/natoms();
            centroid[2] = centroid[2]/natoms();

            return centroid;

        }

        /// Move the center to the centroid of the cluster and forget members.
        void guess_center(){
            center_ = centroid();
            atoms_.clear();
        }

        /// Return the index

        /// Returns the number of atoms in cluster.
        std::size_t natoms(){
            return atoms_.size();
        }

        const std::vector<Atom>&  atoms() const { return  atoms_;}

    private:
        Vector3 center_;
        std::vector<Atom> atoms_;
    }; // KCluster

} // namespace basis
} // namespace mpqc


#endif /* MPQC_BASIS_KCLUSTER_HPP */

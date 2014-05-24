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

#ifndef CHEMISTRY_QC_BASIS_KCLUSTER_HPP
#define CHEMISTRY_QC_BASIS_KCLUSTER_HPP

#include <chemistry/molecule/atom.h>
#include <chemistry/qc/basis/gaussbas.h>
#include <Eigen/Dense>
#include <vector>

namespace mpqc {
namespace TA{
    namespace cluster {
        // A wrapper around sc::Atom which also knows the atoms index in the
        // molecule.
        class ClusterAtom : public sc::Atom {
        public:

            // Takes the atom we want along with its molecular index.
            ClusterAtom(const sc::Atom &atom, std::size_t index) :
                sc::Atom(atom),
                mol_index_(index),
                center_(atom.r(0),atom.r(1),atom.r(2))
            {}

            // return the index of the atom in the molecule.
            std::size_t mol_index() const {return mol_index_;}
            Eigen::Vector3d center() const {return center_;}

        private:
            // Don't use the default constructor
            ClusterAtom();

            std::size_t mol_index_;
            Eigen::Vector3d center_;
        }; // class ClusterAtom
    } // namespace cluster

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
        KCluster(const Vector3 &pos = Vector3(0,0,0)) : center_(pos), atoms_()
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
            Vector3 atom_vec(atom.r(0), atom.r(1), atom.r(2));
            // Computes the length of the vector to the atom from the center.
            return (atom_vec - center_).norm();
        }
        /// Returns the position vector of the center of the cluster.
        const Vector3& center() const { return center_; }

        /// Allows the user to modify the center possibly for looking for
        /// better clusters
        void set_center(Vector3 pos){center_ = pos;}

        /// Returns the centorid of the cluster. Using the "center of charge".
        Vector3 centroid(){

            Vector3 centroid(0,0,0);
            std::size_t total_charge = 0;

            // Loop over each atom and sum the position times the charge.
            for(auto atom : atoms_){
              centroid[0] += atom.r(0) * atom.Z();
              centroid[1] += atom.r(1) * atom.Z();
              centroid[2] += atom.r(2) * atom.Z();
              total_charge += atom.Z();
            }

            // Get the average
            centroid = centroid/double(total_charge);

            // Loop over all of the members of the cluster and total their
            // positions in each diminsion.
            //for(auto i = 0; i < n_atoms; ++i){
            //    centroid[0] += atoms_[i].r(0) * atoms_.[i].Z();
            //    centroid[1] += atoms_[i].r(1);
            //    centroid[2] += atoms_[i].r(2);
            //}

            //// Get the average position in each dimension.
            //centroid[0] = centroid[0]/n_atoms;
            //centroid[1] = centroid[1]/n_atoms;
            //centroid[2] = centroid[2]/n_atoms;

            return centroid;

        }

        /// Sort the atoms by closest to given point
        void sort_atoms(Vector3 sort_point){
            std::stable_sort(atoms_.begin(), atoms_.end(), [=]
                             (const Atom &a, const Atom &b){
                return Vector3(sort_point - a.center()).norm() <
                       Vector3(sort_point - b.center()).norm();
              }
            );
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

} // namespace TA
} // namespace mpqc


#endif /* CHEMISTRY_QC_BASIS_KCLUSTER_HPP */

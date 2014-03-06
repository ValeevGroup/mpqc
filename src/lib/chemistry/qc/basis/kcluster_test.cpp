//
// kcluster_test.cpp
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

#include <chemistry/qc/basis/kcluster.hpp>

#define BOOST_TEST_MODULE test_kcluster
#include <boost/test/included/unit_test.hpp>
#include <vector>

using namespace mpqc;
using namespace TA;
using namespace sc;
using namespace boost::unit_test;

using Katom = KCluster::Atom;
using Vector3 = KCluster::Vector3;

// Test the cluster atoms
BOOST_AUTO_TEST_CASE( cluster_atom_test ){
    Ref<Molecule> mol = new Molecule;
    mol->add_atom(2, 0.0, 0.0, 3.7);
    mol->add_atom(10, 0.0, 0.0, 0.0);

    // Test cluster atom constructor no default constructor
    Katom atom(mol->atom(0), 0);
    BOOST_CHECK_EQUAL(atom.mol_index(), 0);
    BOOST_CHECK(sc::Atom(atom) == mol->atom(0));

    Katom atom2(mol->atom(1), 1);
    BOOST_CHECK_EQUAL(atom2.mol_index(), 1);
    BOOST_CHECK(sc::Atom(atom2) == mol->atom(1));

}

// Test KCluster Constructor
BOOST_AUTO_TEST_CASE(kcluster_consructor_test){
    // Test default constructor
    KCluster cluster;

    // Check that the default center is (0,0,0)
    Vector3 correct_center = Vector3::Zero(3);
    BOOST_CHECK_EQUAL_COLLECTIONS(
                    cluster.center().data(), cluster.center().data()+3,
                    correct_center.data(), correct_center.data()+3
                    );

    // Check the atom vector should be empty with default constuctor.
    BOOST_CHECK_EQUAL(cluster.natoms(), 0);

    // Check constructor with provided vector
    Vector3 correct_center2;
    correct_center2 << 1.0, 2.0, 3.0;
    KCluster cluster2(correct_center2);

    BOOST_CHECK_EQUAL_COLLECTIONS(
                    cluster2.center().data(), cluster2.center().data()+3,
                    correct_center2.data(), correct_center2.data()+3
                    );


}

BOOST_AUTO_TEST_CASE( kcluster_atoms_test ){
    Ref<Molecule> mol = new Molecule;
    mol->add_atom(2, 0.0, 0.0, 3.7);
    mol->add_atom(2, 0.0, 0.0, 0.0);
    Vector3 mol_com(mol->center_of_mass()[0],
                              mol->center_of_mass()[1],
                              mol->center_of_mass()[2]);
    KCluster cluster( mol_com );

    BOOST_CHECK_EQUAL(cluster.natoms(), 0);

    // Check adding sc::Atoms
    for(auto i = 0; i < mol->natom(); ++i){
        cluster.add_atom(mol->atom(i), i);
    }  // Don't have access to individual atoms
    BOOST_CHECK_EQUAL(cluster.natoms(), 2);

    // Check adding ClusterAtoms
    KCluster cluster2( mol_com );
    for(auto i = 0; i < cluster.natoms(); ++i){
        cluster2.add_atom(cluster.atoms()[i]); // Add ClusterAtoms from vector.
    }
    BOOST_CHECK_EQUAL(cluster2.natoms(), 2);

    // Check that atoms are the same as the ones in the molecule
    for(auto i = 0; i < mol->natom(); ++i){
        BOOST_CHECK(cluster.atoms()[i] == cluster2.atoms()[i]);
    }
    for(auto i = 0; i < mol->natom(); ++i){
        BOOST_CHECK(sc::Atom(cluster.atoms()[i]) == mol->atom(i));
    }
}

BOOST_AUTO_TEST_CASE( kcluster_function_test ){
    Ref<Molecule> mol = new Molecule;
    mol->add_atom(2, 0.0, 0.0, 3.7);
    mol->add_atom(2, 0.0, 0.0, 0.0);
    Vector3 mol_com(mol->center_of_mass()[0],
                              mol->center_of_mass()[1],
                              mol->center_of_mass()[2]);
    KCluster cluster( mol_com );
    for(auto i = 0; i < mol->natom(); ++i){
        cluster.add_atom(mol->atom(i), i);
    }

    std::vector<Katom> atoms = cluster.atoms();


    // Check distance from cluster center
    for(auto i = 0; i < mol->natom(); ++i){
        // Get the postition of the atom in an eigen vector
        const Vector3 atom_pos = Eigen::Map<const Vector3>(&(mol->atom(i).r()[0]));
        BOOST_CHECK(cluster.distance(atoms[i]) ==
                    (atom_pos - mol_com).norm()
                   );
    }

    // Center was check in the constuctor test.

    // Check centroid function
    Vector3 centroid(0,0,0);
    for(auto&  atom : atoms){
        centroid += Eigen::Map<Vector3>(&(atom.r()[0]));
    }
    centroid = centroid * (1.0/atoms.size());
    // Make sure the difference between the vectors is zero.
    BOOST_CHECK((centroid - cluster.centroid()).norm() == 0);

    // Finally test guess center
    cluster.guess_center();
    BOOST_CHECK((centroid - cluster.center()).norm() == 0);
    BOOST_CHECK_EQUAL(cluster.natoms(), 0);

}








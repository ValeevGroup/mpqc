#include "cluster_concept.h"

#include <vector>
#include <iostream>

struct Atom {
    Eigen::Vector3d r;

    Atom(Eigen::Vector3d xyz) : r(xyz) {}
};

Eigen::Vector3d center(Atom const &a){
    return a.r;
}

struct Molecule {
    std::vector<Atom> atoms;
    auto begin() const -> decltype(atoms.cbegin()) {
        return atoms.cbegin();
    }
    auto end() const -> decltype(atoms.cend()) {
        return atoms.cend();
    }
};

Eigen::Vector3d center(Molecule const &m){
    Eigen::Vector3d c;
    c.setZero();

    for(auto const &a: m.atoms){
        c[0] += a.r[0];
        c[1] += a.r[1];
        c[2] += a.r[2];
    }

    c /= double(m.atoms.size());

    return c;
}


int main(){
    std::vector<Atom> atoms;
    for(auto i = 0; i < 10; ++i){
        atoms.push_back(Atom({i, 0, 2}));
    }

    Molecule mol; mol.atoms = atoms;

    // Clusterables of atom and molecule
    mpqc::clustering::Clusterable<Atom> a1(atoms[0]);
    mpqc::clustering::Clusterable<Atom> m1(mol);

    // cluster of 2 clusterables 
    mpqc::clustering::Cluster<Atom> c1({a1, m1});
    
    // cluster of 2 clusterables one of which is a cluster
    mpqc::clustering::Cluster<Atom> c2({a1, c1});

    // clusterable of clusterable
    mpqc::clustering::Clusterable<Atom> cluster_ception(m1);

    //
    mpqc::clustering::Clusterable<Atom> cluster_super_ception(c2);
    mpqc::clustering::Cluster<Atom> c3({cluster_super_ception, m1, a1, c1});

    auto all_atoms = c3.flatten();

    for(auto const &a : all_atoms){
        std::cout << a.r.transpose() << std::endl;
    }

    return 0;
}

//
// Created by Chong Peng on 10/21/15.
//

#include <catch.hpp>

#include "../include/tiledarray.h"

#include "../integrals/atomic_integral.h"

#include "../basis/basis_set.h"

#include "../molecule/clustering_functions.h"

using namespace mpqc;


TA::TensorD op(TA::TensorD &&ten){
    return std::move(ten);
}


TEST_CASE("Atomic Integral", "[atomic_integral]"){

    int a = 0;
    char** b;
    auto &world = madness::initialize(a,b);

    std::string xyz_file = "h2o.xyz";
    std::string basis_name = "cc-pVDZ";

    auto clustered_mol = molecule::attach_hydrogens_and_kmeans(molecule::read_xyz(xyz_file).clusterables(),1);

    basis::BasisSet bs(basis_name);
    basis::Basis basis(bs.get_cluster_shells(clustered_mol));

// Atomic Integral
    integrals::AtomicIntegral<TA::TensorD, TA::SparsePolicy> ao_int(world,
                                                                   op,
        std::make_shared<molecule::Molecule>(clustered_mol),
        std::make_shared<basis::Basis>(basis));

//        auto overlap = ao_int.compute_one_electron(L"<κ|λ>");


}

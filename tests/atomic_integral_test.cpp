//
// Created by Chong Peng on 10/21/15.
//

#include <catch.hpp>

#include "../include/tiledarray.h"
#include "../include/libint.h"


#include "../integrals/atomic_integral.h"

#include "../basis/atom_basisset.h"
#include "../basis/basis_set.h"
#include "../basis/basis.h"

#include "../common/namespaces.h"
#include "../molecule/clustering_functions.h"
#include "../molecule/molecule.h"

using namespace mpqc;

TEST_CASE("Atomic Integral", "[atomic_integral]"){

    std::string xyz_file = "h2o.xyz";
    std::string basis_name = "cc-pVDZ";

    auto clustered_mol = molecule::attach_hydrogens_and_kmeans(molecule::read_xyz(xyz_file).clusterables(),1);

    basis::BasisSet bs(basis_name);
    basis::Basis basis(bs.get_cluster_shells(clustered_mol));

    SECTION("Constructor"){
        integrals::AtomicIntegral<TA::Tensor<double>,TA::DensePolicy> ai(std::make_shared<basis::Basis>(basis));
    }


}

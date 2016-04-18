#include "../molecule/atom_based_cluster.h"

#include "../molecule/molecule.h"

#include "../basis/basis_set.h"
#include "../basis/basis.h"

using namespace mpqc;

int main(int argc, char **argv) {
    std::string basis_name = "sto-3g";
    if(argc > 1){
        basis_name = argv[1];
    }

    mol::AtomBasedCluster h20f(mol::Atom({0, -1, 1}, 1, 1), mol::Atom({0, 1, 1}, 1, 1),
                          mol::Atom({0, 0, 0}, 16, 8));

    std::vector<mol::AtomBasedClusterable> clusters{std::move(h20f)};

    mol::Molecule mol(std::move(clusters));

    basis::BasisSet obs_set(basis_name);
    
    // Print flat shells in basis 
    std::cout << "Flat shells in basis\n";
    for(auto const& sh : obs_set.get_flat_shells(mol)){
        std::cout << sh << "\n";
    }

    // Print clustered shells in basis 
    auto counter = 1;
    std::cout << "Clustered shells in basis\n";
    for(auto const& sh_vec : obs_set.get_cluster_shells(mol)){
        std::cout << "Cluster " << counter++ << std::endl;
        for(auto const& sh : sh_vec){
            std::cout << "\t" << sh << std::endl;
        }
    }

    basis::Basis obs(obs_set.get_cluster_shells(mol));
    std::cout << "obs = \n" << obs << std::endl;

    return 0;
}

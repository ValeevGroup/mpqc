#include <memory>
#include <fstream>
#include <algorithm>
#include <iomanip>
#include <chrono>
#include <rapidjson/document.h>
#include <rapidjson/stringbuffer.h>
#include <rapidjson/writer.h>

#include "../include/libint.h"
#include "../include/tiledarray.h"
#include "../include/btas.h"

#include "../utility/make_array.h"
#include "../utility/parallel_print.h"
#include "../utility/parallel_break_point.h"
#include "../utility/array_storage.h"
#include "../utility/time.h"
#include "../utility/json_input.h"

#include "../molecule/atom.h"
#include "../molecule/cluster.h"
#include "../molecule/molecule.h"
#include "../molecule/clustering_functions.h"
#include "../molecule/make_clusters.h"

#include "../basis/atom_basisset.h"
#include "../basis/basis_set.h"
#include "../basis/cluster_shells.h"
#include "../basis/basis.h"

#include "../integrals/integrals.h"

#include "../scf/diagonalize_for_coffs.hpp"
#include "../cc/ccsd_t.h"
//#include "../cc/integral_generator.h"
//#include "../cc/lazy_integral.h"
#include "../cc/ccsd_intermediates.h"
#include "../cc/trange1_engine.h"
#include "../ta_routines/array_to_eigen.h"
#include "../basis/shell_vec_functions.h"

using namespace mpqc;
namespace ints = integrals;


int try_main(int argc, char *argv[], madness::World &world) {


    // parse the input
    rapidjson::Document in;
    parse_input(argc, argv, in);

    std::cout << std::setprecision(15);
    Document cc_in;
    if (in.HasMember("CCSD")){
        cc_in = get_nested(in,"CCSD");
    }
    else if(in.HasMember("CCSD(T)")){
        cc_in = get_nested(in, "CCSD(T)");
    }

    if (!in.HasMember("xyz file") || !in.HasMember("number of bs clusters")
        || !in.HasMember("number of dfbs clusters")
        || !cc_in.HasMember("BlockSize")) {
        if (world.rank() == 0) {
            std::cout << "At a minimum your input file must provide\n";
            std::cout << "\"xyz file\", which is path to an xyz input\n";
            std::cout << "\"number of bs clusters\", which is the number of "
                         "clusters in the obs\n";
            std::cout << "\"number of dfbs clusters\", which is the number of "
                         "clusters in the dfbs\n";
            std::cout << "\"mo block size\", which is the block size for MO "
                         "orbitals\n";
        }
    }


    {

        // Get necessary info
        std::string mol_file = in["xyz file"].GetString();
        int nclusters = in["number of clusters"].GetInt();
        std::size_t blocksize = cc_in["BlockSize"].GetInt();


        // Get basis info
        std::string basis_name = in.HasMember("basis") ? in["basis"].GetString()
                                                       : "cc-pvdz";
        std::string df_basis_name = in.HasMember("df basis")
                                    ? in["df basis"].GetString()
                                    : "cc-pvdz-ri";

        // Get thresh info
        auto threshold = in.HasMember("block sparse threshold")
                         ? in["block sparse threshold"].GetDouble()
                         : 1e-13;

        // get other info
        bool frozen_core = cc_in.HasMember("FrozenCore")
                           ? cc_in["FrozenCore"].GetBool()
                           : false;

        if (world.rank() == 0) {
            std::cout << "Mol file is " << mol_file << std::endl;
            std::cout << "basis is " << basis_name << std::endl;
            std::cout << "df basis is " << df_basis_name << std::endl;
            std::cout << "Using " << nclusters << " clusters"
            << std::endl;
        }

        TiledArray::SparseShape<float>::threshold(threshold);
        tcc::utility::print_par(world, "Sparse threshold is ",
                                TiledArray::SparseShape<float>::threshold(), "\n");

        auto mol = mpqc::molecule::read_xyz(mol_file);
        auto charge = 0;
        auto occ = mol.occupation(charge);
        auto repulsion_energy = mol.nuclear_repulsion();
        auto core_electron = mol.core_electrons();

        tcc::utility::print_par(world, "Nuclear repulsion_energy = ",
                                repulsion_energy, "\n");


        world.gop.fence();

        auto clustered_mol = mpqc::molecule::kmeans(mpqc::molecule::read_xyz(mol_file).clusterables(), nclusters);

        mpqc::basis::BasisSet bs{basis_name};
        mpqc::basis::BasisSet df_bs{df_basis_name};

        std::streambuf *cout_sbuf
                = std::cout.rdbuf(); // Silence libint printing.
        std::ofstream fout("/dev/null");
        std::cout.rdbuf(fout.rdbuf());
        mpqc::basis::Basis basis{bs.get_cluster_shells(clustered_mol)};
        mpqc::basis::Basis df_basis{df_bs.get_cluster_shells(clustered_mol)};
        std::cout.rdbuf(cout_sbuf);

        if (world.rank() == 0) {
            std::cout << "Initial Blocking" << std::endl;
            std::cout << "Basis trange " << std::endl;
            TA::TiledRange1 bs_range = basis.create_trange1();
            std::cout << bs_range << std::endl;
            auto iter = bs_range.begin();
            for (; iter != bs_range.end() - 1; ++iter){
                std::cout << iter->first << " ";
            }
            std::cout << iter->first << " ";
            std::cout << iter->second << std::endl;
            auto minmax_block = cc::minmax_blocksize(bs_range);
            std::cout << minmax_block.first << " " << minmax_block.second << std::endl;
            auto average_block = cc::average_blocksize(bs_range);
            std::cout << "Average: " << average_block << std::endl;
            TA::TiledRange1 dfbs_range = df_basis.create_trange1();
            std::cout << "DF Basis trange " << std::endl;
            std::cout << dfbs_range << std::endl;
            iter = dfbs_range.begin();
            for (; iter != dfbs_range.end() - 1; ++iter){
                std::cout << iter->first << " ";
            }
            std::cout << iter->first << " ";
            std::cout << iter->second << std::endl;
            minmax_block = cc::minmax_blocksize(dfbs_range);
            std::cout << minmax_block.first << " " << minmax_block.second << std::endl;
            average_block = cc::average_blocksize(dfbs_range);
            std::cout << "Average: " << average_block << std::endl;
        }

        std::cout << "Reblock basis, block size set to " << blocksize << std::endl;

        auto new_basis = reblock(basis,cc::reblock_basis,blocksize);
        auto new_range = new_basis.create_trange1();
        std::cout << "Basis trange" << std::endl;
        std::cout << new_range << std::endl;
        auto iter = new_range.begin();
        for (; iter != new_range.end() - 1; ++iter){
            std::cout << iter->first << " ";
        }
        std::cout << iter->first << " ";
        std::cout << iter->second << std::endl;
        auto minmax_block = cc::minmax_blocksize(new_range);
        std::cout << minmax_block.first << " " << minmax_block.second << std::endl;
        auto average_block = cc::average_blocksize(new_range);
        std::cout << "Average: " << average_block << std::endl;

        auto new_dfbasis = reblock(df_basis,cc::reblock_basis,blocksize);
        auto new_dfrange = new_dfbasis.create_trange1();
        std::cout << "DF Basis trange" << std::endl;
        std::cout << new_dfrange << std::endl;
        iter = new_dfrange.begin();
        for (; iter != new_dfrange.end() - 1; ++iter){
            std::cout << iter->first << " ";
        }
        std::cout << iter->first << " ";
        std::cout << iter->second << std::endl;
        minmax_block = cc::minmax_blocksize(new_dfrange);
        std::cout << minmax_block.first << " " << minmax_block.second << std::endl;
        average_block = cc::average_blocksize(new_dfrange);
        std::cout << "Average: " << average_block << std::endl;
//        // start SCF
//        libint2::init();
//
//        const auto bs_array = tcc::utility::make_array(basis, basis);
//
//        // Overlap ints
//        auto overlap_e = ints::make_1body_shr_pool("overlap", basis, clustered_mol);
//        auto S = ints::sparse_integrals(world, overlap_e, bs_array);
//
//        // Overlap ints
//        auto kinetic_e = ints::make_1body_shr_pool("kinetic", basis, clustered_mol);
//        auto T = ints::sparse_integrals(world, kinetic_e, bs_array);
//
//        auto nuclear_e = ints::make_1body_shr_pool("nuclear", basis, clustered_mol);
//        auto V = ints::sparse_integrals(world, nuclear_e, bs_array);
//
//        decltype(T) H;
//        H("i,j") = T("i,j") + V("i,j");
    }

    return 0;
}


int main(int argc, char *argv[]) {

    int rc = 0;

    auto &world = madness::initialize(argc, argv);

    try {

        try_main(argc, argv, world);

    } catch (TiledArray::Exception &e) {
        std::cerr << "!! TiledArray exception: " << e.what() << "\n";
        rc = 1;
    } catch (madness::MadnessException &e) {
        std::cerr << "!! MADNESS exception: " << e.what() << "\n";
        rc = 1;
    } catch (SafeMPI::Exception &e) {
        std::cerr << "!! SafeMPI exception: " << e.what() << "\n";
        rc = 1;
    } catch (std::exception &e) {
        std::cerr << "!! std exception: " << e.what() << "\n";
        rc = 1;
    } catch (...) {
        std::cerr << "!! exception: unknown exception\n";
        rc = 1;
    }


    madness::finalize();
    return rc;
}

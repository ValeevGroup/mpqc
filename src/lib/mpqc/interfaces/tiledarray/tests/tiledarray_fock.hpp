/*
 * tiledarray_ints.hpp
 *
 *  Created on: Aug 14, 2013
 *      Author: drewlewis
 */

#ifndef mpqc_interfaces_tiledarray_tests_tiled_array_fock_hpp
#define mpqc_interfaces_tiledarray_tests_tiled_array_fock_hpp

#include "common.hpp"
#include <mpqc/interfaces/tiledarray/symmscmat.hpp>
#include <chemistry/qc/lcao/soad.h>
#include <Eigen/Dense>

namespace mpqc {
namespace tests {

    TA::Array<double, 2>
    get_overlap(madness::World &world,
                const R<Molecule> &mol,
                const R<Basis> &basis,
                const R<Integral> &int_fac){
        int_fac->set_basis(basis);
        IntPool<R<sc::OneBodyInt> > overlap_pool(int_fac->overlap());
        TA::Array<double, 2> S = Integrals(world, overlap_pool,
                                          tiling::tile_by_atom);
        world.gop.fence();
        return S;
    }

    TA::Array<double, 2>
    get_hcore(madness::World &world,
                const R<Molecule> &mol,
                const R<Basis> &basis,
                const R<Integral> &int_fac){
        int_fac->set_basis(basis);
        IntPool<R<sc::OneBodyInt> > hcore_pool(int_fac->hcore());
        TA::Array<double, 2> H = Integrals(world, hcore_pool,
                                          tiling::tile_by_atom);
        world.gop.fence();
        return H;
    }

    TA::Array<double, 4>
    get_eri(madness::World &world,
                const R<Molecule> &mol,
                const R<Basis> &basis,
                const R<Integral> &int_fac){
        int_fac->set_basis(basis);
        IntPool<R<sc::TwoBodyInt> > eri_pool(int_fac->electron_repulsion());
        TA::Array<double, 4> Eri = Integrals(world, eri_pool,
                                          tiling::tile_by_atom);
        world.gop.fence();
        return Eri;
    }

    TA::Array<double, 2>
    get_eri2(madness::World &world,
                const R<Molecule> &mol,
                const R<Basis> &basis_df,
                const R<Integral> &int_fac){
        int_fac->set_basis(basis_df, basis_df);
        using Eri2Int = sc::TwoBodyTwoCenterInt;
        IntPool<R<Eri2Int> > eri2_pool(int_fac->electron_repulsion2());
        TA::Array<double, 2> Eri2 = Integrals(world, eri2_pool,
                                          tiling::tile_by_atom);
        world.gop.fence();
        return Eri2;
    }

    TA::Array<double, 3>
    get_eri3(madness::World &world,
                const R<Molecule> &mol,
                const R<Basis> &basis,
                const R<Basis> &basis_df,
                const R<Integral> &int_fac){
        int_fac->set_basis(basis, basis, basis_df);
        using Eri3Int = sc::TwoBodyThreeCenterInt;
        IntPool<R<Eri3Int> > eri3_pool(int_fac->electron_repulsion3());
        TA::Array<double, 3> Eri3 = Integrals(world, eri3_pool,
                                          tiling::tile_by_atom);
        world.gop.fence();
        return Eri3;
    }

    TA::Array<double, 2>
    get_soad_guess(madness::World &world,
             const R<Molecule> &mol,
             const R<Basis> &basis){
        using Soad = sc::SuperpositionOfAtomicDensities;

        R<AKeyVal> akv = new AKeyVal;
        akv->assign("molecule", mol.pointer());
        akv->assign("basis", basis.pointer());

        R<Soad> s_guess = new Soad(R<KeyVal>(akv));
        R<sc::SymmSCMatrix>  P = s_guess->ao_density();

        TA::Array<double, 2> D = SymmScMat_To_TiledArray(world, P,
                                                 tiling::tile_by_atom(basis));

        D("i,j") = 0.5 * D("i,j");

        world.gop.fence();
        return D;
    }

    TA::Array<double, 2>
    get_inverse(const TA::Array<double, 2> &array){
        using namespace Eigen;
        MatrixXd Mat = TA::array_to_eigen(array);
        double a = 1.0/(Mat.lpNorm<1>() * Mat.lpNorm<Infinity>());

        TA::Array<double, 2> In = a * array("i,j");

        for(auto i = 0; i < 100; ++i){
            In("i,j") = 2.0 * In("i,j") - In("i,c") * array("c,d") * In("d,j");
        }

        array.world().gop.fence();
        return In;
    }

} // namespace tests
} // namesapce mpqc




#endif /* mpqc_interfaces_tiledarray_tests_tiled_array_fock_hpp */

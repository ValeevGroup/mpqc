/*
 * IntegralEnginePool_Test.cc
 *
 *  Created on: Aug 1, 2013
 *      Author: drewlewis
 */

#include <mpqc/interfaces/tiledarray/array_ints.hpp>
#include <mpqc/integrals/integralenginepool.h>
#include <chemistry/molecule/molecule.h>
#include <chemistry/qc/lcao/soad.h>
#include <chemistry/qc/basis/integral.h>
#include <chemistry/qc/basis/basis.h>
#include <chemistry/qc/basis/split.h>
#include <mpqc/interfaces/tiledarray/symmscmat.hpp>
#include <util/keyval/keyval.h>
#include "common.hpp"
#include <Eigen/Dense>
#include <iostream>
using namespace mpqc;
using namespace mpqc::tests;

int try_main(int argc, char** argv) {

    madness::World& world = madness::initialize(argc, argv);

    // Make Molecule
    if(world.rank() == 0)
        std::cout << "Making Molecule(H2) . . . " << std::endl;
    sc::Ref<sc::Molecule> mol = get_molecule("H2");

    // Make Basis
    if(world.rank() == 0)
        std::cout << "Making basis(STO-3G) . . . " << std::endl;
    sc::Ref<sc::GaussianBasisSet> basis;
    sc::Ref<sc::AssignedKeyVal> akv = new sc::AssignedKeyVal;
    akv->assign("molecule", mol.pointer());
    akv->assign("name", "3-21G");
    basis = new sc::GaussianBasisSet(sc::Ref<sc::KeyVal>(akv));
    if(basis->max_ncontraction() > 1){
        sc::Ref<sc::GaussianBasisSet> split_basis = new sc::SplitBasisSet(basis, basis->name());
        basis = split_basis;
    }
    akv->assign("basis", basis.pointer());


    sc::Ref<sc::GaussianBasisSet> basis_df;
    akv->assign("name", "cc-pVDZ-RI");
    basis_df = new sc::GaussianBasisSet(sc::Ref<sc::KeyVal>(akv));
    if(basis_df->max_ncontraction() > 1){
        sc::Ref<sc::GaussianBasisSet> split_basis_df = new sc::SplitBasisSet(basis_df, basis_df->name());
        basis_df = split_basis_df;
    }
    akv->assign("basis", basis_df.pointer());





    //  Int factory and Engine Pool
    if(world.rank() == 0)
        std::cout << "Making Integral reference . . ." << std::endl;
    sc::Ref<sc::Integral> int_fac = sc::Integral::initial_integral(argc,argv);
    if(int_fac.nonnull())
        sc::Integral::set_default_integral(int_fac);
    int_fac = sc::Integral::get_default_integral()->clone();

    if(world.rank() == 0)
        std::cout << "Setting Integral basis . . ." << std::endl;
    int_fac->set_basis(basis);

    if(world.rank() == 0)
        std::cout << "Making Engine Pool object . . ." << std::endl;

    mpqc::IntegralEnginePool<sc::Ref<sc::OneBodyInt> > hpool(int_fac->hcore());
    mpqc::IntegralEnginePool<sc::Ref<sc::OneBodyInt> > S(int_fac->overlap());
    mpqc::IntegralEnginePool<sc::Ref<sc::TwoBodyInt> > eripool(int_fac->electron_repulsion());

    int_fac->set_basis(basis,basis,basis_df);
    mpqc::IntegralEnginePool<sc::Ref<sc::TwoBodyThreeCenterInt> > eri3pool(int_fac->electron_repulsion3());
    int_fac->set_basis(basis_df,basis_df);
    mpqc::IntegralEnginePool<sc::Ref<sc::TwoBodyTwoCenterInt> > eri2pool(int_fac->electron_repulsion2());


    if(world.rank() == 0)
        std::cout << "Creating One and TwoBody Ints . . . " << std::endl;
    TA::Array<double, 2> S_AO = mpqc::Integrals(world,S, mpqc::tiling::tile_by_atom);
    std::cout << "S_AO = \n" << S_AO << std::endl;
    TA::Array<double, 2> Hcore = mpqc::Integrals(world,
                                           hpool, &mpqc::tiling::tile_by_atom);
    TA::Array<double, 4> Eri = mpqc::Integrals(world, eripool,
                                               &mpqc::tiling::tile_by_atom);
    TA::Array<double, 3> Eri3 = mpqc::Integrals(world, eri3pool,
                                                mpqc::tiling::tile_by_atom);
    TA::Array<double, 2> Eri2 = mpqc::Integrals(world, eri2pool,
                                                mpqc::tiling::tile_by_atom);

    sc::Ref<sc::SuperpositionOfAtomicDensities> soad_guess =
            new sc::SuperpositionOfAtomicDensities(sc::Ref<sc::KeyVal> (akv));
    sc::Ref<sc::SymmSCMatrix> P_guess = soad_guess->ao_density();

    TA::Array<double, 2> D_AO = mpqc::SymmScMat_To_TiledArray(world,
                                         P_guess,
                                         mpqc::tiling::tile_by_atom(basis));
    D_AO("i,j") = D_AO("i,j") * 0.5;

    TA::Array<double, 2> Inv = 0.000001 * Eri2("i,j");
    for(auto i = 0; i < 1000; ++i){
        Inv("i,j") = 2 * Inv("i,j") - Inv("i,c") * Eri2("c,d") * Inv("d,j");
    }
    //std::cout << "Inv = \n" << Inv  << std::endl;


    TA::Array<double, 3> Qpq = Eri3("p,q,Y") * Inv("Y,X");
    TA::Array<double, 2> G_DF = 2.0 *
                     ( Qpq("i,j,X") * ( D_AO("n,m") * Qpq("n,m,X") ) ) -
                     ( Qpq("i,n,X") * (  D_AO("n,m") * Qpq("m,j,X") ) );

    TA::Array<double, 4> ijkl = Qpq("i,j,X") * Qpq("k,l,X");
    std::cout << "Eri = \n" << Eri << std::endl;
    TA::Array<double, 4> diff = Eri("i,j,k,l") - ijkl("i,j,k,l");
    std::cout << "DIff = " << TA::expressions::norminf(diff("i,j,k,l")) << std::endl;

    if(world.rank() == 0){
        std::cout << "Computing intial guess at energy . . . " << std::endl;
    }

    TA::Array<double,2> G = 2 * D_AO("m,n") * Eri("i,j,m,n") -
                                D_AO("m,n") * Eri("i,n,m,j");
    std::cout << "G_AO = \n" << G << std::endl;

    double fake_energy = TA::expressions::dot(
                    2.0 * Hcore("i,j") + G("i,j"),
                    D_AO("i,j"));

    if(world.rank() == 0)
        std::cout << "FAKE ENERGY = " << fake_energy << std::endl;

    world.gop.fence();

    if(world.rank() == 0)
        std::cout << "Test Newton Iteration Method on Overlap" << std::endl;

    Eigen::MatrixXd Seig = TA::array_to_eigen(Eri2);
    double a_1 = 1/(Seig.lpNorm<1>() * Seig.lpNorm<Eigen::Infinity>());
    double norm1 = Seig.lpNorm<1>();
    double normx = Seig.lpNorm<Eigen::Infinity>();
    std::cout << "Norm 1 and Norm infinity = " << norm1 << " " << normx << std::endl;
    std::cout << "Scale = " << a_1  << std::endl;

    return 0;
}

int main(int argc, char** argv){


    try{
        try_main(argc, argv);
    }
    catch(std::exception &e){
        std::cerr << "Exception thrown = " << e.what() << std::endl;
    }
    catch(...){
        std::cerr << "Unknown exception" << std::endl;
    }

    return 0;
}


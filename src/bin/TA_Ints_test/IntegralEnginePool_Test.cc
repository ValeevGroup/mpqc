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
#include <iostream>

int try_main(int argc, char** argv) {

    madness::World& world = madness::initialize(argc, argv);

    // Make Molecule
    if(world.rank() == 0)
        std::cout << "Making Molecule(H2) . . . " << std::endl;
    sc::Ref<sc::Molecule> mol = new sc::Molecule;
    mol->add_atom(8, 0,0,0);
    mol->add_atom(1, 0,1,-1);
    mol->add_atom(1, 0,1,1);

    // Make Basis
    if(world.rank() == 0)
        std::cout << "Making basis(STO-3G) . . . " << std::endl;
    sc::Ref<sc::GaussianBasisSet> basis;
    sc::Ref<sc::AssignedKeyVal> akv = new sc::AssignedKeyVal;
    akv->assign("molecule", mol.pointer());
    akv->assign("name", "STO-3G");
    basis = new sc::GaussianBasisSet(sc::Ref<sc::KeyVal>(akv));
    if(basis->max_ncontraction() > 1){
        sc::Ref<sc::GaussianBasisSet> split_basis = new sc::SplitBasisSet(basis, basis->name());
        basis = split_basis;
    }
    akv->assign("basis", basis.pointer());


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
    mpqc::IntegralEnginePool<sc::Ref<sc::TwoBodyInt> > eripool(
                                               int_fac->electron_repulsion());


    if(world.rank() == 0)
        std::cout << "Creating One and TwoBody Ints . . . " << std::endl;
    TA::Array<double, 2> Hcore = mpqc::Integrals(world,
                                           hpool, &mpqc::tiling::tile_by_atom);
    TA::Array<double, 4> Eri = mpqc::Integrals(world, eripool,
                                               &mpqc::tiling::tile_by_atom);

    sc::Ref<sc::SuperpositionOfAtomicDensities> soad_guess =
            new sc::SuperpositionOfAtomicDensities(sc::Ref<sc::KeyVal> (akv));
    sc::Ref<sc::SymmSCMatrix> P_guess = soad_guess->ao_density();

    TA::Array<double, 2> D_AO_guess = mpqc::SymmScMat_To_TiledArray(world,
                                         P_guess,
                                         mpqc::tiling::tile_by_atom(basis));
    D_AO_guess("i,j") = D_AO_guess("i,j") * 0.5;

    if(world.rank() == 0){
        std::cout << "Computing intial guess at energy . . . " << std::endl;
    }

    TA::Array<double,2> G = D_AO_guess("n,m") * (2 * Eri("i,j,m,n") - Eri("i,n,m,j") );
    double fake_energy = TA::expressions::dot(
                    2.0 * Hcore("i,j") + G("i,j"),
                    D_AO_guess("i,j"));

    if(world.rank() == 0)
        std::cout << "FAKE ENERGY = " << fake_energy << std::endl;

    world.gop.fence();
    world.gop.fence();
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


/*
 * common.hpp
 *
 *  Created on: Aug 14, 2013
 *      Author: drewlewis
 */

#ifndef mpqc_interfaces_tiledarray_tests_common_hpp
#define mpqc_interfaces_tiledarray_tests_common_hpp

#include <chemistry/molecule/molecule.h>
#include <chemistry/qc/basis/basis.h>
#include <chemistry/qc/basis/split.h>
#include <chemistry/qc/basis/integral.h>
#include <mpqc/integrals/integrals.hpp>
#include <mpqc/integrals/integralenginepool.hpp>
#include <tiledarray.h>
#include <mpqc/interfaces/tiledarray/array_ints.hpp>
#include <string>
#include <iostream>

namespace mpqc {
namespace tests {
    // Typedefs and template alaises
    template<typename T> using R = sc::Ref<T>;
    using Molecule = sc::Molecule;
    using Basis = sc::GaussianBasisSet;
    using string = std::string;
    using AKeyVal = sc::AssignedKeyVal;
    using KeyVal = sc::KeyVal;
    using Integral = sc::Integral;
    template<typename T> using IntPool = IntegralEnginePool<T>;

namespace detail {
    // Takes a std::string and will return a ref molecule either water or h2.
    // You may provide your own if you would like.
    R<Molecule>
    get_molecule(const string &mol_name){
        R<Molecule> mol = new Molecule;

        if(mol_name == "H2"){
            mol->add_atom(1, 0,1,-1);
            mol->add_atom(1, 0,1,1);
        }

        else if(mol_name == "H2O"){
            mol->add_atom(8, 0, 0, 0);
            mol->add_atom(1, 0, 1, -1);
            mol->add_atom(1, 0, 1, 1);
        }

        else if (mol_name == "TRange1Test"){
            mol->add_atom( 6,     0,     0,    0);
            mol->add_atom( 9,    -1,    -1,    0);
            mol->add_atom( 1,   0.6,  -0.1,  0.9);
            mol->add_atom(17, -0.75,   1.5,    0);
            mol->add_atom(35,   1.1, -0.18, -1.5);
        }


        return mol;
    }


    // Returns a basis set for a molecule given a string with the basis set name
    // and the ref<molecule> in that the basis is desired for.
    R<Basis>
    get_basis(const string &basis_name, const R<Molecule> &mol){

       R<AKeyVal> akv = new AKeyVal;
       akv->assign("name", basis_name);
       akv->assign("molecule", mol.pointer());

       R<Basis> basis = new Basis(R<KeyVal>(akv));
       if(basis->max_ncontraction() > 1){
           R<Basis> split_basis = new sc::SplitBasisSet(basis);
           basis = split_basis;
       }

       return basis;

    }

}

    // Get molecule.  Figure out which molecule we are getting and then call
    // helper function.
    R<Molecule>
    get_molecule(const string &mol_name = ""){
        R<Molecule> mol;

        if(!mol_name.empty()){
            if(mol_name == "H2"){
                mol = detail::get_molecule(mol_name);
            }
            else if(mol_name == "H2O"){
                mol = detail::get_molecule(mol_name);
            }
            else if (mol_name == "TRange1Test"){
                mol = detail::get_molecule(mol_name);
            }
            else{
                std::cout << "Molecule name not reconized defaulting to H2"
                               << std::endl;
                mol = detail::get_molecule("H2");
            }
        }
        else {
            std::cout << "Molecule name is empty defaulting to H2" << std::endl;
            mol = detail::get_molecule("H2");
        }

        return mol;
    }

    // Get basis set.  If empty basis set get STO-3G
    R<Basis>
    get_basis(const string &basis_name, const R<Molecule> &mol){

       R<Basis> basis;

       if(!basis_name.empty()){
           basis = detail::get_basis(basis_name, mol);
       }
       else{
           std::cout << "Basis name empty defaulting to STO-3G" << std::endl;
           basis = detail::get_basis("STO-3G", mol);
       }

       return basis;

    }

    // Make an integral factory for use in calculations.  This is where
    // We get our engines from.
    R<Integral>
    get_integral_factory(int argc, char** argv){
       R<Integral> integral_fac = Integral::initial_integral(argc, argv);
       if(integral_fac)
           Integral::set_default_integral(integral_fac);
       integral_fac = Integral::get_default_integral()->clone();
       return integral_fac;
    }

} // namespace mpqc
} // namespace tests




#endif /* mpqc_interfaces_tiledarray_tests_common_hpp */

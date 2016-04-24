//
// Created by Chong Peng on 3/2/16.
//

#ifndef MPQC_ATOMIC_INTEGRAL_BASE_H
#define MPQC_ATOMIC_INTEGRAL_BASE_H

#include <string>
#include <vector>
#include <iostream>
#include <cwchar>

#include "../include/tiledarray.h"
#include "../common/namespaces.h"
#include "../molecule/molecule.h"
#include "../basis/basis.h"
#include "../expression/formula.h"
#include "integral_engine_pool.h"
#include "../utility/make_array.h"
#include "task_integrals.h"
#include "make_engine.h"
#include "../expression/orbital_registry.h"

#include <libint2/engine.h>

namespace mpqc {
namespace integrals {

/**
 *
 * \brief base class for AtomicIntegral
 *
 */

class AtomicIntegralBase {
public:

    AtomicIntegralBase() = default;

    /**
     * Constructor
     *
     * @param world reference to madness World
     * @param mol shared pointer to Molecule
     * @param obs shared pointer to OrbitalBasisRegistry
     * @param gtg_params  parameters used in computing f12 integrals
     */

    AtomicIntegralBase(madness::World &world,
                       const std::shared_ptr <molecule::Molecule>& mol,
                       const std::shared_ptr <OrbitalBasisRegistry>& obs,
                       const std::vector <std::pair<double, double>>& gtg_params = std::vector<std::pair<double,double>>()
                        )
    : world_(world), mol_(mol), orbital_basis_registry_(obs), gtg_params_(gtg_params)
    { }

    virtual ~AtomicIntegralBase() = default;

    /// return madness world
    madness::World &get_world() const {
        return world_;
    }

    /// set OrbitalBasisRegistry
    void set_orbital_basis_registry(const std::shared_ptr<OrbitalBasisRegistry>& obs){
        orbital_basis_registry_ = obs;
    }

    /**
     * Given Formula with rank = 4, return DensityFitting formula
     *
     * This function is also used in MolecularIntegral density fitting formula parsing
     *
     * @param formula that has string format (in1 in2 | oper | in3 in4 ) or <in1 in2 | oper | in3 in4 >
     * @return array of wstring with length 3
     *         - string0 ( dfbs | oper | in1 in2 )  or ( dfbs |oper | in1 in3 )
     *         - string1 ( dfbs | oper | dfbs )[inv]
     *         - string2 ( dfbs | oper | in3 in4 )  or ( dfbs | oper | in2 in4 )
     */
    std::array<std::wstring, 3> get_df_formula(const Formula &formula);

protected:

    /// parse operation and return one body engine
    libint2::OneBodyEngine get_one_body_engine(const Operation &operation, int64_t max_nprim, int64_t max_am);

    /// parse operation and  return two body engine kernel
    libint2::MultiplicativeSphericalTwoBodyKernel get_two_body_engine_kernel(const Operation &operation);


    /// parse one body formula and set engine_pool and basis array
    void parse_one_body(const Formula &formula, std::shared_ptr <EnginePool<libint2::OneBodyEngine>> &engine_pool,
                        Barray<2> &bases);

    /// parse two body two center formula and set two body kernel, basis array, max_nprim and max_am
    void parse_two_body_two_center(const Formula &formula, libint2::MultiplicativeSphericalTwoBodyKernel &kernel,
                                   Barray<2> &bases, int64_t &max_nprim, int64_t &max_am);

    /// parse two body three center formula and set two body kernel, basis array, max_nprim and max_am
    void parse_two_body_three_center(const Formula &formula, libint2::MultiplicativeSphericalTwoBodyKernel &kernel,
                                     Barray<3> &bases, int64_t &max_nprim, int64_t &max_am);

    /// parse two body four center formula and set two body kernel, basis array, max_nprim and max_am
    void parse_two_body_four_center(const Formula &formula, libint2::MultiplicativeSphericalTwoBodyKernel &kernel,
                                    Barray<4> &bases, int64_t &max_nprim, int64_t &max_am);

    /**
     *  Given formula with rank = 2 and J or K operation, return the G integral
     *
     *  @param Formula that has string format ( in1 | oper | in2 ), where oper is J or K operation
     *  @return Formula
     *          - Formula that has string format ( in1 in2 | G | obs obs ) for J
     *          - Formula that has string format ( in1 obs | G | in2 obs ) for K, KAlpha, KBeta
     */
    Formula get_jk_formula(const Formula& formula);

    /**
     * Given formula with rank = 2 and J or K operation, return the G integral with DensityFitting
     *
     * @param Formula that has string format ( in1 | oper | in2 ), where oper is J or K operation
     * @return result array of Formula with size 3
     *         - 3 Formula that has string format ( dfbs | G | in1 in2 ) ( dfbs | G | dfbs )[inv] ( dfbs | G | obs obs ) for J
     *         - 3 Formula that has string format ( dfbs | G | in1 obs ) ( dfbs | G | dfbs )[inv] ( dfbs | G | in2 obs ) for K, KAlpha, KBeta
     */
    std::array<Formula,3> get_jk_df_formula(const Formula& formula);

    /**
     * Given formula with rank = 2 and Fock operation, return 3 formula to compute it
     *
     * @param Formula that has string format ( in1 | oper | in2 ), where oper is Fock operation
     * @return array of Formula with size 3
     *         - 3 Formula that has string format ( in1 | H | in2 ) ( in1 | J | in2 ) ( in1 |K | in2 ) for Fock
     *         - 3 Formula that has string format ( in1 | H | in2 ) ( in1 | J | in2 ) ( in1 | KAlpha | in2 ) for FockAlpha
     *         - 3 Formula that has string format ( in1 | H | in2 ) ( in1 | J | in2 ) ( in1 | KBeta | in2 ) for FockBeta
     */
    std::array<Formula,3> get_fock_formula(const Formula& formula);

    /**
     * Given operation that is J or K operation, return the orbital index that maps to density
     * @param operation J or K operation
     *
     * @return result OrbitalIndex
     *         - m for J or K
     *         - m_α for KAlpha
     *         - m_β for KBeta
     */
    OrbitalIndex get_jk_orbital_space(const Operation& operation);

    /// given OrbitalIndex, find the correspoding basis
    std::shared_ptr <basis::Basis> index_to_basis(const OrbitalIndex &index) {
        return orbital_basis_registry_->retrieve(index);
    }

protected:

    madness::World &world_;
    std::shared_ptr <molecule::Molecule> mol_;
    std::shared_ptr<OrbitalBasisRegistry> orbital_basis_registry_;
    std::vector <std::pair<double, double>> gtg_params_;

};
} // end of namespace integral
} // end of namespace mpqc


#endif //MPQC_ATOMIC_INTEGRAL_BASE_H

//
// Created by Chong Peng on 3/2/16.
//

#ifndef TILECLUSTERCHEM_ATOMIC_INTEGRAL_BASE_H
#define TILECLUSTERCHEM_ATOMIC_INTEGRAL_BASE_H

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

namespace mpqc {
namespace integrals {


class AtomicIntegralBase {
public:

    AtomicIntegralBase() = default;


    AtomicIntegralBase(madness::World &world,
                       std::shared_ptr <molecule::Molecule> mol,
                       std::shared_ptr <basis::Basis> obs,
                       std::shared_ptr <basis::Basis> dfbs = nullptr,
                       std::shared_ptr <basis::Basis> auxbs = nullptr,
                       std::vector <std::pair<double, double>> gtg_params = std::vector<std::pair<double,double>>()
                        )
    : world_(world), mol_(mol), obs_(obs), dfbs_(dfbs), abs_(auxbs), gtg_params_(gtg_params)
    {
        if (auxbs != nullptr) {
            ribs_ = std::make_shared<basis::Basis>(std::move(obs->join(*auxbs)));
        }
        else {
            ribs_ = nullptr;
        }
    }

    virtual ~AtomicIntegralBase() = default;


    madness::World &get_world() const {
        return world_;
    }

    const std::shared_ptr <basis::Basis> get_obs() const {
        return obs_;
    }

    const std::shared_ptr <basis::Basis> get_dfbs() const {
        return dfbs_;
    }

    const std::shared_ptr <basis::Basis> get_abs() const {
        return abs_;
    }

    void set_dfbs(const std::shared_ptr <basis::Basis> &dfbs) {
        AtomicIntegralBase::dfbs_ = dfbs;
    }

    void set_abs(const std::shared_ptr <basis::Basis> &abs) {
        AtomicIntegralBase::abs_ = abs;
    }

    std::array<std::wstring, 3> get_df_formula(const Formula &formula);

protected:

    // get one body engine
    libint2::OneBodyEngine get_one_body_engine(const Operation &operation, int64_t max_nprim, int64_t max_am);

    // get two body engine kernel
    libint2::MultiplicativeSphericalTwoBodyKernel get_two_body_engine_kernel(const Operation &operation);


    // parse one body formula and set engine_pool and basis array
    void parse_one_body(const Formula &formula, std::shared_ptr <EnginePool<libint2::OneBodyEngine>> &engine_pool,
                        Barray<2> &bases);

    // parse two body two center formula and set two body kernel and basis array
    void parse_two_body_two_center(const Formula &formula, libint2::MultiplicativeSphericalTwoBodyKernel &kernel,
                                   Barray<2> &bases, int64_t &max_nprim, int64_t &max_am);

    // parse two body three center formula and set two body kernel and basis array
    void parse_two_body_three_center(const Formula &formula, libint2::MultiplicativeSphericalTwoBodyKernel &kernel,
                                     Barray<3> &bases, int64_t &max_nprim, int64_t &max_am);

    // parse two body four center formula and set two body kernel and basis array
    void parse_two_body_four_center(const Formula &formula, libint2::MultiplicativeSphericalTwoBodyKernel &kernel,
                                    Barray<4> &bases, int64_t &max_nprim, int64_t &max_am);


    Formula get_jk_formula(const Formula& formula);

    std::array<Formula,3> get_jk_df_formula(const Formula& formula);

    std::array<Formula,4> get_fock_formula(const Formula& formula);

    OrbitalIndex get_jk_orbital_space(const Operation& operation);

    std::shared_ptr <basis::Basis> index_to_basis(const OrbitalIndex &index) {
        if (index.index() == OrbitalIndex::Index::obs) {
            return obs_;
        }
        else if (index.index() == OrbitalIndex::Index::abs) {
            return abs_;
        }
        else if (index.index() == OrbitalIndex::Index::dfbs) {
            return dfbs_;
        }
        else if (index.index() == OrbitalIndex::Index::ribs) {
            return ribs_;
        }
        else {
            throw std::runtime_error("Wrong Index!");
        }
    }

protected:

    madness::World &world_;
    std::shared_ptr <molecule::Molecule> mol_;
    std::shared_ptr <basis::Basis> obs_;
    std::shared_ptr <basis::Basis> dfbs_;
    std::shared_ptr <basis::Basis> abs_;
    std::shared_ptr <basis::Basis> ribs_;
    std::vector <std::pair<double, double>> gtg_params_;

};
} // end of namespace integral
} // end of namespace mpqc


#endif //TILECLUSTERCHEM_ATOMIC_INTEGRAL_BASE_H

//
// Created by Chong Peng on 10/14/15.
//

#ifndef TILECLUSTERCHEM_ATOMIC_INTEGRAL_H
#define TILECLUSTERCHEM_ATOMIC_INTEGRAL_H

#include <string>



#include"../common/namespaces.h"
#include "../include/tiledarray.h"
#include "../basis/basis.h"
#include "../expression/formula.h"
#include "../utility/make_array.h"
#include "integral_engine_pool.h"
#include "task_integrals.h"
#include "../molecule/molecule.h"
#include "make_engine.h"

namespace mpqc{
namespace integrals{

    //TODO Atomic Integral class
    template<typename Tile, typename Policy>
    class AtomicIntegralBase {
    public:
        using TArray2 = TA::Array <double, 2, Tile, Policy>;
        using TArray3 = TA::Array <double, 3, Tile, Policy>;
        using TArray4 = TA::Array <double, 4, Tile, Policy>;


        AtomicIntegralBase(madness::World& world,
                       std::shared_ptr<molecule::Molecule> mol,
                       std::shared_ptr<basis::Basis> obs,
                       std::shared_ptr<basis::Basis> dfbs = nullptr,
                       std::shared_ptr<basis::Basis> auxbs = nullptr) :
               world_(world), mol_(mol), obs_(obs), dfbs_(dfbs), abs_(auxbs)  { }


        madness::World& get_world() const {
            return world_;
        }

        const std::shared_ptr<basis::Basis> get_obs() const {
            return obs_;
        }

        const std::shared_ptr<basis::Basis> get_dfbs() const {
            return dfbs_;
        }

        const std::shared_ptr<basis::Basis> get_abs() const {
            return abs_;
        }

        void set_dfbs(const std::shared_ptr<basis::Basis> &dfbs) {
            AtomicIntegralBase::dfbs_ = dfbs;
        }

        void set_abs(const std::shared_ptr<basis::Basis> &abs) {
            AtomicIntegralBase::abs_ = abs;
        }

        TArray2 compute_one_electron(const std::wstring& formula){
            Formula tmp(formula);
            return compute_one_electron(std::move(tmp));
        }
        TArray4 compute_two_electron(const std::wstring& formula){
            Formula tmp(formula);
            return compute_one_electron(std::move(tmp));
        }

        virtual TArray2 compute_one_electron(const Formula& formula){}
        virtual TArray4 compute_two_electron(const Formula& formula){}

    private:
        std::shared_ptr<basis::Basis> index_to_basis(const OrbitalIndex& index){
            if(index.index() == OrbitalIndex::Index::obs){
                return obs_;
            }
            else if(index.index() == OrbitalIndex::Index::abs){
                return abs_;
            }
            else if(index.index() == OrbitalIndex::Index::dfbs){
                return  dfbs_;
            }
            else{
                throw std::runtime_error("Wrong Index!");
            }
        }

    private:

        madness::World& world_;
        std::shared_ptr<molecule::Molecule> mol_;
        std::shared_ptr<basis::Basis> obs_;
        std::shared_ptr<basis::Basis> dfbs_;
        std::shared_ptr<basis::Basis> abs_;

    };


//    template <typename Tile, typename Policy>
//    typename AtomicIntegralBase<Tile,Policy>::TArray2 AtomicIntegralBase<Tile,Policy>::compute_one_electron(const Formula& formula) {
//        auto bra_indexs = formula.left_index();
//        auto ket_indexs = formula.right_index();
//
//        TA_ASSERT(bra_indexs.size() == 1);
//        TA_ASSERT(ket_indexs.size() == 1);
//
//        auto bra_index = bra_indexs[0];
//        auto ket_index = ket_indexs[0];
//
//        TA_ASSERT(bra_index.is_ao());
//        TA_ASSERT(ket_index.is_ao());
//
//        auto bra_basis = index_to_basis(bra_index);
//        auto ket_basis = index_to_basis(ket_index);
//
//        TA_ASSERT(bra_basis != nullptr);
//        TA_ASSERT(ket_basis != nullptr);
//
//        // convert operation to libint operator
//        auto operation = formula.operation();
//        libint2::OneBodyEngine::operator_type itype;
//        tcc::integrals::q_vector q;
//        if (operation == Formula::Operation::Overlap) {
//            itype = libint2::OneBodyEngine::overlap;
//        } else if (operation == Formula::Operation::Kinetic) {
//            itype = libint2::OneBodyEngine::kinetic;
//        } else if (operation == Formula::Operation::Nuclear) {
//            itype = libint2::OneBodyEngine::nuclear;
//            q = tcc::integrals::make_q(*mol_);
//        } else {
//            throw std::runtime_error("Invalid One Body Operation");
//        }
//
//        auto max_nprim = std::max(bra_basis->max_nprim(), ket_basis->max_nprim());
//        auto max_am = std::max(bra_basis->max_am(), ket_basis->max_am());
//
//        libint2::OneBodyEngine engine(itype, max_nprim, static_cast<int>(max_am),0);
//
//        if(itype == libint2::OneBodyEngine::nuclear){
//            engine.set_params(std::move(q));
//        }
//
//        auto engine_pool = tcc::integrals::make_pool(engine);
//        const auto bs_array = tcc::utility::make_array(*bra_basis, *ket_basis);
//
//        auto ta_pass_through = [](TA::TensorD &&ten){
//            return std::move(ten);
//        };
//
//        TArray2 result = mpqc::integrals::dense_integrals(world_,engine_pool,bs_array,ta_pass_through);
//
//        return result;
//    }

    //TODO R12 Integral
    //TODO geminal parametors

    // R12 Atomic Integral

//    template<typename Tile, typename Policy>
//    class R12AtomicIntegral : public AtomicIntegralBase<Tile,Policy>{
//    public:
//        typedef TA::Array <double, 2, Tile, Policy> TArray2;
//        typedef TA::Array <double, 3, Tile, Policy> TArray3;
//        typedef TA::Array <double, 4, Tile, Policy> TArray4;
//
//    private:
//         geminal parametors
//    };

    }
}


#endif //TILECLUSTERCHEM_ATOMIC_INTEGRAL_H

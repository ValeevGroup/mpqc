//
// Created by Chong Peng on 10/14/15.
//

#ifndef TILECLUSTERCHEM_ATOMIC_INTEGRAL_H
#define TILECLUSTERCHEM_ATOMIC_INTEGRAL_H

#include <string>
#include <vector>


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


    template<typename Tile, typename Policy>
    class AtomicIntegralBase {
    public:
        using TArray2 = TA::Array <double, 2, Tile, Policy>;
        using TArray3 = TA::Array <double, 3, Tile, Policy>;
        using TArray4 = TA::Array <double, 4, Tile, Policy>;

        AtomicIntegralBase() = default;

        AtomicIntegralBase(madness::World& world,
                       std::shared_ptr<molecule::Molecule> mol,
                       std::shared_ptr<basis::Basis> obs,
                       std::shared_ptr<basis::Basis> dfbs = nullptr,
                       std::shared_ptr<basis::Basis> auxbs = nullptr,
                       std::vector<std::pair<double,double>> gtg_params = std::vector<std::pair<double,double>>() ) :
               world_(world), mol_(mol), obs_(obs), dfbs_(dfbs), abs_(auxbs), gtg_params_(gtg_params)  { }

        virtual ~AtomicIntegralBase() = default;


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

        // compute integrals that has two dimension
        virtual TArray2 compute2(const std::wstring& formula_string) = 0;
        // compute integrals that has three dimension
        virtual TArray3 compute3(const std::wstring& formula_stirng) = 0;
        // compute integrals that has four dimension
        virtual TArray4 compute4(const std::wstring& formula_string) = 0;


    protected:
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

    protected:

        madness::World& world_;
        std::shared_ptr<molecule::Molecule> mol_;
        std::shared_ptr<basis::Basis> obs_;
        std::shared_ptr<basis::Basis> dfbs_;
        std::shared_ptr<basis::Basis> abs_;
        std::vector<std::pair<double,double>> gtg_params_;

    };

    /// Atomic Integral Class
    //// Op is a function
    /// Op will take TA::TensorD as argument and return Tile

    template<typename Tile, typename Policy>
    class AtomicIntegral : public AtomicIntegralBase<Tile,Policy>{

    public:
        using TArray2 = typename AtomicIntegralBase<Tile,Policy>::TArray2;
        using TArray3 = typename AtomicIntegralBase<Tile,Policy>::TArray3;
        using TArray4 = typename AtomicIntegralBase<Tile,Policy>::TArray4;
        typedef Tile (*Op)(TA::TensorD &&);

        AtomicIntegral() = default;

        AtomicIntegral(madness::World& world,
                Op op,
                std::shared_ptr<molecule::Molecule> mol,
                std::shared_ptr<basis::Basis> obs,
                std::shared_ptr<basis::Basis> dfbs = nullptr,
                std::shared_ptr<basis::Basis> auxbs = nullptr,
                std::vector<std::pair<double,double>> gtg_params = std::vector<std::pair<double,double>>()
        ) : AtomicIntegralBase<Tile,Policy>(world,mol,obs,dfbs,auxbs,gtg_params), op_(op){}


        virtual ~AtomicIntegral() = default;


        TArray2 compute2(const std::wstring& );
        TArray3 compute3(const std::wstring& );
        TArray4 compute4(const std::wstring& );

    private:
        // compute integral for sparse policy
        template <typename E, unsigned long N, typename U = Policy>
        TA::Array<double, N,Tile,typename std::enable_if<std::is_same<U,TA::SparsePolicy>::value, TA::SparsePolicy>::type> compute_integrals(
                madness::World& world, E const &engine, Barray<N> const &bases)
        {
            auto result = mpqc::integrals::sparse_integrals(world,engine,bases,op_);
            return result;
        }

        // compute integral for dense policy
        template <typename E, unsigned long N, typename U = Policy>
        TA::Array<double, N,Tile,typename std::enable_if<std::is_same<U,TA::DensePolicy>::value, TA::DensePolicy>::type> compute_integrals(
                madness::World& world, E const &engine, Barray<N> const &bases)
        {
            auto result = mpqc::integrals::dense_integrals(world,engine,bases,op_);
            return result;
        }

    private:
        Op op_;

    };


    // TODO working with chemical and physical notation
    template <typename Tile, typename Policy>
    typename AtomicIntegral<Tile,Policy>::TArray2 AtomicIntegral<Tile,Policy>::compute2(const std::wstring& formula_string) {

        Formula formula(formula_string);


        // use one body engine
        if(formula.operation().is_onebody()){

            auto bra_indexs = formula.left_index();
            auto ket_indexs = formula.right_index();

            TA_ASSERT(bra_indexs.size() == 1);
            TA_ASSERT(ket_indexs.size() == 1);

            auto bra_index = bra_indexs[0];
            auto ket_index = ket_indexs[0];

            TA_ASSERT(bra_index.is_ao());
            TA_ASSERT(ket_index.is_ao());

            auto bra_basis = this->index_to_basis(bra_index);
            auto ket_basis = this->index_to_basis(ket_index);

            TA_ASSERT(bra_basis != nullptr);
            TA_ASSERT(ket_basis != nullptr);

            auto max_nprim = std::max(bra_basis->max_nprim(), ket_basis->max_nprim());
            auto max_am = std::max(bra_basis->max_am(), ket_basis->max_am());
            auto bs_array = tcc::utility::make_array(*bra_basis, *ket_basis);

            // convert operation to libint operator
            auto operation = formula.operation();
            libint2::OneBodyEngine::operator_type itype;
            std::vector<std::pair<double, std::array<double, 3>>> q;
            if (operation.get_operation() == Operation::Operations::Overlap) {
                itype = libint2::OneBodyEngine::overlap;
            } else if (operation.get_operation() == Operation::Operations::Kinetic) {
                itype = libint2::OneBodyEngine::kinetic;
            } else if (operation.get_operation() == Operation::Operations::Nuclear) {
                itype = libint2::OneBodyEngine::nuclear;
                q = make_q(*(this->mol_));
            } else {
                throw std::runtime_error("Invalid One Body Operation");
            }

            libint2::OneBodyEngine engine(itype, max_nprim, static_cast<int>(max_am),0);


            if(itype == libint2::OneBodyEngine::nuclear){
                engine.set_params(std::move(q));
            }

            auto engine_pool = make_pool(engine);

            auto result = compute_integrals(this->world_,engine_pool,bs_array);
            return result;
        }
        // use two body engine
        else if(formula.operation().is_twobody()){
            TA_USER_ASSERT(formula.notation() == Formula::Notation::Chemical, "Two Body Two Center Integral Must Use Chemical Notation");

            auto bra_indexs = formula.left_index();
            auto ket_indexs = formula.right_index();

            TA_ASSERT(bra_indexs.size() == 0);
            TA_ASSERT(ket_indexs.size() == 2);

            auto ket_index0 = ket_indexs[0];
            auto ket_index1 = ket_indexs[1];

            TA_ASSERT(ket_index0.is_ao());
            TA_ASSERT(ket_index1.is_ao());

            auto ket_basis0 = this->index_to_basis(ket_index0);
            auto ket_basis1 = this->index_to_basis(ket_index1);

            TA_ASSERT(ket_basis0 != nullptr);
            TA_ASSERT(ket_basis1 != nullptr);

            auto max_nprim = std::max(ket_basis0->max_nprim(), ket_basis1->max_nprim());
            auto max_am = std::max(ket_basis0->max_am(), ket_basis1->max_am());
            auto bs_array = tcc::utility::make_array(*ket_basis0, *ket_basis1);
            auto operation = formula.operation();
            if (operation.get_operation() == Operation::Operations::Coulomb) {
                libint2::TwoBodyEngine<libint2::Coulomb> engine(max_nprim, static_cast<int>(max_am));
                auto engine_pool = make_pool(engine);
                auto result = compute_integrals(this->world_,engine_pool,bs_array);
                return result;
            }
            else if(operation.get_operation()== Operation::Operations::cGTGCoulomb) {

                if(this->gtg_params_.empty()){
                    throw std::runtime_error("Gaussian Type Genminal Parameters are empty!");
                }

                libint2::TwoBodyEngine<libint2::cGTG_times_Coulomb> engine(max_nprim, static_cast<int>(max_am),0,std::numeric_limits<double>::epsilon(),this->gtg_params_);
                auto engine_pool = make_pool(engine);
                auto result = compute_integrals(this->world_,engine_pool,bs_array);
                return result;
            }
            else if(operation.get_operation() == Operation::Operations::cGTG){

                if(this->gtg_params_.empty()){
                    throw std::runtime_error("Gaussian Type Genminal Parameters are empty!");
                }

                libint2::TwoBodyEngine<libint2::cGTG> engine(max_nprim, static_cast<int>(max_am),0,std::numeric_limits<double>::epsilon(),this->gtg_params_);
                auto engine_pool = make_pool(engine);
                auto result = compute_integrals(this->world_,engine_pool,bs_array);
                return result;
            }
            else if(operation.get_operation() == Operation::Operations::cGTG2){

                if(this->gtg_params_.empty()){
                    throw std::runtime_error("Gaussian Type Genminal Parameters are empty!");
                }

                libint2::TwoBodyEngine<libint2::DelcGTG_square> engine(max_nprim, static_cast<int>(max_am),0,std::numeric_limits<double>::epsilon(),this->gtg_params_);
                auto engine_pool = make_pool(engine);
                auto result = compute_integrals(this->world_,engine_pool,bs_array);
                return result;
            }
            else {
                throw std::runtime_error("Invalid Two Body Operation");
            }

        }

    }

    template <typename Tile, typename Policy>
    typename AtomicIntegral<Tile,Policy>::TArray3 AtomicIntegral<Tile,Policy>::compute3(const std::wstring& formula_string) {
        Formula formula(formula_string);

        TA_USER_ASSERT(formula.notation() == Formula::Notation::Chemical, "Three Center Integral Must Use Chemical Notation");

        auto bra_indexs = formula.left_index();
        auto ket_indexs = formula.right_index();

        TA_ASSERT(bra_indexs.size() == 1);
        TA_ASSERT(ket_indexs.size() == 2);


        TA_ASSERT(bra_indexs[0].is_ao());
        TA_ASSERT(ket_indexs[0].is_ao());
        TA_ASSERT(ket_indexs[1].is_ao());

        auto bra_basis0 = this->index_to_basis(bra_indexs[0]);
        auto ket_basis0 = this->index_to_basis(ket_indexs[0]);
        auto ket_basis1 = this->index_to_basis(ket_indexs[1]);

        TA_ASSERT(bra_basis0 != nullptr);
        TA_ASSERT(ket_basis0 != nullptr);
        TA_ASSERT(ket_basis1 != nullptr);

        auto max_nprim = std::max({bra_basis0->max_nprim(), ket_basis0->max_nprim()
                                          ,ket_basis1->max_nprim()});
        auto max_am = std::max({bra_basis0->max_am(), ket_basis0->max_am(),
                                ket_basis1->max_am()});

        auto bs_array = tcc::utility::make_array(*bra_basis0, *ket_basis0, *ket_basis1);

        // convert operation to libint operator
        auto operation = formula.operation();
        if (operation.get_operation() == Operation::Operations::Coulomb) {
            libint2::TwoBodyEngine<libint2::Coulomb> engine(max_nprim, static_cast<int>(max_am));
            auto engine_pool = make_pool(engine);
            auto result = compute_integrals(this->world_,engine_pool,bs_array);
            return result;
        }
        else if(operation.get_operation()== Operation::Operations::cGTGCoulomb) {

            if(this->gtg_params_.empty()){
                throw std::runtime_error("Gaussian Type Genminal Parameters are empty!");
            }

            libint2::TwoBodyEngine<libint2::cGTG_times_Coulomb> engine(max_nprim, static_cast<int>(max_am),0,std::numeric_limits<double>::epsilon(),this->gtg_params_);
            auto engine_pool = make_pool(engine);
            auto result = compute_integrals(this->world_,engine_pool,bs_array);
            return result;
        }
        else if(operation.get_operation() == Operation::Operations::cGTG){

            if(this->gtg_params_.empty()){
                throw std::runtime_error("Gaussian Type Genminal Parameters are empty!");
            }

            libint2::TwoBodyEngine<libint2::cGTG> engine(max_nprim, static_cast<int>(max_am),0,std::numeric_limits<double>::epsilon(),this->gtg_params_);
            auto engine_pool = make_pool(engine);
            auto result = compute_integrals(this->world_,engine_pool,bs_array);
            return result;
        }
        else if(operation.get_operation() == Operation::Operations::cGTG2){

            if(this->gtg_params_.empty()){
                throw std::runtime_error("Gaussian Type Genminal Parameters are empty!");
            }

            libint2::TwoBodyEngine<libint2::DelcGTG_square> engine(max_nprim, static_cast<int>(max_am),0,std::numeric_limits<double>::epsilon(),this->gtg_params_);
            auto engine_pool = make_pool(engine);
            auto result = compute_integrals(this->world_,engine_pool,bs_array);
            return result;
        }
        else {
            throw std::runtime_error("Invalid Two Body Operation");
        }
    }

    template <typename Tile, typename Policy>
    typename AtomicIntegral<Tile,Policy>::TArray4 AtomicIntegral<Tile,Policy>::compute4(const std::wstring& formula_string) {
        Formula formula(formula_string);

        auto bra_indexs = formula.left_index();
        auto ket_indexs = formula.right_index();

        TA_ASSERT(bra_indexs.size() == 2);
        TA_ASSERT(ket_indexs.size() == 2);


        TA_ASSERT(bra_indexs[0].is_ao());
        TA_ASSERT(ket_indexs[0].is_ao());
        TA_ASSERT(bra_indexs[1].is_ao());
        TA_ASSERT(ket_indexs[1].is_ao());

        auto bra_basis0 = this->index_to_basis(bra_indexs[0]);
        auto ket_basis0 = this->index_to_basis(ket_indexs[0]);
        auto bra_basis1 = this->index_to_basis(bra_indexs[1]);
        auto ket_basis1 = this->index_to_basis(ket_indexs[1]);

        TA_ASSERT(bra_basis0 != nullptr);
        TA_ASSERT(ket_basis0 != nullptr);
        TA_ASSERT(bra_basis1 != nullptr);
        TA_ASSERT(ket_basis1 != nullptr);

        auto max_nprim = std::max({bra_basis0->max_nprim(), bra_basis1->max_nprim()
                                          ,ket_basis0->max_nprim(), ket_basis1->max_nprim()});
        auto max_am = std::max({bra_basis0->max_am(), bra_basis1->max_am(),
                                ket_basis0->max_am(),ket_basis1->max_am()});

        std::array<basis::Basis,4> bs_array;
        if (formula.notation() == Formula::Notation::Chemical){
            bs_array = {*bra_basis0, *bra_basis1, *ket_basis0, *ket_basis1};
        }
        else if(formula.notation() == Formula::Notation::Physical){
            bs_array = {*bra_basis0, *ket_basis0, *bra_basis1, *ket_basis1};
        }

        // convert operation to libint operator
        auto operation = formula.operation();
        if (operation.get_operation() == Operation::Operations::Coulomb) {
            libint2::TwoBodyEngine<libint2::Coulomb> engine(max_nprim, static_cast<int>(max_am));
            auto engine_pool = make_pool(engine);
            auto result = compute_integrals(this->world_,engine_pool,bs_array);
            return result;
        }
        else if(operation.get_operation()== Operation::Operations::cGTGCoulomb) {

            if(this->gtg_params_.empty()){
                throw std::runtime_error("Gaussian Type Genminal Parameters are empty!");
            }

            libint2::TwoBodyEngine<libint2::cGTG_times_Coulomb> engine(max_nprim, static_cast<int>(max_am),0,std::numeric_limits<double>::epsilon(),this->gtg_params_);
            auto engine_pool = make_pool(engine);
            auto result = compute_integrals(this->world_,engine_pool,bs_array);
            return result;
        }
        else if(operation.get_operation() == Operation::Operations::cGTG){

            if(this->gtg_params_.empty()){
                throw std::runtime_error("Gaussian Type Genminal Parameters are empty!");
            }

            libint2::TwoBodyEngine<libint2::cGTG> engine(max_nprim, static_cast<int>(max_am),0,std::numeric_limits<double>::epsilon(),this->gtg_params_);
            auto engine_pool = make_pool(engine);
            auto result = compute_integrals(this->world_,engine_pool,bs_array);
            return result;
        }
        else if(operation.get_operation() == Operation::Operations::cGTG2){

            if(this->gtg_params_.empty()){
                throw std::runtime_error("Gaussian Type Genminal Parameters are empty!");
            }

            libint2::TwoBodyEngine<libint2::DelcGTG_square> engine(max_nprim, static_cast<int>(max_am),0,std::numeric_limits<double>::epsilon(),this->gtg_params_);
            auto engine_pool = make_pool(engine);
            auto result = compute_integrals(this->world_,engine_pool,bs_array);
            return result;
        }
        else {
            throw std::runtime_error("Invalid Two Body Operation");
        }
    }
    //TODO R12 Integral

    // R12 Atomic Integral

    }
}


#endif //TILECLUSTERCHEM_ATOMIC_INTEGRAL_H

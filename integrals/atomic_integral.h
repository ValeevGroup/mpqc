//
// Created by Chong Peng on 10/14/15.
//

#ifndef TILECLUSTERCHEM_ATOMIC_INTEGRAL_H
#define TILECLUSTERCHEM_ATOMIC_INTEGRAL_H

#include <string>
#include <vector>
#include <iostream>
#include <cwchar>

#include "../common/namespaces.h"
#include "../ta_routines/sqrt_inv.h"
#include "../include/tiledarray.h"
#include "../basis/basis.h"
#include "../expression/formula.h"
#include "../expression/formula_registry.h"
#include "integral_engine_pool.h"
#include "task_integrals.h"
#include "../molecule/molecule.h"
#include "make_engine.h"
#include "../utility/make_array.h"
#include "../utility/wcout_utf8.h"
#include "../ta_routines/array_to_eigen.h"

namespace mpqc{
namespace integrals{


    //TODO return expression instead of array?
    class AtomicIntegralBase {
    public:

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

    protected:

        // get one body engine
        libint2::OneBodyEngine get_one_body_engine(const Operation& operation, int64_t max_nprim, int64_t max_am);

        // get two body engine kernel
        libint2::MultiplicativeSphericalTwoBodyKernel get_two_body_engine_kernel(const Operation &operation);


        // parse one body formula and set engine_pool and basis array
        void parse_one_body(const Formula& formula, std::shared_ptr<EnginePool<libint2::OneBodyEngine>>& engine_pool, Barray<2>& bases);

        // parse two body two center formula and set two body kernel and basis array
        void parse_two_body_two_center(const Formula& formula, libint2::MultiplicativeSphericalTwoBodyKernel& kernel, Barray<2>& bases, int64_t& max_nprim, int64_t& max_am);

        // parse two body three center formula and set two body kernel and basis array
        void parse_two_body_three_center(const Formula& formula, libint2::MultiplicativeSphericalTwoBodyKernel& kernel, Barray<3>& bases, int64_t& max_nprim, int64_t& max_am);

        // parse two body four center formula and set two body kernel and basis array
        void parse_two_body_four_center(const Formula& formula, libint2::MultiplicativeSphericalTwoBodyKernel& kernel, Barray<4>& bases, int64_t& max_nprim, int64_t& max_am);



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


    libint2::OneBodyEngine  AtomicIntegralBase::get_one_body_engine(const Operation &operation, int64_t max_nprim, int64_t max_am){

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

        return engine;
    }

    libint2::MultiplicativeSphericalTwoBodyKernel AtomicIntegralBase::get_two_body_engine_kernel(
            const Operation &operation)
    {
        libint2::MultiplicativeSphericalTwoBodyKernel kernel;
        if (operation.get_operation() == Operation::Operations::Coulomb) {
            kernel = libint2::Coulomb;
        }
        else if(operation.get_operation()== Operation::Operations::cGTGCoulomb) {
            kernel = libint2::cGTG_times_Coulomb;
        }
        else if(operation.get_operation() == Operation::Operations::cGTG){
            kernel = libint2::cGTG;
        }
        else if(operation.get_operation() == Operation::Operations::cGTG2){
            kernel = libint2::DelcGTG_square;
        }
        else {
            throw std::runtime_error("Invalid Two Body Operation");
        }
        return kernel;
    }


    void AtomicIntegralBase::parse_one_body(const Formula& formula, std::shared_ptr<EnginePool<libint2::OneBodyEngine>>& engine_pool, Barray<2>& bases){

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
        bases = tcc::utility::make_array(*bra_basis, *ket_basis);

        // convert operation to libint operator
        auto operation = formula.operation();
        engine_pool = make_pool(get_one_body_engine(operation,max_nprim,max_am));
    }

    void AtomicIntegralBase::parse_two_body_two_center(const Formula& formula, libint2::MultiplicativeSphericalTwoBodyKernel& kernel, Barray<2>& bases, int64_t& max_nprim, int64_t& max_am){

        TA_USER_ASSERT(formula.notation() == Formula::Notation::Chemical, "Two Body Two Center Integral Must Use Chemical Notation");

        auto bra_indexs = formula.left_index();
        auto ket_indexs = formula.right_index();

        TA_ASSERT(bra_indexs.size() == 1);
        TA_ASSERT(ket_indexs.size() == 1);

        auto bra_index0 = bra_indexs[0];
        auto ket_index0 = ket_indexs[0];

        TA_ASSERT(ket_index0.is_ao());
        TA_ASSERT(bra_index0.is_ao());

        auto bra_basis0 = this->index_to_basis(bra_index0);
        auto ket_basis0 = this->index_to_basis(ket_index0);

        TA_ASSERT(bra_basis0 != nullptr);
        TA_ASSERT(ket_basis0 != nullptr);

        max_nprim = std::max(bra_basis0->max_nprim(), ket_basis0->max_nprim());
        max_am = std::max(bra_basis0->max_am(), ket_basis0->max_am());
        bases = tcc::utility::make_array(*bra_basis0, *ket_basis0);

        auto operation = formula.operation();
        kernel = get_two_body_engine_kernel(operation);

    }

    void AtomicIntegralBase::parse_two_body_three_center(const Formula& formula, libint2::MultiplicativeSphericalTwoBodyKernel& kernel, Barray<3>& bases, int64_t& max_nprim, int64_t& max_am){

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

        max_nprim = std::max({bra_basis0->max_nprim(), ket_basis0->max_nprim()
                                     ,ket_basis1->max_nprim()});
        max_am = std::max({bra_basis0->max_am(), ket_basis0->max_am(),
                           ket_basis1->max_am()});

        bases = tcc::utility::make_array(*bra_basis0, *ket_basis0, *ket_basis1);

        // convert operation to libint operator
        auto operation = formula.operation();
        kernel = get_two_body_engine_kernel(operation);

    }

    void AtomicIntegralBase::parse_two_body_four_center(const Formula& formula, libint2::MultiplicativeSphericalTwoBodyKernel& kernel, Barray<4>& bases, int64_t& max_nprim, int64_t& max_am){
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

        max_nprim = std::max({bra_basis0->max_nprim(), bra_basis1->max_nprim()
                                     ,ket_basis0->max_nprim(), ket_basis1->max_nprim()});
        max_am = std::max({bra_basis0->max_am(), bra_basis1->max_am(),
                           ket_basis0->max_am(),ket_basis1->max_am()});

        if (formula.notation() == Formula::Notation::Chemical){
            bases = {*bra_basis0, *bra_basis1, *ket_basis0, *ket_basis1};
        }
        else if(formula.notation() == Formula::Notation::Physical){
            bases = {*bra_basis0, *ket_basis0, *bra_basis1, *ket_basis1};
        }

        // convert operation to libint operator
        auto operation = formula.operation();

        kernel = get_two_body_engine_kernel(operation);
    }

    /// Atomic Integral Class
    //// Op is a function
    /// Op will take TA::TensorD as argument and return Tile

    template<typename Tile, typename Policy>
    class AtomicIntegral : public AtomicIntegralBase{

    public:
        using TArray = TA::DistArray<Tile,Policy>;

        typedef Tile (*Op)(TA::TensorD &&);

        template<unsigned int N, typename E>
        using IntegralBuilder = integrals::IntegralBuilder<N,E,Op>;

        AtomicIntegral() = default;

        AtomicIntegral(madness::World& world,
                Op op,
                std::shared_ptr<molecule::Molecule> mol,
                std::shared_ptr<basis::Basis> obs,
                std::shared_ptr<basis::Basis> dfbs = nullptr,
                std::shared_ptr<basis::Basis> auxbs = nullptr,
                std::vector<std::pair<double,double>> gtg_params = std::vector<std::pair<double,double>>()
        ) : AtomicIntegralBase(world,mol,obs,dfbs,auxbs,gtg_params), op_(op), formula_registry_(){}

        virtual ~AtomicIntegral() = default;

        TArray compute(const std::wstring& );

        TArray compute_direct(const std::wstring& );

    protected:

        // compute two body integral
        template <typename Basis>
        TArray compute_two_body_integral(const libint2::MultiplicativeSphericalTwoBodyKernel& kernel, const Basis& basis, int64_t max_nprim, int64_t max_am);

        // compute integrals that has two dimension
        TArray compute2(const Formula& formula_string);
        // compute integrals that has three dimension
        TArray compute3(const Formula& formula_stirng);
        // compute integrals that has four dimension
        TArray compute4(const Formula& formula_string);


    private:

        //TODO direct integral
        //TODO return expression
        //TODO Screener for different type of integral
        // compute integral for sparse policy
        template <typename E, unsigned long N, typename U = Policy>
        TA::Array<double, N,Tile,typename std::enable_if<std::is_same<U,TA::SparsePolicy>::value, TA::SparsePolicy>::type> compute_integrals(
                madness::World& world, E const &engine, Barray<N> const &bases)
        {
            auto sreener = std::make_shared<Screener>(Screener{});
            auto result = mpqc::integrals::sparse_integrals(world,engine,bases,sreener,op_);
            return result;
        }

        // compute integral for dense policy
        template <typename E, unsigned long N, typename U = Policy>
        TA::Array<double, N,Tile,typename std::enable_if<std::is_same<U,TA::DensePolicy>::value, TA::DensePolicy>::type> compute_integrals(
                madness::World& world, E const &engine, Barray<N> const &bases)
        {
            auto sreener = std::make_shared<Screener>(Screener{});
            auto result = mpqc::integrals::dense_integrals(world,engine,bases,sreener,op_);
            return result;
        }

    private:
        FormulaRegistry<TArray> formula_registry_;
        Op op_;

    };


    template <typename Tile, typename Policy>
    template<typename Basis>
    typename AtomicIntegral<Tile, Policy>::TArray AtomicIntegral<Tile, Policy>::compute_two_body_integral(
            const libint2::MultiplicativeSphericalTwoBodyKernel& kernel, const Basis& bs_array, int64_t max_nprim,
            int64_t max_am) {

        typename AtomicIntegral<Tile,Policy>::TArray result;

        if(kernel == libint2::Coulomb){
            libint2::TwoBodyEngine<libint2::Coulomb> engine(max_nprim, static_cast<int>(max_am));
            auto engine_pool = make_pool(engine);
            result = compute_integrals(this->world_,engine_pool,bs_array);
        }
        else{
            if(this->gtg_params_.empty()){
                throw std::runtime_error("Gaussian Type Genminal Parameters are empty!");
            }

            if(kernel == libint2::cGTG){
                libint2::TwoBodyEngine<libint2::cGTG> engine(max_nprim, static_cast<int>(max_am),0,std::numeric_limits<double>::epsilon(),this->gtg_params_);
                auto engine_pool = make_pool(engine);
                result = compute_integrals(this->world_,engine_pool,bs_array);
            }
            else if(kernel == libint2::cGTG_times_Coulomb){
                libint2::TwoBodyEngine<libint2::cGTG_times_Coulomb> engine(max_nprim, static_cast<int>(max_am),0,std::numeric_limits<double>::epsilon(),this->gtg_params_);
                auto engine_pool = make_pool(engine);
                result = compute_integrals(this->world_,engine_pool,bs_array);
            }
            else if(kernel == libint2::DelcGTG_square){
                libint2::TwoBodyEngine<libint2::DelcGTG_square> engine(max_nprim, static_cast<int>(max_am),0,std::numeric_limits<double>::epsilon(),this->gtg_params_);
                auto engine_pool = make_pool(engine);
                result = compute_integrals(this->world_,engine_pool,bs_array);
            }

        }
        return result;
    }

    template <typename Tile, typename Policy>
    typename AtomicIntegral<Tile,Policy>::TArray AtomicIntegral<Tile,Policy>::compute(const std::wstring& formula_string) {

        Formula formula(formula_string);

        auto iter = formula_registry_.find(formula);

        if(iter != formula_registry_.end()){
            return iter->second;
        }else{

            if(formula.rank() == 2){
                auto result =  compute2(formula);
                formula_registry_.insert(formula, result);
                return result;
            }
            else if(formula.rank() == 3){
                auto result =  compute3(formula);
                formula_registry_.insert(formula, result);
                return result;
            }
            else if(formula.rank() == 4){
                auto result =  compute4(formula);
                formula_registry_.insert(formula, result);
                return result;
            }
        }

    }

    template <typename Tile, typename Policy>
    typename AtomicIntegral<Tile,Policy>::TArray AtomicIntegral<Tile,Policy>::compute2(const Formula& formula) {

        Barray<2> bs_array;

        // use one body engine
        if(formula.operation().is_onebody()){

            std::shared_ptr<EnginePool<libint2::OneBodyEngine>> engine_pool;
            parse_one_body(formula,engine_pool,bs_array);
            auto result = compute_integrals(this->world_,engine_pool,bs_array);

            std::cout << "Computed One Body Integral: ";
            wcout_utf8(formula.formula_string());
            std::cout << std::endl;
            return result;
        }
        // use two body engine
        else if(formula.operation().is_twobody()){

            int64_t max_nprim, max_am;
            libint2::MultiplicativeSphericalTwoBodyKernel kernel;

            parse_two_body_two_center(formula,kernel,bs_array,max_nprim,max_am);
            TA::DistArray<Tile,Policy> result = compute_two_body_integral( kernel, bs_array, max_nprim, max_am);

            std::cout << "Computed Twobody Two Center Integral: ";
            wcout_utf8(formula.formula_string());
            std::cout << std::endl;
            return result;
        }

    }

    template <typename Tile, typename Policy>
    typename AtomicIntegral<Tile,Policy>::TArray AtomicIntegral<Tile,Policy>::compute3(const Formula& formula) {

        libint2::MultiplicativeSphericalTwoBodyKernel kernel;
        Barray<3> bs_array;

        int64_t max_nprim, max_am;
        parse_two_body_three_center(formula,kernel,bs_array,max_nprim,max_am);

        TA::DistArray<Tile,Policy> result = compute_two_body_integral( kernel, bs_array, max_nprim, max_am);

        std::cout << "Computed Twobody Three Center Integral: ";
        wcout_utf8(formula.formula_string());
        std::cout << std::endl;
        return result;
    }

    template <typename Tile, typename Policy>
    typename AtomicIntegral<Tile,Policy>::TArray AtomicIntegral<Tile,Policy>::compute4(const Formula& formula) {

        if(formula.operation().has_option(Operation::Options::DensityFitting)){
            TA_ASSERT(this->dfbs_!= nullptr);


            // convert operation to libint operator
            auto operation = formula.operation();
            if (operation.get_operation() == Operation::Operations::Coulomb) {

                std::wstring two_center_formula_string = L"( Κ|G| Λ )";
                typename AtomicIntegral<Tile,Policy>::TArray two_center = this->compute2(two_center_formula_string);
                auto two_center_eigen = tcc::array_ops::array_to_eigen(two_center);

                MatrixD Leig = Eig::LLT<MatrixD>(two_center_eigen).matrixL();
                MatrixD two_center_inverse_sqrt_eigen = Leig.inverse();

                auto tr_V = two_center.trange().data()[0];
                auto two_center_inverse_sqrt = tcc::array_ops::eigen_to_array<Tile>(two_center.get_world(), two_center_inverse_sqrt_eigen, tr_V, tr_V);


                std::wstring three_center_formula_string;
                if (formula.notation() == Formula::Notation::Chemical){
                    auto ket0 = formula.right_index()[0].name();
                    auto ket1 = formula.right_index()[1].name();
                    three_center_formula_string = L"( Κ|G| " + ket0 + L" " + ket1 + L")";
                }
                else if(formula.notation() == Formula::Notation::Physical){
                    auto ket0 = formula.left_index()[0].name();
                    auto ket1 = formula.right_index()[0].name();
                    three_center_formula_string = L"( Κ|G| " + ket0 + L" " + ket1 + L")";
                }

                auto three_center = compute3(three_center_formula_string);

                three_center("X,i,j") = two_center_inverse_sqrt("K,X")*three_center("K,i,j");

                typename AtomicIntegral<Tile,Policy>::TArray result;

                if (formula.notation() == Formula::Notation::Chemical){
                    result("i,j,k,l") = three_center("K,i,j")*three_center("K,k,l");
                }
                else if(formula.notation() == Formula::Notation::Physical){
                    result("i,k,j,l") = three_center("K,i,j")*three_center("K,k,l");
                }
                std::cout << "Computed Two Body Four Center Density-Fitting Integral: ";
                wcout_utf8(formula.formula_string());
                std::cout << std::endl;
                return result;
            }
            else if(operation.get_operation()== Operation::Operations::cGTGCoulomb) {

                if(this->gtg_params_.empty()){
                    throw std::runtime_error("Gaussian Type Genminal Parameters are empty!");
                }

                std::wstring two_center_formula_string = L"( Κ|GR| Λ )";
                typename AtomicIntegral<Tile,Policy>::TArray two_center = this->compute2(two_center_formula_string);
//                std::cout << two_center;
                auto two_center_eigen = tcc::array_ops::array_to_eigen(two_center);

                MatrixD Leig = Eig::LLT<MatrixD>(two_center_eigen).matrixL();
                MatrixD two_center_inverse_sqrt_eigen = Leig.inverse();

                auto tr_V = two_center.trange().data()[0];
                auto two_center_inverse_sqrt = tcc::array_ops::eigen_to_array<Tile>(two_center.get_world(), two_center_inverse_sqrt_eigen, tr_V, tr_V);

                std::wstring three_center_formula_string;
                if (formula.notation() == Formula::Notation::Chemical){
                    auto ket0 = formula.right_index()[0].name();
                    auto ket1 = formula.right_index()[1].name();
                    three_center_formula_string = L"( Κ|GR| " + ket0 + L" " + ket1 + L")";
                }
                else if(formula.notation() == Formula::Notation::Physical){
                    auto ket0 = formula.left_index()[0].name();
                    auto ket1 = formula.right_index()[0].name();
                    three_center_formula_string = L"( Κ|GR| " + ket0 + L" " + ket1 + L")";
                }

                auto three_center = compute3(three_center_formula_string);

                three_center("X,i,j") = two_center_inverse_sqrt("X,K")*three_center("K,i,j");

                typename AtomicIntegral<Tile,Policy>::TArray result;

                if (formula.notation() == Formula::Notation::Chemical){
                    result("i,j,k,l") = three_center("K,i,j")*three_center("K,k,l");
                }
                else if(formula.notation() == Formula::Notation::Physical){
                    result("i,k,j,l") = three_center("K,i,j")*three_center("K,k,l");
                }
                std::cout << "Computed Two Body Four Center Density-Fitting Integral: ";
                wcout_utf8(formula.formula_string());
                std::cout << std::endl;
                return result;

            }
            else if(operation.get_operation() == Operation::Operations::cGTG){

                if(this->gtg_params_.empty()){
                    throw std::runtime_error("Gaussian Type Genminal Parameters are empty!");
                }
                std::wstring two_center_formula_string = L"( Κ|R| Λ )";
                typename AtomicIntegral<Tile,Policy>::TArray two_center = this->compute2(two_center_formula_string);
//                std::cout << two_center;
                auto two_center_eigen = tcc::array_ops::array_to_eigen(two_center);

                MatrixD Leig = Eig::LLT<MatrixD>(two_center_eigen).matrixL();
                MatrixD two_center_inverse_sqrt_eigen = Leig.inverse();

                auto tr_V = two_center.trange().data()[0];
                auto two_center_inverse_sqrt = tcc::array_ops::eigen_to_array<Tile>(two_center.get_world(), two_center_inverse_sqrt_eigen, tr_V, tr_V);

                std::wstring three_center_formula_string;
                if (formula.notation() == Formula::Notation::Chemical){
                    auto ket0 = formula.right_index()[0].name();
                    auto ket1 = formula.right_index()[1].name();
                    three_center_formula_string = L"( Κ|R| " + ket0 + L" " + ket1 + L")";
                }
                else if(formula.notation() == Formula::Notation::Physical){
                    auto ket0 = formula.left_index()[0].name();
                    auto ket1 = formula.right_index()[0].name();
                    three_center_formula_string = L"( Κ|R| " + ket0 + L" " + ket1 + L")";
                }

                auto three_center = compute3(three_center_formula_string);

                three_center("X,i,j") = two_center_inverse_sqrt("X,K")*three_center("K,i,j");

                typename AtomicIntegral<Tile,Policy>::TArray result;

                if (formula.notation() == Formula::Notation::Chemical){
                    result("i,j,k,l") = three_center("K,i,j")*three_center("K,k,l");
                }
                else if(formula.notation() == Formula::Notation::Physical){
                    result("i,k,j,l") = three_center("K,i,j")*three_center("K,k,l");
                }
                std::cout << "Computed Two Body Four Center Density-Fitting Integral: ";
                wcout_utf8(formula.formula_string());
                std::cout << std::endl;
                return result;

            }
            else if(operation.get_operation() == Operation::Operations::cGTG2){

                if(this->gtg_params_.empty()){
                    throw std::runtime_error("Gaussian Type Genminal Parameters are empty!");
                }

                std::wstring two_center_formula_string = L"( Κ|R2| Λ )";
                typename AtomicIntegral<Tile,Policy>::TArray two_center = this->compute2(two_center_formula_string);
//                std::cout << two_center;
                auto two_center_eigen = tcc::array_ops::array_to_eigen(two_center);

                MatrixD Leig = Eig::LLT<MatrixD>(two_center_eigen).matrixL();
                MatrixD two_center_inverse_sqrt_eigen = Leig.inverse();

                auto tr_V = two_center.trange().data()[0];
                auto two_center_inverse_sqrt = tcc::array_ops::eigen_to_array<Tile>(two_center.get_world(), two_center_inverse_sqrt_eigen, tr_V, tr_V);

                std::wstring three_center_formula_string;
                if (formula.notation() == Formula::Notation::Chemical){
                    auto ket0 = formula.right_index()[0].name();
                    auto ket1 = formula.right_index()[1].name();
                    three_center_formula_string = L"( Κ|R2| " + ket0 + L" " + ket1 + L")";
                }
                else if(formula.notation() == Formula::Notation::Physical){
                    auto ket0 = formula.left_index()[0].name();
                    auto ket1 = formula.right_index()[0].name();
                    three_center_formula_string = L"( Κ|R2| " + ket0 + L" " + ket1 + L")";
                }

                auto three_center = compute3(three_center_formula_string);

                three_center("X,i,j") = two_center_inverse_sqrt("X,K")*three_center("K,i,j");

                typename AtomicIntegral<Tile,Policy>::TArray result;

                if (formula.notation() == Formula::Notation::Chemical){
                    result("i,j,k,l") = three_center("K,i,j")*three_center("K,k,l");
                }
                else if(formula.notation() == Formula::Notation::Physical){
                    result("i,k,j,l") = three_center("K,i,j")*three_center("K,k,l");
                }
                std::cout << "Computed Two Body Four Center Density-Fitting Integral: ";
                wcout_utf8(formula.formula_string());
                std::cout << std::endl;
                return result;
            }
            else {
                throw std::runtime_error("Invalid Two Body Operation");
            }
        }
        else{

            int64_t max_nprim, max_am;
            libint2::MultiplicativeSphericalTwoBodyKernel kernel;
            Barray<4> bs_array;

            parse_two_body_four_center(formula,kernel,bs_array,max_nprim,max_am);

            TA::DistArray<Tile,Policy> result = compute_two_body_integral( kernel, bs_array, max_nprim, max_am);

            std::cout << "Computed Twobody Four Center Integral: ";
            wcout_utf8(formula.formula_string());
            std::cout << std::endl;
            return result;
        }

    }
    //TODO R12 Integral

}
}


#endif //TILECLUSTERCHEM_ATOMIC_INTEGRAL_H

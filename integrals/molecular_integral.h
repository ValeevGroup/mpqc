//
// Created by Chong Peng on 1/7/16.
//

#ifndef TILECLUSTERCHEM_MOLECULAR_INTEGRAL_H
#define TILECLUSTERCHEM_MOLECULAR_INTEGRAL_H

#include <string>
#include <vector>

#include "../include/tiledarray.h"
#include "../common/namespaces.h"
#include "../expression/formula.h"
#include "../expression/formula_registry.h"
#include "atomic_integral.h"
#include "../expression/orbital_registry.h"
#include "../ta_routines/diagonal_array.h"
#include "../scf/mo_build.h"


namespace mpqc{
namespace integrals{

// TODO MO transform that minimize operations by permutation
    /**
     * \brief Molecule Integral computation class
     *  This class computes molecular integral using  Formula
     *
     *  compute(formula) return TArray object
     *  (formula) return TArray expression
     *
     */
    template <typename Tile, typename Policy>
    class MolecularIntegral{
    public:
        using TArray = TA::DistArray<Tile, Policy>;
        using AtomicIntegral = AtomicIntegral<Tile,Policy>;

        /**
         *  Constructor
         *  @param atomic_integral  reference to AtomicIntegral class
         *  @param orbital_space_registry  shared pointer to OrbitalSpaceRegistry, which contain AO to MO coefficients
         *  @param formula_registry  FormulaRegistry used to store computed integral
         *  @param in rapidjson Document object
         *
         *
         *  Options in Input
         *  @param AccurateTime, bool, control if use fence in timing, default false
         */
        MolecularIntegral(AtomicIntegral &atomic_integral,
                          const std::shared_ptr<OrbitalSpaceRegistry<TArray>> orbital_space_registry,
                          const FormulaRegistry<TArray> &formula_registry,
                          const rapidjson::Document& in = rapidjson::Document()
                        )
                : world_(atomic_integral.get_world()), atomic_integral_(atomic_integral),
                  orbital_space_registry_(orbital_space_registry),
                  mo_formula_registry_(formula_registry)
        {
            atomic_integral_.set_orbital_space_registry(orbital_space_registry);
            parse_input(in);
        }


        /**
         *  Constructor
         *  @param atomic_integral  reference to AtomicIntegral class
         *  @param orbital_space_registry  shared pointer to OrbitalSpaceRegistry, which contain AO to MO coefficients
         *
         */
        MolecularIntegral(AtomicIntegral &atomic_integral,
                          const std::shared_ptr<OrbitalSpaceRegistry<TArray>> orbital_space_registry,
                          const rapidjson::Document& in = rapidjson::Document()
                        )
                : world_(atomic_integral.get_world()), atomic_integral_(atomic_integral),
                  orbital_space_registry_(orbital_space_registry),
                  mo_formula_registry_()
        {
            atomic_integral_.set_orbital_space_registry(orbital_space_registry);
            parse_input(in);
        }

        /// return reference to madness::World
        madness::World &get_world() const {
            return world_;
        }

        /// return reference to AtomicIntegral object
        AtomicIntegral& atomic_integral() const {
            return atomic_integral_;
        }

        /// return reference to AtomicIntegral object
        AtomicIntegral& atomic_integral() {
            return atomic_integral_;
        }

        /// wrapper to operator() function in AtomicIntegral
        TA::expressions::TsrExpr<TArray,true> atomic_integral (const std::wstring& str){
            return std::move(atomic_integral_(str));
        };

        /// reutrn shared pointer to OrbitalSpaceRegistry
        const std::shared_ptr<OrbitalSpaceRegistry<TArray>> orbital_space() const {
            return orbital_space_registry_;
        }

        /// return reference to FormulaRegistry
        const FormulaRegistry<TArray> &registry() const {
            return mo_formula_registry_;
        }

        /// return reference to FormulaRegistry
        FormulaRegistry<TArray> &registry() {
            return mo_formula_registry_;
        }

        /// wrapper to compute function
        TArray compute(const std::wstring& );


        /**
         *  compute integral by Formula
         *  this function will look into registry first
         *  if Formula computed, it will return it from registry
         *  if not, it will compute it
         */
        TArray compute(const Formula&);

        /// compute with str and return expression
        TA::expressions::TsrExpr<TArray,true> operator() (const std::wstring& str){
            auto formula = Formula(str);
            TArray array = compute(formula);
            auto& result = mo_formula_registry_.retrieve(formula);
            return result(formula.to_ta_expression());
        };

        /// remove all formula that has operation oper_str in mo_registry and ao_registry
        void remove_operation_all(madness::World& world, const std::wstring& oper_str){

            Operation operation(oper_str);
            Operation::Operations oper = operation.oper();

            mo_formula_registry_.remove_operation(world, oper);
            atomic_integral().registry().remove_operation(world, oper);
        }

    private:

        /// function to parse input
        void parse_input(const rapidjson::Document &in) {
            if(in.IsObject()){
                accurate_time_ = in.HasMember("AccurateTime") ? in["AccurateTime"].GetBool() : false;
            }
            else{
                accurate_time_ = false;
            }

            utility::print_par(world_, "\nConstructing Molecular Integral Class \n");
            utility::print_par(world_, "AccurateTime: " , accurate_time_, "\n");
            std::cout << std::endl;
        }
        /// compute integrals that has two dimension
        TArray compute2(const Formula& formula_string);

        /// compute integrals that has three dimension
        TArray compute3(const Formula& formula_string);

        /// compute integrals that has four dimension
        TArray compute4(const Formula& formula_string);

    private:

        /// find the corresponding AO formula, if index is already AO, it will be ignored
        Formula mo_to_ao(const Formula& formula);

        /// assert all index in formula are in MO
        void assert_all_mo(const Formula &formula);

    private:

        madness::World& world_;
        AtomicIntegral& atomic_integral_;
        std::shared_ptr<OrbitalSpaceRegistry<TArray>> orbital_space_registry_;
        FormulaRegistry<TArray> mo_formula_registry_;
        bool accurate_time_;
    };

    template <typename Tile, typename Policy>
    typename MolecularIntegral<Tile,Policy>::TArray MolecularIntegral<Tile,Policy>::compute2(const Formula &formula_string) {

        double time = 0.0;
        mpqc_time::t_point time0;
        mpqc_time::t_point time1;

        TArray result;
        // Identity matrix
        if(formula_string.operation().oper() == Operation::Operations::Identity){

            time0 = mpqc_time::now(world_,accurate_time_);

            auto left_index1 = formula_string.left_index()[0];
            auto right_index1 = formula_string.right_index()[0];
            auto left1 = orbital_space_registry_->retrieve(left_index1);
            auto right1 = orbital_space_registry_->retrieve(right_index1);

            //TODO better way to make model for diagonal matrix
            TArray tmp;
            tmp("i,j") = left1("k,i")*right1("k,j");

            // create diagonal array
            result = array_ops::create_diagonal_matrix(tmp,1.0);

            time1 = mpqc_time::now(world_,accurate_time_);
            time+= mpqc_time::duration_in_s(time0,time1);

            utility::print_par(world_, "Computed Identity: ");
            utility::wprint_par(world_, formula_string.formula_string());
            double size = utility::array_size(result);
            utility::print_par(world_," Size: ", size, " GB");
            utility::print_par(world_," Time: ", time, " s\n");
            return result;
        }


        // get AO
        auto ao_formula = mo_to_ao(formula_string);
        auto ao_integral = atomic_integral_.compute(ao_formula);

        time0 = mpqc_time::now(world_,accurate_time_);
        // convert to MO
        result = ao_integral;
        // get coefficient
        auto left_index1 = formula_string.left_index()[0];
        if(left_index1.is_mo()){
            auto left1 = orbital_space_registry_->retrieve(left_index1);
            result("i,r") = result("p,r")*left1("p,i");
        }
        auto right_index1 = formula_string.right_index()[0];
        if(right_index1.is_mo()){
            auto right1 = orbital_space_registry_->retrieve(right_index1);
            result("p,k") = result("p,r")*right1("r,k");
        }

        time1 = mpqc_time::now(world_,accurate_time_);
        time+= mpqc_time::duration_in_s(time0,time1);
        utility::print_par(world_, "Transformed MO Integral: ");
        utility::wprint_par(world_, formula_string.formula_string());
        double size = utility::array_size(result);
        utility::print_par(world_," Size: ", size, " GB");
        utility::print_par(world_," Time: ", time, " s\n");

        return result;
    }

    template <typename Tile, typename Policy>
    typename MolecularIntegral<Tile,Policy>::TArray MolecularIntegral<Tile,Policy>::compute3(const Formula &formula_string) {

        double time = 0.0;
        mpqc_time::t_point time0;
        mpqc_time::t_point time1;

        TArray result;
        // get AO
        auto ao_formula = mo_to_ao(formula_string);
        auto ao_integral = atomic_integral_.compute(ao_formula);

        // convert to MO, only convert the right side
        time0 = mpqc_time::now(world_,accurate_time_);

        // get coefficient
        auto right_index1 = formula_string.right_index()[0];
        if (right_index1.is_mo()) {
            auto right1 = orbital_space_registry_->retrieve(right_index1);
            result("K,i,q") = ao_integral("K,p,q") * right1("p,i");
        }
        auto right_index2 = formula_string.right_index()[1];
        if (right_index2.is_mo()) {
            auto right2 = orbital_space_registry_->retrieve(right_index2);
            result("K,p,j") = result("K,p,q") * right2("q,j");
        }

        time1 = mpqc_time::now(world_,accurate_time_);
        time+= mpqc_time::duration_in_s(time0,time1);

        utility::print_par(world_, "Transformed MO Integral: ");
        utility::wprint_par(world_, formula_string.formula_string());
        double size = utility::array_size(result);
        utility::print_par(world_," Size: ", size, " GB");
        utility::print_par(world_," Time: ", time, " s\n");

        return result;
    };

    template <typename Tile, typename Policy>
    typename MolecularIntegral<Tile,Policy>::TArray MolecularIntegral<Tile,Policy>::compute4(const Formula &formula_string) {

        double time = 0.0;
        mpqc_time::t_point time0;
        mpqc_time::t_point time1;
        TArray result;
        if(formula_string.operation().has_option(Operation::Options::DensityFitting)){

            // get df formula
            auto df_formulas = atomic_integral_.get_df_formula(formula_string);

            auto notation = formula_string.notation();
            // compute integral
            TArray left = compute(df_formulas[0]);

            TArray right = compute(df_formulas[2]);

            TArray center = atomic_integral_.compute(df_formulas[1]);

            time0 = mpqc_time::now(world_,accurate_time_);

            if(notation==Formula::Notation::Chemical){
                result("i,j,k,l") = left("q,i,j")*center("q,p")*right("p,k,l");
            }
            else{
                result("i,k,j,l") = left("q,i,j")*center("q,p")*right("p,k,l");
            }

            time1 = mpqc_time::now(world_,accurate_time_);
            time+= mpqc_time::duration_in_s(time0,time1);

        }else {
            // get AO
            auto ao_formula = mo_to_ao(formula_string);
            auto ao_integral = atomic_integral_.compute(ao_formula);

            // convert to MO
            time0 = mpqc_time::now(world_,accurate_time_);

            // get coefficient
            auto left_index1 = formula_string.left_index()[0];
            if (left_index1.is_mo()) {
                auto left1 = orbital_space_registry_->retrieve(left_index1);
                result("i,q,r,s") = ao_integral("p,q,r,s") * left1("p,i");
            }

            auto left_index2 = formula_string.left_index()[1];
            if (left_index2.is_mo()) {
                auto left2 = orbital_space_registry_->retrieve(left_index2);
                result("p,i,r,s") = result("p,q,r,s") * left2("q,i");
            }

            auto right_index1 = formula_string.right_index()[0];
            if (right_index1.is_mo()) {
                auto right1 = orbital_space_registry_->retrieve(right_index1);
                result("p,q,i,s") = result("p,q,r,s") * right1("r,i");
            }
            auto right_index2 = formula_string.right_index()[1];
            if (right_index2.is_mo()) {
                auto right2 = orbital_space_registry_->retrieve(right_index2);
                result("p,q,r,i") = result("p,q,r,s") * right2("s,i");
            }

            time1 = mpqc_time::now(world_,accurate_time_);
            time+= mpqc_time::duration_in_s(time0,time1);
        }

        utility::print_par(world_, "Transformed MO Integral: ");
        utility::wprint_par(world_, formula_string.formula_string());
        double size = utility::array_size(result);
        utility::print_par(world_," Size: ", size, " GB");
        utility::print_par(world_," Time: ", time, " s\n");

        return result;
    }

template <typename Tile, typename Policy>
Formula MolecularIntegral<Tile,Policy>::mo_to_ao(const Formula &formula) {

    std::vector<OrbitalIndex> ao_left_index, ao_right_index;

    auto left_index = formula.left_index();
    for(const auto& index : left_index){

        // find the correspoding ao index
        if(index.is_mo()){
            auto ao_index = orbital_space_registry_->retrieve(index).ao_key();
            ao_left_index.push_back(ao_index);
        }
        // if already ao, do nothing
        else{
            ao_left_index.push_back(index);
        }
    }

    auto right_index = formula.right_index();
    for(const auto& index : right_index){
        // find the correspoding ao index
        if(index.is_mo()){
            auto ao_index = orbital_space_registry_->retrieve(index).ao_key();
            ao_right_index.push_back(ao_index);
        }
        // if already ao, do nothing
        else{
            ao_right_index.push_back(index);
        }
    }

    // set formula with ao index
    auto ao_formula = formula;
    ao_formula.set_left_index(ao_left_index);
    ao_formula.set_right_index(ao_right_index);

    return ao_formula;
}

template <typename Tile, typename Policy>
void MolecularIntegral<Tile,Policy>::assert_all_mo(const Formula &formula) {

    auto left = formula.left_index();
    for(auto& index : left){
        TA_ASSERT(index.is_mo());
    }

    auto right = formula.right_index();
    for(auto& index : right){
        TA_ASSERT(index.is_mo());
    }

}
template <typename Tile, typename Policy>
typename MolecularIntegral<Tile,Policy>::TArray MolecularIntegral<Tile,Policy>::compute(const std::wstring &formula_string) {
    Formula formula(formula_string);
    return compute(formula);
}

template <typename Tile, typename Policy>
typename MolecularIntegral<Tile,Policy>::TArray MolecularIntegral<Tile,Policy>::compute(const Formula& formula) {

    auto iter = mo_formula_registry_.find(formula);

    TArray result;

    if(iter != mo_formula_registry_.end()){
        result = iter->second;
        utility::print_par(world_,"Retrived MO Integral: ");
        utility::wprint_par(world_, formula.formula_string());
        double size = utility::array_size(result);
        utility::print_par(world_," Size: ", size, " GB\n");
    }else{

        if(formula.rank() == 2){
            result =  compute2(formula);
            mo_formula_registry_.insert(formula, result);
        }
        else if(formula.rank() == 3){
            result =  compute3(formula);
            mo_formula_registry_.insert(formula, result);
        }
        else if(formula.rank() == 4){
            result =  compute4(formula);
            mo_formula_registry_.insert(formula, result);
        }
    }

    // make sure all processes obtained result and insert formula
//    world_.gop.fence();

    return result;
}
} // namespace integral
} // namespace mpqc

#endif //TILECLUSTERCHEM_MOLECULAR_INTEGRAL_H

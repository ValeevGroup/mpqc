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
#include "../expression/orbital_space_registry.h"
#include "../ta_routines/diagonal_array.h"


namespace mpqc{
namespace integrals{

    template <typename Tile, typename Policy>
    class MolecularIntegral{
    public:
        using TArray = TA::DistArray<Tile, Policy>;
        using AtomicIntegral = AtomicIntegral<Tile,Policy>;

        MolecularIntegral(AtomicIntegral &atomic_integral,
                          const std::shared_ptr<OrbitalSpaceRegistry<TArray>> orbital_space_registry,
                          const FormulaRegistry<TArray> &formula_registry)
                : world_(atomic_integral.get_world()), atomic_integral_(atomic_integral),
                  orbital_space_registry_(orbital_space_registry),
                  mo_formula_registry_(formula_registry)
        {
            atomic_integral_.set_orbital_space_registry(orbital_space_registry);
        }

        MolecularIntegral(AtomicIntegral &atomic_integral,
                          const std::shared_ptr<OrbitalSpaceRegistry<TArray>> orbital_space_registry)
                : world_(atomic_integral.get_world()), atomic_integral_(atomic_integral),
                  orbital_space_registry_(orbital_space_registry),
                  mo_formula_registry_()
        {
            atomic_integral_.set_orbital_space_registry(orbital_space_registry);
        }


        AtomicIntegral& atomic_integral() const {
            return atomic_integral_;
        }

        const std::shared_ptr<OrbitalSpaceRegistry<TArray>> orbital_space() const {
            return orbital_space_registry_;
        }

        const FormulaRegistry<TArray> &registry() const {
            return mo_formula_registry_;
        }

        TArray compute(const std::wstring& );
        TArray compute(const Formula&);

        TA::expressions::TsrExpr<TArray,true> operator() (const std::wstring& str){
            auto formula = Formula(str);
            TArray array = compute(formula);
            auto& result = mo_formula_registry_.retrieve(formula);
            return result(formula.to_ta_expression());
        };

    private:

        //TODO more operation F, J, K...
        // compute integrals that has two dimension
        TArray compute2(const Formula& formula_string);
        // compute integrals that has four dimension
        TArray compute4(const Formula& formula_string);

    private:

        Formula mo_to_ao(const Formula& formula);

        void assert_all_mo(const Formula &formula);

    private:

        madness::World& world_;
        AtomicIntegral& atomic_integral_;
        std::shared_ptr<OrbitalSpaceRegistry<TArray>> orbital_space_registry_;
        FormulaRegistry<TArray> mo_formula_registry_;
    };

    template <typename Tile, typename Policy>
    typename MolecularIntegral<Tile,Policy>::TArray MolecularIntegral<Tile,Policy>::compute2(const Formula &formula_string) {

        double time = 0.0;

        TArray result;
        // Identity matrix
        if(formula_string.operation().oper() == Operation::Operations::Identity){

            auto time0 = mpqc_time::fenced_now(world_);
            auto left_index1 = formula_string.left_index()[0];
            auto right_index1 = formula_string.right_index()[0];
            auto left1 = orbital_space_registry_->retrieve(left_index1);
            auto right1 = orbital_space_registry_->retrieve(right_index1);

            //TODO better way to make model for diagonal matrix
            TArray tmp;
            tmp("i,j") = left1("k,i")*right1("k,j");

            // create diagonal array
            result = ta_routines::create_diagonal_matrix(tmp,1.0);

            auto time1 = mpqc_time::fenced_now(world_);
            time+= mpqc_time::duration_in_s(time0,time1);
            utility::print_par(world_, "Computed Identity: ");
            utility::wprint_par(world_, formula_string.formula_string());
            utility::print_par(world_," Time: ", time, " s");
            double size = utility::array_size(result);
            utility::print_par(world_," Size: ", size, " GB\n");
            return result;
        }


        // get AO
        auto ao_formula = mo_to_ao(formula_string);
        auto ao_integral = atomic_integral_.compute(ao_formula);

        auto time0 = mpqc_time::fenced_now(world_);
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

        auto time1 = mpqc_time::fenced_now(world_);
        time+= mpqc_time::duration_in_s(time0,time1);
        utility::print_par(world_, "Transformed MO Integral: ");
        utility::wprint_par(world_, formula_string.formula_string());
        utility::print_par(world_," Time: ", time, " s");
        double size = utility::array_size(result);
        utility::print_par(world_," Size: ", size, " GB\n");

        return result;
    }

//TODO better inverse of two center
    template <typename Tile, typename Policy>
    typename MolecularIntegral<Tile,Policy>::TArray MolecularIntegral<Tile,Policy>::compute4(const Formula &formula_string) {

        double time = 0.0;
        TArray result;
        if(formula_string.operation().has_option(Operation::Options::DensityFitting)){
            // get AO
            auto ao_formula = mo_to_ao(formula_string);

            // get df formula
            auto df_formulas = atomic_integral_.get_df_formula(ao_formula);

            auto notation = formula_string.notation();
            // compute integral
            TArray left;
            {
                auto left_ao = atomic_integral_.compute(df_formulas[0]);

                auto time0 = mpqc_time::fenced_now(world_);

                left("i,j,k") = left_ao("i,j,k");

                OrbitalIndex left_index1;
                if(notation == Formula::Notation::Chemical){
                    left_index1 = formula_string.left_index()[1];
                }
                else{
                    left_index1 = formula_string.right_index()[0];
                }
                if (left_index1.is_mo()) {
                    auto left1 = orbital_space_registry_->retrieve(left_index1);
                    left("p,i,r") = left("p,q,r") * left1("q,i");
                }

                auto left_index2 = formula_string.left_index()[0];
                if (left_index2.is_mo()) {
                    auto left2 = orbital_space_registry_->retrieve(left_index2);
                    left("p,q,i") = left("p,q,r") * left2("r,i");
                }

                auto time1 = mpqc_time::fenced_now(world_);
                time+= mpqc_time::duration_in_s(time0,time1);
            }

            TArray right;
            {
                auto right_ao = atomic_integral_.compute(df_formulas[2]);

                auto time0 = mpqc_time::fenced_now(world_);
                right("i,j,k") = right_ao("i,j,k");

                OrbitalIndex right_index1;

                if(notation==Formula::Notation::Chemical){
                    right_index1 = formula_string.right_index()[0];
                }
                else{
                    right_index1 = formula_string.left_index()[1];
                }

                if (right_index1.is_mo()) {
                    auto right1 = orbital_space_registry_->retrieve(right_index1);
                    right("p,i,r") = right("p,q,r") * right1("q,i");
                }

                auto right_index2 = formula_string.right_index()[1];
                if (right_index2.is_mo()) {
                    auto right2 = orbital_space_registry_->retrieve(right_index2);
                    right("p,q,i") = right("p,q,r") * right2("r,i");
                }
                auto time1 = mpqc_time::fenced_now(world_);
                time+= mpqc_time::duration_in_s(time0,time1);
            }

            auto center = atomic_integral_.compute(df_formulas[1]);

            auto time0 = mpqc_time::fenced_now(world_);

            if(notation==Formula::Notation::Chemical){
                result("i,j,k,l") = left("q,j,i")*center("q,p")*right("p,k,l");
            }
            else{
                result("i,k,j,l") = left("q,j,i")*center("q,p")*right("p,k,l");
            }

            auto time1 = mpqc_time::fenced_now(world_);
            time+= mpqc_time::duration_in_s(time0,time1);

        }else {
            // get AO
            auto ao_formula = mo_to_ao(formula_string);
            auto ao_integral = atomic_integral_.compute(ao_formula);

            // convert to MO
            auto time0 = mpqc_time::fenced_now(world_);
            result = ao_integral;

            // get coefficient
            auto left_index1 = formula_string.left_index()[0];
            if (left_index1.is_mo()) {
                auto left1 = orbital_space_registry_->retrieve(left_index1);
                result("i,q,r,s") = result("p,q,r,s") * left1("p,i");
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

            auto time1 = mpqc_time::fenced_now(world_);
            time+= mpqc_time::duration_in_s(time0,time1);
        }

        utility::print_par(world_, "Transformed MO Integral: ");
        utility::wprint_par(world_, formula_string.formula_string());
        utility::print_par(world_," Time: ", time, " s");
        double size = utility::array_size(result);
        utility::print_par(world_," Size: ", size, " GB\n");

        return result;
    }

//TODO fix the string name
template <typename Tile, typename Policy>
Formula MolecularIntegral<Tile,Policy>::mo_to_ao(const Formula &formula) {

    std::vector<OrbitalIndex> ao_left_index, ao_right_index;

    auto left_index = formula.left_index();
    for(const auto& index : left_index){
        ao_left_index.push_back(index.mo_to_ao());
    }

    auto right_index = formula.right_index();
    for(const auto& index : right_index){
        ao_right_index.push_back(index.mo_to_ao());
    }

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
        return iter->second;
    }else{

        if(formula.rank() == 2){
            result =  compute2(formula);
            mo_formula_registry_.insert(formula, result);
            return result;
        }
        else if(formula.rank() == 4){
            result =  compute4(formula);
            mo_formula_registry_.insert(formula, result);
            return result;
        }
    }

    return result;
}
} // namespace integral
} // namespace mpqc

#endif //TILECLUSTERCHEM_MOLECULAR_INTEGRAL_H

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
                : atomic_integral_(atomic_integral),
                  orbital_space_registry_(orbital_space_registry),
                  mo_formula_registry_(formula_registry)
        {
            atomic_integral_.set_orbital_space_registry(orbital_space_registry);
        }

        MolecularIntegral(AtomicIntegral &atomic_integral,
                          const std::shared_ptr<OrbitalSpaceRegistry<TArray>> orbital_space_registry)
                : atomic_integral_(atomic_integral),
                  orbital_space_registry_(orbital_space_registry),
                  mo_formula_registry_()
        {
            atomic_integral_.set_orbital_space_registry(orbital_space_registry);
        }


        AtomicIntegral& atomic_integral() const {
            return atomic_integral_;
        }

        const FormulaRegistry<TArray> &registry() const {
            return mo_formula_registry_;
        }

        TArray compute(const std::wstring& );

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

        AtomicIntegral& atomic_integral_;
        std::shared_ptr<OrbitalSpaceRegistry<TArray>> orbital_space_registry_;
        FormulaRegistry<TArray> mo_formula_registry_;
    };

    template <typename Tile, typename Policy>
    typename MolecularIntegral<Tile,Policy>::TArray MolecularIntegral<Tile,Policy>::compute2(const Formula &formula_string) {

        // get AO
        auto ao_formula = mo_to_ao(formula_string);
        auto ao_integral = atomic_integral_.compute(ao_formula);

        // convert to MO
        TArray result = ao_integral;
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

        return result;
    }

//TODO better inverse of two center
    template <typename Tile, typename Policy>
    typename MolecularIntegral<Tile,Policy>::TArray MolecularIntegral<Tile,Policy>::compute4(const Formula &formula_string) {

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
            }

            TArray right;
            {
                auto right_ao = atomic_integral_.compute(df_formulas[2]);
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
            }

            auto center = compute(df_formulas[1]);
            //inverse two center integral
            auto center_eig = array_ops::array_to_eigen(center);
            MatrixD L_inv_eig = MatrixD(Eig::LLT<MatrixD>(center_eig).matrixL()).inverse();
            center_eig = L_inv_eig.transpose() * L_inv_eig;
            auto tr_center = center.trange().data()[0];
            center = array_ops::eigen_to_array<TA::TensorD>(center.get_world(), center_eig, tr_center, tr_center);


            TArray result;
            if(notation==Formula::Notation::Chemical){
                result("i,j,k,l") = left("q,j,i")*center("q,p")*right("p,k,l");
            }
            else{
                result("i,k,j,l") = left("q,j,i")*center("q,p")*right("p,k,l");
            }
            return result;


        }else {
            // get AO
            auto ao_formula = mo_to_ao(formula_string);
            auto ao_integral = atomic_integral_.compute(ao_formula);

            // convert to MO
            TArray result = ao_integral;

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

            return result;
        }
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

    auto iter = mo_formula_registry_.find(formula);

    if(iter != mo_formula_registry_.end()){
        return iter->second;
    }else{

        if(formula.rank() == 2){
            auto result =  compute2(formula);
            mo_formula_registry_.insert(formula, result);
            return result;
        }
        else if(formula.rank() == 4){
            auto result =  compute4(formula);
            mo_formula_registry_.insert(formula, result);
            return result;
        }
    }
}
} // namespace integral
} // namespace mpqc

#endif //TILECLUSTERCHEM_MOLECULAR_INTEGRAL_H

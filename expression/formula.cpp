//
// Created by Chong Peng on 10/15/15.
//

#include <algorithm>
#include <iostream>


#include "formula.h"

#include <boost/algorithm/string.hpp>
#include <TiledArray/error.h>

namespace mpqc{


    Formula::Formula(std::wstring formula_string) {
        if (formula_string.empty()) {
            throw std::runtime_error("Empty Formula!");
        }
        formula_ = formula_string;

        // detect the brackets <> or ()
        auto left_mark = [](const wchar_t letter) {
            return letter == L'<' || letter == L'(';
        };

        auto right_mark = [](const wchar_t letter) {
            return letter == L'>' || letter == L')';
        };

        // find formula between < > or ( )
        auto bra_symbol = std::find_if(formula_string.cbegin(), formula_string.cend(), left_mark);
        auto ket_symbol = std::find_if(formula_string.cbegin(), formula_string.cend(), right_mark);

        TA_ASSERT(bra_symbol < ket_symbol);

        if (bra_symbol == formula_string.cend() || ket_symbol == formula_string.cend()) {
            throw std::runtime_error("Formula should start with < or ( and end with > or ). \n");
        }

        TA_ASSERT( (*bra_symbol==L'<'&&*ket_symbol==L'>') || (*bra_symbol==L'(')&&*ket_symbol==L')');

        if(*bra_symbol==L'<'){
            notation_ = Notation::Physical;
        }else{
            notation_ = Notation::Chemical;
        }

        std::wstring formula(bra_symbol + 1, ket_symbol);

//        std::wcout << formula_string;


        // find option between [ and ]
        auto option_left = std::find(formula_string.cbegin(), formula_string.cend(), L'[');
        auto option_right = std::find(formula_string.cbegin(), formula_string.cend(), L']');

        std::wstring option;

        if( (option_left!=formula_string.cend()) && (option_right!=formula_string.cend())){
            TA_ASSERT(option_left < option_right);
            option = std::wstring(option_left+1,option_right);
        }


        // split the formula by |
        std::vector<std::wstring> split_formula;
        boost::split(split_formula, formula, boost::is_any_of(L"|"), boost::token_compress_on);

        if (split_formula.size() == 2) {
            std::wstring left_formula = std::move(split_formula[0]);
            std::wstring right_formula = std::move(split_formula[1]);

            operation_ = Operation(L" ",option);
            left_index_ = check_orbital_index(left_formula);
            right_index_ = check_orbital_index(right_formula);

        }
        else if (split_formula.size() == 3) {
            std::wstring left_formula = std::move(split_formula[0]);
            std::wstring operation = std::move(split_formula[1]);
            std::wstring right_formula = std::move(split_formula[2]);

            operation_ = Operation(operation, option);
            left_index_ = check_orbital_index(left_formula);
            right_index_ = check_orbital_index(right_formula);

        }
        else {
            throw std::runtime_error("Formula in wrong length. \n");
        }

        // error detecting

        // 1. density fitting only for four center code
        if(operation_.has_option(Operation::Options::DensityFitting) && (left_index_.size()!=2) && (right_index_.size()!=2)){
            throw std::runtime_error("Density Fitting is only available for Four Center Integral!");
        }

        // 2. one body operation
        if(operation_.is_onebody() && (left_index_.size()!=1) && (right_index_.size()!=1)){
            throw std::runtime_error("One body Operator with Wrong Index Size!");
        }

        // 3. more than three index
        if( (left_index_.size()>=3) || (right_index_.size()>=3)){
            throw std::runtime_error("Wrong Number of Index!");
        }

    }

    std::size_t Formula::rank() const {
        return (left_index_.size() + right_index_.size());
    }

    std::vector<OrbitalIndex> Formula::check_orbital_index(std::wstring index_array) {

        std::vector<std::wstring> split_index;
        boost::split(split_index,index_array,boost::is_any_of(L" "),boost::token_compress_on);

        if(split_index.empty()){
            return std::vector<OrbitalIndex>();
        }else{
            std::vector<OrbitalIndex> result;

            for(auto index : split_index){
                if(!index.empty()){
                    result.emplace_back(index);
                }
            }
            return result;
        }
    }

}

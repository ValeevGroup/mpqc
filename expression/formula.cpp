//
// Created by Chong Peng on 10/15/15.
//

#include "formula.h"
#include <boost/algorithm/string.hpp>

namespace mpqc{

    const std::unordered_map<std::wstring, Formula::Operation> Formula::string_to_operation = {
            {L"T", Operation::Kinetic},
            {L"V", Operation::Nuclear},
            {L"G", Operation::Coulomb},
            {L"R", Operation::cGTG },
            {L"GR", Operation::cGTGCoulomb },
            {L"R2", Operation::cGTG2 }
    };

    Formula::Formula(std::wstring formula) {
        if(formula.empty()){
            throw std::runtime_error("Empty Formula!");
        }
        formula_ = formula;
        // detect the brackets <>
        if(!(formula.front() == '<' && formula.back() == '>')){
            throw std::runtime_error("Formula should start with < and end with >. \n");
        }
        // remove brackets <>
        formula.pop_back();
        formula.erase(formula.begin());

        // split the string by |
        std::vector<std::wstring> split_formula;
        boost::split(split_formula,formula,boost::is_any_of("|"), boost::token_compress_on);

        if(split_formula.size() == 2){
            std::wstring left_formula = std::move(split_formula[0]);
            std::wstring right_formula = std::move(split_formula[1]);
            operation_ = Operation::Overlap;
        }
        else if(split_formula.size() == 3){
            std::wstring left_formula = std::move(split_formula[0]);
            std::wstring operation = std::move(split_formula[1]);
            std::wstring right_formula = std::move(split_formula[2]);

            operation_ = check_operation(operation);
            left_index_ = check_orbital_index(left_formula);
            right_index_ = check_orbital_index(right_formula);

        }
        else {
            throw std::runtime_error("Formula in wrong length. \n");
        }

    }

    Formula::Operation Formula::check_operation(std::wstring oper) {

        if(oper.empty()){
            throw std::runtime_error("Empty Operation!. \n");
        }

        boost::trim(oper);

        auto iter = string_to_operation.find(oper);

        if(iter == string_to_operation.end()){
            throw std::runtime_error("Operation is invalid operation!  \n");
        }

        auto operation = iter->second;

        return operation;
    }

    std::vector <OrbitalIndex> Formula::check_orbital_index(std::wstring index_array) {

        std::vector<std::wstring> split_index;
        boost::split(split_index,index_array,boost::is_any_of(" "),boost::token_compress_on);

        std::vector<OrbitalIndex> result;

        for(auto index : split_index){
            result.emplace_back(index);
        }

        return result;

    }
}

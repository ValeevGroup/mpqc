//
// Created by Chong Peng on 10/15/15.
//

#include "formula.h"
#include <boost/algorithm/string.hpp>

namespace mpqc{

    const std::unordered_map<std::string, Formula::Operation> Formula::string_to_operation = {
            {"K", Operation::Kinetic},
            {"V", Operation::Nuclear},
            {"G", Operation::Coulomb},
            {"R", Operation::cGTG },
            {"GR", Operation::cGTGCoulomb },
            {"R2", Operation::cGTG2 }
    };

    Formula::Formula(std::string formula) {
        if(formula.empty()){
            throw std::runtime_error("Empty Formula!");
        }
        formula_ = formula;
        // detect the brackets <>
        if(!(formula.front() == '<' && formula.back() == '>')){
            throw std::runtime_error("Formula should start with < and end with >. \n" + formula);
        }
        // remove brackets <>
        formula.pop_back();
        formula.erase(formula.begin());

        // split the string by |
        std::vector<std::string> split_formula;
        boost::split(split_formula,formula,boost::is_any_of("|"), boost::token_compress_on);

        if(split_formula.size() == 2){
            std::string left_formula = std::move(split_formula[0]);
            std::string right_formula = std::move(split_formula[1]);
            operation_ = Operation::Overlap;
        }
        else if(split_formula.size() == 3){
            std::string left_formula = std::move(split_formula[0]);
            std::string operation = std::move(split_formula[1]);
            std::string right_formula = std::move(split_formula[2]);

            operation_ = check_operation(operation);
            left_index_ = check_orbital_index(left_formula);
            right_index_ = check_orbital_index(right_formula);

        }
        else {
            throw std::runtime_error("Formula in wrong length. \n" + formula_);
        }

    }

    Formula::Operation Formula::check_operation(std::string oper) {

        if(oper.empty()){
            throw std::runtime_error("Empty Operation!. \n" + formula_);
        }

        auto iter = string_to_operation.find(oper);

        if(iter == string_to_operation.end()){
            throw std::runtime_error(oper + " is invalid operation!  \n");
        }

        auto operation = iter->second;

        return operation;
    }

    std::vector <OrbitalIndex> Formula::check_orbital_index(std::string index_array) {

        std::vector<std::string> split_index;
        boost::split(split_index,index_array,boost::is_any_of(" "),boost::token_compress_on);

        std::vector<OrbitalIndex> result;

        for(auto index : split_index){
            result.emplace_back(index);
        }

        return result;

    }
}

//
// Created by Chong Peng on 10/15/15.
//

#include "formula.h"
#include <boost/algorithm/string.hpp>

namespace mpqc{

    const std::unordered_map<std::wstring, Formula::Operation> Formula::one_body_operation = {
            {L"", Operation::Overlap},
            {L"T", Operation::Kinetic},
            {L"V", Operation::Nuclear}
    };

    const std::unordered_map<std::wstring, Formula::Operation> Formula::two_body_operation = {
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
        if(!(formula.front() == L'<' && formula.back() == L'>')){
            throw std::runtime_error("Formula should start with < and end with >. \n");
        }
        // remove brackets <>
        formula.pop_back();
        formula.erase(formula.begin());

        // split the string by |
        std::vector<std::wstring> split_formula;
        boost::split(split_formula,formula,boost::is_any_of(L"|"), boost::token_compress_on);

        if(split_formula.size() == 2){
            std::wstring left_formula = std::move(split_formula[0]);
            std::wstring right_formula = std::move(split_formula[1]);

            operation_ = Operation::Overlap;
            left_index_ = check_orbital_index(left_formula);
            right_index_ = check_orbital_index(right_formula);

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

        if (oper.empty()) {
            throw std::runtime_error("Empty Operation!. \n");
        }

        boost::trim(oper);

        auto iter1 = one_body_operation.find(oper);

        if (iter1 == one_body_operation.end()) {
            auto iter2 = two_body_operation.find(oper);

            if(iter2 != two_body_operation.end()){
                return iter2->second;
            }else{
                throw std::runtime_error("Invalid Operation");
            }
        }
        else{
            return iter1->second;
        }

    }

    std::vector<OrbitalIndex> Formula::check_orbital_index(std::wstring index_array) {

        std::vector<std::wstring> split_index;
        boost::split(split_index,index_array,boost::is_any_of(L" "),boost::token_compress_on);

        std::vector<OrbitalIndex> result;

        for(auto index : split_index){
            result.emplace_back(index);
        }

        return result;

    }

    bool Formula::is_onebody() const {

        for (auto it = one_body_operation.begin(); it != one_body_operation.end(); ++it ){
            if (it->second == operation_){
                return true;
            }
        }
        return false;
    }

    bool Formula::is_twobody() const {

        for (auto it = two_body_operation.begin(); it != two_body_operation.end(); ++it ){
            if (it->second == operation_){
                return true;
            }
        }
        return false;
    }
}

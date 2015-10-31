//
// Created by Chong Peng on 10/31/15.
//

#include "operation.h"

#include <boost/algorithm/string.hpp>


namespace mpqc{
    using Operations = mpqc::Operation::Operations;

    const std::unordered_map<std::wstring, Operations> Operation::one_body_operation = {
        {L"", Operations::Overlap},
        {L"T", Operations::Kinetic},
        {L"V", Operations::Nuclear}
    };

    const std::unordered_map<std::wstring, Operations> Operation::two_body_operation = {
        {L"G", Operations::Coulomb},
        {L"R", Operations::cGTG },
        {L"GR", Operations::cGTGCoulomb },
        {L"R2", Operations::cGTG2 }
    };

    bool Operation::is_onebody() const {

        for (auto it = one_body_operation.begin(); it != one_body_operation.end(); ++it ){
            if (it->second == operation_){
                return true;
            }
        }
        return false;
    }

    bool Operation::is_twobody() const {

        for (auto it = two_body_operation.begin(); it != two_body_operation.end(); ++it ){
            if (it->second == operation_){
                return true;
            }
        }
        return false;
    }

    Operation::Operation(std::wstring oper) {
        if (oper.empty()) {
            throw std::runtime_error("Empty Operation!. \n");
        }

        boost::trim(oper);

        auto iter1 = one_body_operation.find(oper);

        if (iter1 == one_body_operation.end()) {
            auto iter2 = two_body_operation.find(oper);

            if(iter2 != two_body_operation.end()){
                operation_ = iter2->second;
            }else{
                throw std::runtime_error("Invalid Operation");
            }
        }
        else{
            operation_ = iter1->second;
        }
    }
}

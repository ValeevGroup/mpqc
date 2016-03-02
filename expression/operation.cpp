//
// Created by Chong Peng on 10/31/15.
//

#include "operation.h"

#include <boost/algorithm/string.hpp>


namespace mpqc{
    using Operations = mpqc::Operation::Operations;
    using Options = mpqc::Operation::Options;

    const std::unordered_map<std::wstring, Operations> Operation::one_body_operation = {
        {L"", Operations::Overlap},
        {L"T", Operations::Kinetic},
        {L"V", Operations::Nuclear}
    };

    const std::unordered_map<std::wstring, Operations> Operation::two_body_operation = {
        {L"G", Operations::Coulomb},
        {L"R", Operations::cGTG },
        {L"GR", Operations::cGTGCoulomb },
        {L"R2", Operations::cGTG2 },
        {L"dR2", Operations::DelcGTG2}
    };

    const std::unordered_map<std::wstring, Operations> Operation::fock_operation = {
        {L"J", Operations::J},
        {L"K", Operations::K },
        {L"K(α)", Operations::KAlpha },
        {L"K(β)", Operations::KBeta },
        {L"F", Operations::Fock},
        {L"F(α)", Operations::FockAlpha},
        {L"F(β)", Operations::FockBeta},
    };

    const std::map<Operations, std::wstring> Operation::operation_to_string = {
            {Operations::Overlap, L""},
            {Operations::Kinetic, L"T"},
            {Operations::Nuclear, L"V"},
            {Operations::Coulomb, L"G"},
            {Operations::cGTG, L"R"},
            {Operations::cGTG2, L"R2"},
            {Operations::cGTGCoulomb, L"GR" },
            {Operations::DelcGTG2, L"dR2"}
    };

    const std::unordered_map<std::wstring, Options> Operation::option = {
            {L"df", Options::DensityFitting},
            {L"inv", Options::Inverse },
            {L"inv_sq", Options::InverseSquareRoot }
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

    bool Operation::is_fock() const {
        return (operation_==Operations::Fock || operation_==Operations::FockAlpha || operation_==Operations::FockBeta);
    }

    bool Operation::is_jk() const {
        return (operation_== Operations::J || operation_==Operations::K || operation_==Operations::KAlpha ||
                operation_==Operations::KBeta);
    }

    bool Operation::is_r12() const {
        return (operation_==Operations::cGTG2 || operation_ == Operations::cGTG2 || operation_==Operations::cGTGCoulomb
                || operation_==Operations::DelcGTG2);
    }

    bool Operation::has_option(Options op) const {
        auto df = std::find(options_.cbegin(),options_.cend(),op);
        return (df != options_.cend());
    }

    Operation::Operation(std::wstring oper, std::wstring opt) {
        // parse operation
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

                auto iter3 = fock_operation.find(oper);
                if(iter3 != fock_operation.end()){
                    operation_ = iter3->second;
                }
                else{
                    throw std::runtime_error("Invalid Operation");
                }
            }

        }
        else{
            operation_ = iter1->second;
        }

        // parse option

        // split string by , and space
        if(!opt.empty()){
            std::vector<std::wstring> split_option;
            boost::split(split_option,opt,boost::is_any_of(L", "), boost::token_compress_on);

            std::vector<Options> result;

            for(const auto& op : split_option){
                auto iter = option.find(op);
                if(iter == option.end()){
                    throw std::runtime_error("Invalid Option");
                }
                else{
                    result.push_back(iter->second);
                }
            }

            std::sort(result.begin(),result.end());
            options_ = result;
        }
    }

    bool Operation::operator==(const Operation &other) const {

        bool same_operation = (this->operation_ == other.operation_);

        bool same_option = false;

        if(options_.size() == other.options_.size()){
            same_option = std::equal(options_.begin(), options_.end(), other.options_.begin());
        }

        return same_operation && same_option;
    }

    bool Operation::operator<(const Operation &other) const {
        return other.oper() == this->oper() ? this->options() < other.options() : this->oper() <
                                                                                  other.oper();
    }

}

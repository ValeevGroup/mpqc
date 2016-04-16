//
// Created by Chong Peng on 10/31/15.
//

#include "operation.h"

#include <boost/algorithm/string.hpp>
#include <string>
#include <memory>


namespace mpqc{
    using Operations = mpqc::Operation::Operations;
    using Options = mpqc::Operation::Options;

    const std::unordered_map<std::wstring, Operations> Operation::one_body_operation = {
        {L"", Operations::Overlap},
        {L"T", Operations::Kinetic},
        {L"V", Operations::Nuclear},
        {L"H", Operations::Core},
        {L"I", Operations::Identity }
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
        {L"hJ", Operations::hJ},
        {L"K", Operations::K },
        {L"K(α)", Operations::KAlpha },
        {L"K(β)", Operations::KBeta },
        {L"F", Operations::Fock},
        {L"F(α)", Operations::FockAlpha},
        {L"F(β)", Operations::FockBeta}
    };

    const std::map<Operations, std::wstring> Operation::oper_to_string = {
            {Operations::Overlap, L""},
            {Operations::Kinetic, L"T"},
            {Operations::Nuclear, L"V"},
            {Operations::Core, L"H"},
            {Operations::Coulomb, L"G"},
            {Operations::cGTG, L"R"},
            {Operations::cGTG2, L"R2"},
            {Operations::cGTGCoulomb, L"GR" },
            {Operations::DelcGTG2, L"dR2"},
            {Operations::J, L"J"},
            {Operations::hJ, L"hJ"},
            {Operations::K, L"K"},
            {Operations::KAlpha, L"K(α)"},
            {Operations::KBeta, L"K(β)"},
            {Operations::Fock, L"F"},
            {Operations::FockAlpha, L"F(α)"},
            {Operations::FockBeta, L"F(β)" },
            {Operations::Identity, L"I"}
    };

    const std::map<Options, std::wstring> Operation::option_to_string = {
            {Options::DensityFitting, L"df" },
            {Options::Inverse, L"inv"},
            {Options::InverseSquareRoot, L"inv_sqr"}
    };

    const std::unordered_map<std::wstring, Options> Operation::option = {
            {L"df", Options::DensityFitting},
            {L"inv", Options::Inverse },
            {L"inv_sqr", Options::InverseSquareRoot }
    };

    bool Operation::is_onebody() const {

        for (auto it = one_body_operation.begin(); it != one_body_operation.end(); ++it ){
            if (it->second == oper_){
                return true;
            }
        }
        return false;
    }

    bool Operation::is_twobody() const {

        for (auto it = two_body_operation.begin(); it != two_body_operation.end(); ++it ){
            if (it->second == oper_){
                return true;
            }
        }
        return false;
    }

    bool Operation::is_fock() const {
        return (oper_==Operations::Fock || oper_==Operations::FockAlpha || oper_==Operations::FockBeta);
    }

    bool Operation::is_jk() const {
        return (oper_== Operations::J || oper_==Operations::K || oper_==Operations::KAlpha ||
                oper_==Operations::KBeta);
    }

    bool Operation::is_r12() const {
        return (oper_==Operations::cGTG2 || oper_ == Operations::cGTG || oper_==Operations::cGTGCoulomb
                || oper_==Operations::DelcGTG2);
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
                oper_ = iter2->second;
            }else{

                auto iter3 = fock_operation.find(oper);
                if(iter3 != fock_operation.end()){
                    oper_ = iter3->second;
                }
                else{
                    throw std::runtime_error("Invalid Operation");
                }
            }

        }
        else{
            oper_ = iter1->second;
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

        bool same_operation = (this->oper_ == other.oper_);

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

    const std::wstring Operation::oper_string() const {
        const auto result = oper_to_string.find(oper_);
        return result->second;
    }

    const std::wstring Operation::option_string() const {
        std::wstring result;
        if(options_.empty()){
            return  result;
        }
        for(const auto& option : options_){
            result += option_to_string.find(option)->second + L",";
        }
        result = L"[" + result;
        result.back() = L']';
        return result;
    }
}

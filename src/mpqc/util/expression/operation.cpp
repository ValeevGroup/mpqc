//
// Created by Chong Peng on 10/31/15.
//


#include <boost/algorithm/string.hpp>
#include <string>
#include <memory>

#include <mpqc/util/expression/operation.h>

namespace mpqc{
    using Operations = mpqc::Operation::Operators;
    using Options = mpqc::Operation::Options;

    const std::unordered_map<std::wstring, Operations> Operation::one_body_operation = {
        {L"", Operators::Overlap},
        {L"T", Operators::Kinetic},
        {L"V", Operators::Nuclear},
        {L"H", Operators::Core},
        {L"I", Operators::Identity }
    };

    const std::unordered_map<std::wstring, Operations> Operation::two_body_operation = {
        {L"G", Operators::Coulomb},
        {L"R", Operators::cGTG },
        {L"GR", Operators::cGTGCoulomb },
        {L"R2", Operators::cGTG2 },
        {L"dR2", Operators::DelcGTG2}
    };

    const std::unordered_map<std::wstring, Operations> Operation::fock_operation = {
        {L"J", Operators::J},
        {L"hJ", Operators::hJ},
        {L"K", Operators::K },
        {L"K(α)", Operators::KAlpha },
        {L"K(β)", Operators::KBeta },
        {L"F", Operators::Fock},
        {L"F(α)", Operators::FockAlpha},
        {L"F(β)", Operators::FockBeta}
    };

    const std::map<Operations, std::wstring> Operation::oper_to_string = {
            {Operators::Overlap, L""},
            {Operators::Kinetic, L"T"},
            {Operators::Nuclear, L"V"},
            {Operators::Core, L"H"},
            {Operators::Coulomb, L"G"},
            {Operators::cGTG, L"R"},
            {Operators::cGTG2, L"R2"},
            {Operators::cGTGCoulomb, L"GR" },
            {Operators::DelcGTG2, L"dR2"},
            {Operators::J, L"J"},
            {Operators::hJ, L"hJ"},
            {Operators::K, L"K"},
            {Operators::KAlpha, L"K(α)"},
            {Operators::KBeta, L"K(β)"},
            {Operators::Fock, L"F"},
            {Operators::FockAlpha, L"F(α)"},
            {Operators::FockBeta, L"F(β)" },
            {Operators::Identity, L"I"}
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
        return (oper_==Operators::Fock || oper_==Operators::FockAlpha || oper_==Operators::FockBeta);
    }

    bool Operation::is_jk() const {
        return (oper_== Operators::J || oper_==Operators::K || oper_==Operators::KAlpha ||
                oper_==Operators::KBeta);
    }

    bool Operation::is_r12() const {
        return (oper_==Operators::cGTG2 || oper_ == Operators::cGTG || oper_==Operators::cGTGCoulomb
                || oper_==Operators::DelcGTG2);
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

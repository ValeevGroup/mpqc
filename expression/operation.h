//
// Created by Chong Peng on 10/31/15.
//

#ifndef TILECLUSTERCHEM_OPERATION_H
#define TILECLUSTERCHEM_OPERATION_H

#include <unordered_map>
#include <vector>
#include <string>
#include <map>

namespace mpqc{

    class Operation{
    public:
        enum class Operations{
            Overlap = 0,
            Kinetic = 1,
            Nuclear = 2,
            Coulomb = 3,
            cGTG = 4,
            cGTG2 = 5,
            cGTGCoulomb = 6,
            DelcGTG2 = 7,
            J = 8,
            K = 9,
            KAlpha = 10,
            KBeta = 11,
            F = 12,
            FAlpha = 13,
            FBeta = 14
        };

        enum class Options{
            DensityFitting = 0,
            InverseSquareRoot = 1
        };

        static const std::unordered_map<std::wstring, Operations> one_body_operation;
        static const std::unordered_map<std::wstring, Operations> two_body_operation;
        static const std::map<Operations,std::wstring> operation_to_string;
        static const std::unordered_map<std::wstring, Options> option;

        Operation() = default;
        Operation(Operation const &) = default;
        Operation(Operation &&) = default;
        Operation& operator=(Operation const &) = default;
        Operation& operator=(Operation &&) = default;

        Operation(std::wstring operation, std::wstring option = L"");

        const std::vector<Options> get_options() const {
            return options_;
        }

        const Operations &get_operation() const {
            return operation_;
        }

        const std::wstring string() const{
            const auto result = operation_to_string.find(operation_);
            return result->second;
        }

        bool has_option(Options op) const;

        bool is_onebody() const;
        bool is_twobody() const;
        bool is_r12() const;

        bool operator==(Operation const & other) const;

        bool operator!=(Operation const & other) const{
            return !(*this==other);
        }

        bool operator<(const Operation& other) const;

    private:

        Operations operation_;
        std::vector<Options> options_;

    };


}



#endif //TILECLUSTERCHEM_OPERATION_H

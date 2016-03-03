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
            Fock = 12,
            FockAlpha = 13,
            FockBeta = 14
        };

        enum class Options{
            DensityFitting = 0,
            Inverse = 1,
            InverseSquareRoot = 2
        };

        static const std::unordered_map<std::wstring, Operations> one_body_operation;
        static const std::unordered_map<std::wstring, Operations> two_body_operation;
        static const std::unordered_map<std::wstring, Operations> fock_operation;
        static const std::map<Operations,std::wstring> oper_to_string;
        static const std::map<Options,std::wstring> option_to_string;
        static const std::unordered_map<std::wstring, Options> option;

        Operation() = default;
        Operation(Operation const &) = default;
        Operation(Operation &&) = default;
        Operation& operator=(Operation const &) = default;
        Operation& operator=(Operation &&) = default;

        Operation(std::wstring operation, std::wstring option = L"");

        const std::vector<Options> options() const {
            return options_;
        }

        const Operations &oper() const {
            return oper_;
        }

        Operations &oper() {
            return oper_;
        }

        void set_oper(const Operations &oper) {
            Operation::oper_ = oper;
        }

        const std::wstring oper_string() const;
        const std::wstring option_string() const;

        bool has_option(Options op) const;

        bool is_onebody() const;
        bool is_twobody() const;
        bool is_fock() const;
        bool is_jk() const;
        bool is_r12() const;

        bool operator==(Operation const & other) const;

        bool operator!=(Operation const & other) const{
            return !(*this==other);
        }

        bool operator<(const Operation& other) const;

    private:

        Operations oper_;
        std::vector<Options> options_;

    };


}



#endif //TILECLUSTERCHEM_OPERATION_H

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

    /**
     * \brief Class to represent operation used in Formula
     *
     * Dictionary of wstring to operations
     *  - "" -> Overlap
     *  - T -> Kinetic
     *  - V -> Nuclear
     *  - H -> Hcore
     *  - I -> Indentity
     *  - G -> Coulomb(libint::Coulomb)
     *  - R -> cGTG(libint::cGTG)
     *  - GR -> cGTGCoulomb(libint::cGTG_times_Coulomb)
     *  - R2 -> cGTG2(libint::cGTG, parameters will get squared)
     *  - dR2 -> DelcGTG2(libint::DelcGTG_square)
     *  - J -> J(Coulomb)
     *  - hJ -> hJ (H + J)
     *  - K -> K(Exchange)
     *  - K(α) -> KAlpha(Exchange for Alpha Spin)
     *  - K(β) -> KBeta(Exchange for Beta Spin)
     *  - F -> Fock
     *  - F(α) -> FockAlpha(Fock for Alpha Spin)
     *  - F(β) -> FockBeta(Fock for Beta Spin)
     *
     *  Dictionary of wstring to options
     *  - df -> DensityFitting
     *  - inv -> Inverse
     *  - inv_sqr -> InverseSquareRoot
     */

    class Operation{
    public:
        /**
         *  Operations types
         */
        enum class Operations{
            Overlap = 0,
            Kinetic = 1,
            Nuclear = 2,
            Core = 3,
            Coulomb = 4,
            cGTG = 5,
            cGTG2 = 6,
            cGTGCoulomb = 7,
            DelcGTG2 = 8,
            J = 9,
            K = 10,
            KAlpha = 11,
            KBeta = 12,
            Fock = 13,
            FockAlpha = 14,
            FockBeta = 15,
            hJ = 16,
            Identity = 17
        };

        /**
         * Options types
         */
        enum class Options{
            DensityFitting = 0,
            Inverse = 1,
            InverseSquareRoot = 2
        };

        /**
         * maps of string to operations and option
         * vice versa
         */
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

        /**
         * Constructor
         *
         * @param operation  string of operation
         * @param option string of a list of option separated by ",", it will sort options so that they
         * always come in the same order
         */
        Operation(std::wstring operation, std::wstring option = L"");

        /// return options options_
        const std::vector<Options> options() const {
            return options_;
        }

        /// return operation oper_
        const Operations &oper() const {
            return oper_;
        }

        /// return operation oper_
        Operations &oper() {
            return oper_;
        }

        /// set operation oper_
        void set_oper(const Operations &oper) {
            Operation::oper_ = oper;
        }

        /// return string that correspond to oper_
        const std::wstring oper_string() const;

        /// return string that correspond to options_ wraped in []
        const std::wstring option_string() const;

        /// return true if have Options op
        bool has_option(Options op) const;

        /// check if oper_ is one body operation
        bool is_onebody() const;

        /// check if oper_ is two body operation
        bool is_twobody() const;

        /// check if oper_ is fock operation
        bool is_fock() const;

        /// check if oper_ is J or K operation
        bool is_jk() const;

        /// check if oper_ is R12 operation
        bool is_r12() const;

        /// equality check by comparing oper_ and options_
        bool operator==(Operation const & other) const;

        /// equality check by comparing oper_ and options_
        bool operator!=(Operation const & other) const{
            return !(*this==other);
        }

        /// comparison by comparing oper_ and options_
        bool operator<(const Operation& other) const;

    private:

        Operations oper_;
        std::vector<Options> options_;

    };


}



#endif //TILECLUSTERCHEM_OPERATION_H

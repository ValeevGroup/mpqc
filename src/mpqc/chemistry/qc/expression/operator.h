//
// Created by Chong Peng on 10/31/15.
//

#ifndef SRC_MPQC_UTIL_EXPRESSION_OPERATOR_H_
#define SRC_MPQC_UTIL_EXPRESSION_OPERATOR_H_

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
     *  - I -> Identity
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

    class Operator {
    public:
        /**
         *  Operator types
         */
        enum class Type {
            __first_1body_operator = 0,
            Identity = 0,
            Overlap = 1,
            Kinetic = 2,
            Nuclear = 3,
            Core = 4,
            __last_1body_operator = 4,
            __first_fock_operator = 32,
            J = 32,
            K = 33,
            KAlpha = 34,
            KBeta = 35,
            Fock = 36,
            FockAlpha = 37,
            FockBeta = 38,
            hJ = 39,
            __last_fock_operator = 39,
            __first_2body_operator = 128,
            Coulomb = 128,
            cGTG = 129,
            cGTG2 = 130,
            cGTGCoulomb = 131,
            DelcGTG2 = 132,
            __last_2body_operator = 132
        };

        /**
         * Option types
         */
        enum class Option {
            DensityFitting = 0,
            Inverse = 1,
            InverseSquareRoot = 2
        };

        /**
         * maps of string to operations and option
         * vice versa
         */
        static const std::unordered_map<std::wstring, Type> one_body_operation;
        static const std::unordered_map<std::wstring, Type> two_body_operation;
        static const std::unordered_map<std::wstring, Type> fock_operation;
        static const std::map<Type,std::wstring> oper_to_string;
        static const std::map<Option,std::wstring> option_to_string;
        static const std::unordered_map<std::wstring, Option> string_to_option;

        Operator() = default;
        Operator(Operator const &) = default;
        Operator(Operator &&) = default;
        Operator& operator=(Operator const &) = default;
        Operator& operator=(Operator &&) = default;

        /**
         * Constructor
         *
         * @param operation  string of operation
         * @param option string of a list of option separated by ",", it will sort options so that they
         * always come in the same order
         */
        Operator(std::wstring operation, std::wstring option = L"");

        /// return options options_
        const std::vector<Option> option() const {
            return option_;
        }

        /// return operation oper_
        const Type &type() const {
            return type_;
        }

        /// set operation oper_
        void set_type(const Type &oper) {
            Operator::type_ = oper;
        }

        /// return string that correspond to oper_
        const std::wstring oper_string() const;

        /// return string that correspond to options_ wraped in []
        const std::wstring option_string() const;

        /// return true if have Option op
        bool has_option(Option op) const;

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
        bool operator==(Operator const & other) const;

        /// equality check by comparing oper_ and options_
        bool operator!=(Operator const & other) const{
            return !(*this==other);
        }

        /// comparison by comparing oper_ and options_
        bool operator<(const Operator& other) const;

    private:

        Type type_;
        std::vector<Option> option_;

    };

}

#endif  // SRC_MPQC_UTIL_EXPRESSION_OPERATOR_H_

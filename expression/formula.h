//
// Created by Chong Peng on 10/15/15.
//

#ifndef TILECLUSTERCHEM_FORMULA_H
#define TILECLUSTERCHEM_FORMULA_H

#include <vector>

#include "greek_to_english_name.h"
#include "orbital_index.h"
#include "operation.h"

using mpqc::OrbitalIndex;
using mpqc::Operation;

namespace mpqc{

    /**
     * \brief Formula class that represent quantum chemistry equations
     *
     * format for formula

        - Physical Notation
        <index1 index2|operation|index3 index4>[option1,option2]

        - Chemical Notation
        (index1 index2|operation|index3 index4)[option]

        see OrbitalIndex for description of index
        see Operation for description of operation and option

     */
    class Formula{
    public:

        /// Types of Notation
        enum class Notation {Chemical=0, Physical=1};

        Formula() = default;
        Formula(Formula const &) = default;
        Formula(Formula &&) = default;
        Formula& operator=(Formula const &) = default;
        Formula& operator=(Formula &&) = default;

        /**
         *  Constructor
         *  @param formula a string that has formula format
         */
        Formula(std::wstring formula);

        /// reconstruct the formula in string formula format
        std::wstring formula_string() const;

        /// return left_index
        const std::vector<OrbitalIndex> &left_index() const {
            return left_index_;
        }

        /// return left_index
        std::vector<OrbitalIndex> &left_index() {
            return left_index_;
        }

        /// return right_index
        const std::vector<OrbitalIndex> &right_index() const {
            return right_index_;
        }

        /// return right_index
        std::vector<OrbitalIndex> &right_index() {
            return right_index_;
        }

        /// set right_index
        void set_right_index(const std::vector<OrbitalIndex> &right_index) {
            Formula::right_index_ = right_index;
        }

        /// set left_index
        void set_left_index(const std::vector<OrbitalIndex> &left_index) {
            Formula::left_index_ = left_index;
        }

        /// set operation
        void set_operation(const Operation &operation) {
            Formula::operation_ = operation;
        }

        /// set notation
        void set_notation(const Notation &notation) {
            Formula::notation_ = notation;
        }

        /// return operation
        const Operation &operation() const {
            return operation_;
        }

        /// return operation
        Operation &operation(){
            return operation_;
        }

        /// return notation
        const Notation &notation() const {
            return notation_;
        }

        /// check if formula has index in left_index and right_index
        bool has_index(const OrbitalIndex& index) const;

        /// dimension of formula(2, 3 or 4)
        std::size_t rank() const;

        /// comparison by comparing operation, notation, left_index and right_index
        bool operator<(const Formula& other) const;

        /// check equality by comparing operation, left_index, right_index and notation
        bool operator==(const Formula& other) const;

        /// check equality by comparing operation, left_index, right_index and notation
        bool operator!=(const Formula& other) const {
            return !(*this==other);
        }

        /// convert to TA expression string format
        std::string to_ta_expression() const;

    private:

        /// parse the index on one side
        std::vector<OrbitalIndex> check_orbital_index(std::wstring index_array);

    private:

        Operation operation_;
        Notation  notation_;
        std::vector<OrbitalIndex> left_index_;
        std::vector<OrbitalIndex> right_index_;
    };
}


#endif //TILECLUSTERCHEM_FORMULA_H

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

    /* format for formula

        <index1 index2|operation|index3 index4>[option1,option2]

        (index1 index2|operation|index3 index4)[option]

     */
    class Formula{
    public:
        enum class Notation {Chemical, Physical};


        Formula() = default;
        Formula(Formula const &) = default;
        Formula(Formula &&) = default;
        Formula& operator=(Formula const &) = default;
        Formula& operator=(Formula &&) = default;

        Formula(std::wstring formula);

        std::wstring formula_string() const;

        const std::vector<OrbitalIndex> &left_index() const {
            return left_index_;
        }

        std::vector<OrbitalIndex> &left_index() {
            return left_index_;
        }

        const std::vector<OrbitalIndex> &right_index() const {
            return right_index_;
        }

        void set_right_index(const std::vector<OrbitalIndex> &right_index) {
            Formula::right_index_ = right_index;
        }

        void set_left_index(const std::vector<OrbitalIndex> &left_index) {
            Formula::left_index_ = left_index;
        }

        std::vector<OrbitalIndex> &right_index() {
            return right_index_;
        }

        const Operation &operation() const {
            return operation_;
        }

        const Notation &notation() const {
            return notation_;
        }

        std::size_t rank() const;

        bool operator<(const Formula& other) const;

        bool operator==(const Formula& other) const;

        bool operator!=(const Formula& other) const {
            return !(*this==other);
        }

        // convert to TA expression string format
        std::string to_ta_expression() const;

    private:

        std::vector<OrbitalIndex> check_orbital_index(std::wstring index_array);

    private:

        Operation operation_;
        Notation  notation_;
        std::vector<OrbitalIndex> left_index_;
        std::vector<OrbitalIndex> right_index_;
    };
}


#endif //TILECLUSTERCHEM_FORMULA_H

//
// Created by Chong Peng on 10/15/15.
//

#ifndef TILECLUSTERCHEM_FORMULA_H
#define TILECLUSTERCHEM_FORMULA_H

#include <vector>
#include <unordered_map>

#include "orbital_index.h"

using mpqc::OrbitalIndex;

namespace mpqc{

    /* format for formula

        <index1 index2|operation|index3 index4>

     */
    class Formula{
    public:
        enum class Operation{Overlap, Kinetic, Nuclear, Coulomb, cGTG, cGTGCoulomb, cGTG2};

        static const std::unordered_map<std::string, Operation> string_to_operation;

        Formula() = default;
        Formula(Formula const &) = default;
        Formula(Formula &&) = default;
        Formula& operator=(Formula const &) = default;
        Formula& operator=(Formula &&) = default;

        const std::string &formula() const {
            return formula_;
        }

        const std::vector<OrbitalIndex> &left_index() const {
            return left_index_;
        }

        const std::vector<OrbitalIndex> &right_index() const {
            return right_index_;
        }

        const Operation &operation() const {
            return operation_;
        }

        Formula(std::string formula);

    private:

        Operation check_operation(std::string oper);
        std::vector<OrbitalIndex> check_orbital_index(std::string index_array);

    private:

        std::string formula_;
        Operation operation_;
        std::vector<OrbitalIndex> left_index_;
        std::vector<OrbitalIndex> right_index_;
    };
}


#endif //TILECLUSTERCHEM_FORMULA_H

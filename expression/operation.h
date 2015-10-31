//
// Created by Chong Peng on 10/31/15.
//

#ifndef TILECLUSTERCHEM_OPERATION_H
#define TILECLUSTERCHEM_OPERATION_H

#include <unordered_map>
#include <string>

namespace mpqc{

    class Operation{
    public:
        enum class Operations{Overlap, Kinetic, Nuclear, Coulomb, cGTG, cGTGCoulomb, cGTG2};

        static const std::unordered_map<std::wstring, Operations> one_body_operation;
        static const std::unordered_map<std::wstring, Operations> two_body_operation;


        Operation() = default;
        Operation(Operation const &) = default;
        Operation(Operation &&) = default;
        Operation& operator=(Operation const &) = default;
        Operation& operator=(Operation &&) = default;

        Operation(std::wstring operation);

        const Operations &get_operation() const {
            return operation_;
        }

        bool is_onebody() const;
        bool is_twobody() const;


    private:

        Operations operation_;

    };


}



#endif //TILECLUSTERCHEM_OPERATION_H

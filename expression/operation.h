//
// Created by Chong Peng on 10/31/15.
//

#ifndef TILECLUSTERCHEM_OPERATION_H
#define TILECLUSTERCHEM_OPERATION_H

#include <unordered_map>
#include <vector>
#include <string>

namespace mpqc{

    class Operation{
    public:
        enum class Operations{Overlap, Kinetic, Nuclear, Coulomb, cGTG, cGTGCoulomb, cGTG2};

        enum class Options{DensityFitting};

        static const std::unordered_map<std::wstring, Operations> one_body_operation;
        static const std::unordered_map<std::wstring, Operations> two_body_operation;
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

        bool has_option(Options op) const;

        bool is_onebody() const;
        bool is_twobody() const;
        bool is_r12() const;

        bool operator==(Operation const & other) const;

    private:

        Operations operation_;
        std::vector<Options> options_;

    };


}



#endif //TILECLUSTERCHEM_OPERATION_H

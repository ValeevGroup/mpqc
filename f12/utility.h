//
// Created by Chong Peng on 10/21/15.
//

#ifndef TILECLUSTERCHEM_UTILITY_H
#define TILECLUSTERCHEM_UTILITY_H

#include <string>
#include "../include/eigen.h"
#include <vector>

namespace mpqc{
namespace f12{

    double basis_to_f12exponent(const std::string& basis_name);

    std::vector<std::pair<double,double>> stg_ng_fit(std::size_t n, double zeta);

    std::vector<std::pair<double,double>> gtg_params_squared(const std::vector<std::pair<double,double>>& pragmas);

    // Gaussian Type Geminal parameters
    struct GTGParams{

        double exponent;
        int n_fit;

        GTGParams(double zeta, int n = 6) : exponent(zeta), n_fit(n) { }
        GTGParams(std::string basis_name, int n = 6) : exponent(basis_to_f12exponent(basis_name)), n_fit(n) {}

        std::vector<std::pair<double,double>> compute()
        {
            return stg_ng_fit(n_fit,exponent);
        };
    };


}
}
#endif //TILECLUSTERCHEM_UTILITY_H

//
// Created by Chong Peng on 10/21/15.
//

#ifndef TILECLUSTERCHEM_UTILITY_H
#define TILECLUSTERCHEM_UTILITY_H

#include <string>
#include <Eigen/Dense>
#include <vector>

namespace mpqc{
namespace f12{

    double basis_to_f12exponent(const std::string& basis_name);

    std::vector<std::pair<double,double>> stg_ng_fit(int n, double zeta);


}
}
#endif //TILECLUSTERCHEM_UTILITY_H

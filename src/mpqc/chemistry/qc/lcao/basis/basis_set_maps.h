
#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_BASIS_BASIS_SET_MAPS_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_BASIS_BASIS_SET_MAPS_H_

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_BASIS_BASIS_SET_MAPS_H_
#include <map>
#include <string>

namespace mpqc {
namespace lcao {
static std::map<std::string, unsigned int> elements{
    {"H", 1},   {"He", 2},  {"Li", 3},  {"Be", 4},  {"B", 5},   {"C", 6},
    {"N", 7},   {"O", 8},   {"F", 9},   {"Ne", 10}, {"Na", 11}, {"Mg", 12},
    {"Al", 13}, {"Si", 14}, {"P", 15},  {"S", 16},  {"Cl", 17}, {"Ar", 18},
    {"K", 19},  {"Ca", 20}, {"Sc", 21}, {"Ti", 22}, {"V", 23},  {"Cr", 24},
    {"Mn", 25}, {"Fe", 26}, {"Co", 27}, {"Ni", 28}, {"Cu", 29}, {"Zn", 30},
    {"Ga", 31}, {"Ge", 32}, {"As", 33}, {"Se", 34}, {"Br", 35}, {"Kr", 36},
};

static std::map<std::string, unsigned int> ang_mo_map{
    {"S", 0}, {"P", 1}, {"D", 2}, {"F", 3}, {"G", 4}, {"H", 5}, {"I", 6}};

}  // namespace lcao
}  // namespace mpqc

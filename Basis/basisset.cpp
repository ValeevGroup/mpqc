#include "basisset.h"
#include <fstream>
#include <sstream>

void BasisSet::read_basis(std::string const &s) const {
    std::ifstream basis_file(s);
    std::string line;
    while (std::getline(basis_file, line)) {
        std::stringstream ss(line);
        std::string elem;
        ss >> elem;
    }
}

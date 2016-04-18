/*
 * Example test to ensure that I can compile with libint basis headers
 */

#include <libint2.hpp>

int main(int argc, char *argv[]) {
    libint2::initialize();

    libint2::Atom h1{1, 0.0, 0.0, 0.0};
    libint2::Atom h2{1, 0.0, 0.0, 1.0};
    std::vector<libint2::Atom> mol = {std::move(h1), std::move(h2)};

    libint2::BasisSet sto3g("STO-3G", mol);

    std::cout << "STO-3G Basis set for H2 = \n" ;
    for(auto const &shell : sto3g){
        std::cout << "\t" << shell << std::endl;
    }
        
    libint2::finalize();

    return 0;
}

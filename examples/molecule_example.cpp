#include "../molecule/cluster.h"
#include "../molecule/atom.h"

#include "../molecule/molecule.h"

#include <iostream>

using namespace mpqc;
using namespace molecule;

int main(int argc, char **argv) {
    Cluster c12(Atom({0, 0, 0}, 1, 1), Atom({0, 0, 1}, 1, 1));
    Cluster c34(Atom({0, 0, 2}, 1, 1), Atom({0, 0, 3}, 1, 1));

    return 0;
}

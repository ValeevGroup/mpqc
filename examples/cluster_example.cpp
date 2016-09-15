#include "../molecule/cluster.h"
#include "../molecule/atom.h"

using namespace mpqc;
using namespace molecule;

int main(int argc, char **argv) {
    Atom a1({0, 0, 0}, 1, 1);
    Atom a2({0, 0, 1}, 1, 1);
    Atom a3({0, 0, 2}, 1, 1);
    Atom a4({0, 0, 3}, 1, 1);
    Atom a5({0, 0, 4}, 1, 1);
    Atom a6({0, 0, 5}, 1, 1);

    std::cout << "Atom 1 = " << a1 << std::endl;
    std::cout << "Atom 2 = " << a2 << std::endl;
    std::cout << "Atom 3 = " << a3 << std::endl;
    std::cout << "Atom 4 = " << a4 << std::endl;
    std::cout << "Atom 5 = " << a5 << std::endl;
    std::cout << "Atom 6 = " << a6 << std::endl;

    // Add clusterable construction
    Cluster c12;
    c12.add_clusterable(std::move(a1));
    c12.add_clusterable(std::move(a2));

    // Varadic template construction
    Cluster c34(std::move(a3), std::move(a4));

    // Construct with a vector
    std::vector<Clusterable> a56;
    a56.emplace_back(std::move(a5));
    a56.emplace_back(std::move(a6));
    Cluster c56(std::move(a56));

    // We can print clusters
    std::cout << "Cluster 12 = " << c12 << std::endl;
    std::cout << "Cluster 34 = " << c34 << std::endl;
    std::cout << "Cluster 56 = " << c56 << std::endl;

    // We can compute centers for clusters
    c12.update_center();
    c34.update_center();
    c56.update_center();

    std::cout << "Cluster 12 center = " << c12.center().transpose()
              << std::endl;
    std::cout << "Cluster 34 center = " << c34.center().transpose()
              << std::endl;
    std::cout << "Cluster 56 center = " << c56.center().transpose()
              << std::endl;

    // We can nest clusters
    Cluster nc1234(std::move(c12), std::move(c34));
    std::cout << "\nNested Clusters 1 2 3 4 = " << nc1234 << std::endl;

    // Finally get its center
    nc1234.update_center();
    std::cout << "Nested Clusters 1 2 3 4 center = "
              << nc1234.center().transpose() << std::endl;

    return 0;
}

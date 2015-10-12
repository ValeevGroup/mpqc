#include "catch.hpp"

#include "../include/eigen.h"
#include "../common/typedefs.h"

#include "../clustering/kmeans.h"

#include <iostream>

using namespace mpqc;

namespace clustering_test {

// Dummy class cluster for testing k-means.
struct Cluster {
    Vec3D center;
    std::vector<Vec3D> elements;

    auto begin() const -> decltype(elements.begin()) {
        return elements.begin();
    }

    auto end() const -> decltype(elements.end()) { return elements.end(); }
};

void remove_clusterables(Cluster &cluster) { cluster.elements.clear(); }

void attach_clusterable(Cluster &cluster, Vec3D const &vec) {
    cluster.elements.emplace_back(vec);
}

void update_center(Cluster &cluster) {
    Vec3D new_center = {0.0, 0.0, 0.0};
    for (auto const &elem : cluster.elements) {
        new_center += elem;
    }
    new_center /= double(cluster.elements.size());

    cluster.center = new_center;
}

void set_center(Cluster &cluster, Vec3D const &vec) { cluster.center = vec; }

Vec3D const &center(Cluster const &cluster) { return cluster.center; }

std::ostream &operator<<(std::ostream &os, Cluster const &cluster) {
    os << "Cluster(" << cluster.elements.size() << "):\n";
    os << "\tcenter: " << cluster.center.transpose() << "\n";
    for (auto const &vec : cluster.elements) {
        os << "\tvec:" << vec.transpose() << std::endl;
    }

    return os;
}

} // namespace clustering test

namespace Eigen {
Vec3D center(Vec3D const &vec) { return vec; }
}


TEST_CASE("k-means can cluster clusterables", "[k-means, clustering]") {

    namespace ctest = clustering_test;
    using Cluster = ctest::Cluster;

    auto default_kmeans = clustering::Kmeans();
    auto vectors = std::vector<Vec3D>(
          {Vec3D{0.0, 0.0, 0.0}, Vec3D{1.0, 0.0, 0.0}, Vec3D{0.0, 1.0, 0.0},
           Vec3D{2.0, 2.0, 2.0}, Vec3D{3.0, 2.0, 2.0}, Vec3D{2.0, 3.0, 2.0}});

    SECTION("k-means can form different number of clusters") {
        REQUIRE_THROWS(default_kmeans.cluster<Cluster>(vectors, -1));
        REQUIRE_THROWS(default_kmeans.cluster<Cluster>(vectors, 0));
        auto one_clusters = default_kmeans.cluster<Cluster>(vectors, 1);
        REQUIRE(one_clusters.size() == 1ul);
        REQUIRE(clustering::kmeans_objective(one_clusters)
                == Approx(20.666666667));

        auto two_clusters = default_kmeans.cluster<Cluster>(vectors, 2);
        REQUIRE(two_clusters.size() == 2ul);
        REQUIRE(clustering::kmeans_objective(two_clusters)
                == Approx(2.666666667));

        auto three_clusters = default_kmeans.cluster<Cluster>(vectors, 3);
        REQUIRE(three_clusters.size() == 3ul);
        REQUIRE(clustering::kmeans_objective(three_clusters)
                == Approx(1.8333333333));

        auto six_clusters = default_kmeans.cluster<Cluster>(vectors, 6);
        REQUIRE(six_clusters.size() == 6ul);
        REQUIRE(clustering::kmeans_objective(six_clusters) == 0.0);

        REQUIRE_THROWS(default_kmeans.cluster<Cluster>(vectors, 7));
    }
}

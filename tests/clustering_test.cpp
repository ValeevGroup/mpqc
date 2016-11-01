#include "catch.hpp"

#include "mpqc/math/external/eigen/eigen.h"


#include "../clustering/kmeans.h"

using namespace mpqc;

namespace clustering_test {

// Dummy class cluster for testing k-means.
struct Cluster {
    Vector3d center;
    std::vector<Vector3d> elements;

    auto begin() const -> decltype(elements.begin()) {
        return elements.begin();
    }

    auto end() const -> decltype(elements.end()) { return elements.end(); }
};

void remove_clusterables(Cluster &cluster) { cluster.elements.clear(); }

void attach_clusterable(Cluster &cluster, Vector3d const &vec) {
    cluster.elements.emplace_back(vec);
}

void update_center(Cluster &cluster) {
    Vector3d new_center = {0.0, 0.0, 0.0};
    for (auto const &elem : cluster.elements) {
        new_center += elem;
    }
    new_center /= double(cluster.elements.size());

    cluster.center = new_center;
}

void set_center(Cluster &cluster, Vector3d const &vec) { cluster.center = vec; }

Vector3d const &center(Cluster const &cluster) { return cluster.center; }

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
Vector3d center(Vector3d const &vec) { return vec; }
}


TEST_CASE("k-means can cluster clusterables", "[k-means, clustering]") {

    namespace ctest = clustering_test;
    using Cluster = ctest::Cluster;

    auto default_kmeans = clustering::Kmeans();
    auto vectors = std::vector<Vector3d>(
          {Vector3d{0.0, 0.0, 0.0}, Vector3d{1.0, 0.0, 0.0}, Vector3d{0.0, 1.0, 0.0},
           Vector3d{2.0, 2.0, 2.0}, Vector3d{3.0, 2.0, 2.0}, Vector3d{2.0, 3.0, 2.0}});

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

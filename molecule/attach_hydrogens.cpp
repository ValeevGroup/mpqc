#include "attach_hydrogens.h"
#include "cluster.h"
#include "../include/tbb.h"
#include <algorithm>
#include <iostream>
#include <unordered_map>

namespace tcc {
namespace molecule {
namespace clustering {

// Put every clusterable into its own cluster.
std::vector<Cluster>
individual_clustering(std::vector<Clusterable> const &c_ables) {
    std::vector<Cluster> cs;
    cs.reserve(c_ables.size());
    for (auto const &c_able : c_ables) {
        Cluster c;
        c.add_clusterable(c_able);
        cs.push_back(std::move(c));
    }
    return cs;
}

// A map to keep track of which heavy the hydrogens belog to.
// Similar to unordered multimap, but I used a vector instead.
template <typename It>
std::vector<std::pair<It, std::vector<It>>>
init_ownership_map(std::vector<Clusterable> const &clusterables,
                   const It first_heavy) {
    auto end = clusterables.end();
    std::vector<std::pair<It, std::vector<It>>> h_owners;
    // need decltype due to constness issues.
    h_owners.reserve(std::distance<decltype(end)>(first_heavy, end));

    // Make one entry for each heavy.
    for (auto it = first_heavy; it != end; ++it) {
        h_owners.emplace_back(std::make_pair(it, std::vector<It>{}));
    }

    return h_owners;
}

template <typename It>
void assign_hydrogens_to_heavies(
    std::vector<Clusterable> const &cs, It first_heavy,
    std::vector<std::pair<It, std::vector<It>>> &h_owners) {

    auto end = cs.end();
    using range = tbb::blocked_range<It>;

    /* TODO find way to make parallel loop deterministic 
    // Lambda finds the closes heavy to each hydrogen then calculates the
    // heavies position in the owner vector and finally adds the hydrogen iter
    // to the list of hydrogens owned by that heavy.
    tbb::spin_mutex storage_lock;
    auto h_attacher = [&](range const &r) {
        for (auto h_it = r.begin(); h_it != r.end(); ++h_it) {

            auto nearest_heavy
                = std::min_element(first_heavy, end, [&](Clusterable const &a,
                                                         Clusterable const &b) {
                    double a_dist = (a.center() - h_it->center()).norm();
                    double b_dist = (b.center() - h_it->center()).norm();
                    return a_dist < b_dist;
                });

            auto heavy_ordinal = std::distance(first_heavy, nearest_heavy);

            tbb::spin_mutex::scoped_lock lock(storage_lock);
            h_owners[heavy_ordinal].second.push_back(h_it);
        }
    };

    tbb::parallel_for(range(cs.begin(), first_heavy), h_attacher);
    */ 

    // Use serial version until above is fixed
    for (auto it = cs.begin(); it != first_heavy; ++it) {
        auto nearest_heavy = std::min_element(
            first_heavy, end, [&it](Clusterable const &a, Clusterable const &b) {
                double a_dist = (a.center() - it->center()).norm();
                double b_dist = (b.center() - it->center()).norm();
                return a_dist < b_dist;
            });

        auto heavy_ordinal = std::distance(first_heavy, nearest_heavy);
        h_owners[heavy_ordinal].second.push_back(it);
    }
}

template <typename It>
std::vector<Cluster>
attach_to_owners(std::vector<std::pair<It, std::vector<It>>> const &h_owners) {
    std::vector<Cluster> cs;
    cs.reserve(h_owners.size());

    // For each heavy make a cluster and add the heavy + all hydrogens to that
    // cluster, finally add cluster to vector of clusters.
    for (auto const &h_pair : h_owners) {
        Cluster c;
        c.add_clusterable(*h_pair.first);
        for (auto h_it : h_pair.second) { // just pointers so copy is ok.
            c.add_clusterable(*h_it);
        }
        cs.push_back(std::move(c));
    }

    for (auto &c : cs) {
        c.compute_com();
    }

    return cs;
}


std::vector<Cluster> attach_hydrogens::
operator()(std::vector<Clusterable> clusterables) {
    const auto end = clusterables.end();

    // TODO ensure test this loop to ensure it doesn't give different orders on 
    // different nodes
    // 
    // tbb::parallel_sort(
    //     clusterables.begin(), end,
    //     [](Clusterable &a, Clusterable &b) { return a.charge() <
    //     b.charge();
    //     });

    // serial stable sort to divide the range into hydrogens and heavies
    std::stable_sort(clusterables.begin(), end,
                     [](Clusterable const &a, Clusterable const &b) {
        return a.charge() < b.charge();
    });

    // find the location of the first heavy atom.
    decltype(clusterables.cbegin()) first_heavy = std::upper_bound(
        clusterables.begin(), end, 1,
        [](int val, Clusterable const &c) { return val < c.charge(); });

    // Check for no heavies or no hydrogens.
    if (first_heavy == end || first_heavy == clusterables.begin()) {
        return individual_clustering(clusterables);
    }

    // Figure out which hydrogens go with which heavies.
    auto h_owners = init_ownership_map(clusterables, first_heavy);
    assign_hydrogens_to_heavies(clusterables, first_heavy, h_owners);

    // Create and return vector of clusters with hydrogens attached.
    return attach_to_owners(h_owners);
}

} // namespace clustering
} // namespace molecule
} // namespace tcc

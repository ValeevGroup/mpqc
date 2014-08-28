#ifndef COMMON_H
#define COMMON_H

#include "cluster_concept.h"
#include "molecule_fwd.h"

/**
 * center_of_mass returns an position vector which is the center of mass
 * of the input types with centers.
 */
position_t center_of_mass(const std::vector<Clusterable> cs, double mass);

inline double diff_squaredNorm(const position_t a, const position_t b) {
    double out = 0.0;
    for (unsigned int i = 0; i < 3; ++i) {
        double temp = a[i] - b[i];
        out += temp * temp;
    }
    return out;
}

#endif // COMMON_H

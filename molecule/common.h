#ifndef COMMON_H
#define COMMON_H

#include "cluster_concept.h"
#include "molecule_fwd.h"

/**
 * center_of_mass returns an position vector which is the center of mass
 * of the input types with centers.
 */
position_t center_of_mass(const std::vector<Clusterable> cs, double mass);


#endif // COMMON_H

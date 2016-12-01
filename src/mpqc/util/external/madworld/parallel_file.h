//
// Created by Chong Peng on 3/24/16.
//

#ifndef MPQC4_SRC_MPQC_UTIL_EXTERNAL_MADWORLD_PARALLEL_FILE_H_
#define MPQC4_SRC_MPQC_UTIL_EXTERNAL_MADWORLD_PARALLEL_FILE_H_

#include <madness/world/world.h>
#include <madness/world/worldgop.h>

namespace mpqc {
namespace utility {

void parallel_read_file(madness::World &world, const std::string &filename,
                        char *&buffer);

void parallel_read_file(madness::World &world, const std::string &filename,
                        std::stringstream &output);

}  // end of namespace utility
}  // end of namespace mpqc

#endif  // MPQC4_SRC_MPQC_UTIL_EXTERNAL_MADWORLD_PARALLEL_FILE_H_

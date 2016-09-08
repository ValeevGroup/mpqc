//
// Created by Chong Peng on 8/26/16.
//

#ifndef MPQC_TILE_CONVERT_H
#define MPQC_TILE_CONVERT_H

#include <../include/tiledarray.h>

namespace mpqc{
namespace ta_routines{


inline TA::TensorD TATensorDPassThrough(TA::TensorD &&ten) { return std::move(ten); }


}
}



#endif //MPQC_TILE_CONVERT_H

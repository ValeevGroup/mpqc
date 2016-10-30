//
// Created by Chong Peng on 8/26/16.
//

#ifndef MPQC_TILE_CONVERT_H
#define MPQC_TILE_CONVERT_H

#include <../include/tiledarray.h>
#include <../common/namespaces.h>

namespace mpqc{
namespace ta_routines{

struct TensorDPassThrough {
  TA::TensorD operator()(TA::TensorD &&ten) const {
    return std::move(ten);
  }
};

struct TensorZPassThrough {
  TA::TensorZ operator()(TA::TensorZ &&ten) const {
    return std::move(ten);
    }
};

}
}



#endif //MPQC_TILE_CONVERT_H

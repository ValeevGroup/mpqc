//
// Created by Chong Peng on 8/26/16.
//

#include <tile_convert.h>

TA::TensorD ta_tensor_pass_through(TA::TensorD &&ten)
{
  return std::move(ten);
}

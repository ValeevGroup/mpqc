/*
 * print.cpp
 *
 *  Created on: Nov 1, 2016
 *      Author: evaleev
 */


#include "mpqc/util/misc/print.h"

namespace mpqc{
namespace util{

// print progress
void print_progress(std::size_t lowprogress, std::size_t upprogress,
                    std::vector<std::size_t>& progress_points) {
  std::size_t size = progress_points.size();
  std::size_t total = progress_points[size - 1];
  for (auto i = lowprogress; i < upprogress; i++) {
    if (std::find(progress_points.begin(), progress_points.end(), i) !=
        progress_points.end()) {
      int percent = 10 * (double(i) / total);
      percent *= 10;
      ExEnv::out0() << percent << "% done." << std::endl;
    }
  }
}






} // end of namespace util
} // end of namespace mpqc




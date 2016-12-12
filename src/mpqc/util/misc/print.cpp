/*
 * print.cpp
 *
 *  Created on: Nov 1, 2016
 *      Author: evaleev
 */

#include <cmath>
#include <iostream>
#include <vector>

#include <algorithm>

#include "mpqc/util/misc/print.h"

// print progress
void print_progress(std::size_t lowprogress, std::size_t upprogress,
                    std::vector<std::size_t>& progress_points) {
  std::size_t total = progress_points.size();
  if (total < 10) {
    std::size_t divide = total;
    std::size_t percent = 10 * double(upprogress) / divide;
    std::cout << 10*percent << "% done." << std::endl;
  } else {
    for (int i = lowprogress; i < upprogress; i++) {
      if (std::find(progress_points.begin(), progress_points.end(), i) !=
          progress_points.end()) {
        int percent = 100 * double(i) / total;
        std::cout << percent << "% done." << std::endl;
      }
    }
  }
}

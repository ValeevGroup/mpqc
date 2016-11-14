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
void print_progress(int lowprogress, int upprogress, int total) {
  int divide = 10;
  if (total < 10) {
    divide = total;
    int percent = 100 * double(upprogress) / divide;
    std::cout << percent << "% done." << std::endl;
  } else {
    int increase = std::round(double(total) / divide);
    std::vector<int> progress_points;
    for (int i = 0; i < total; i += increase) {
      progress_points.push_back(i);
    }

    for (int i = lowprogress; i < upprogress; i++) {
      if (std::find(progress_points.begin(), progress_points.end(), i) !=
          progress_points.end()) {
        int percent = 100 * double(i) / total;
        std::cout << percent << "% done." << std::endl;
      }
    }
  }
}

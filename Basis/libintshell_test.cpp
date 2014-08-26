#include "libintshell.h"
#include <iostream>

int main(int argc, char *argv[]) {
  auto s = libint2::Shell{
          {5.033151300, 1.169596100, 0.380389000},
                 {{0, false, {-0.09996723, 0.39951283, 0.70011547}}, {
                     1, false, {0.15591627, 0.60768372, 0.39195739}}},
                 {{0.0, 0.0, 0.0}}};
  std::cout << "Shell s = \n" << s << std::endl;
  return 0;
}

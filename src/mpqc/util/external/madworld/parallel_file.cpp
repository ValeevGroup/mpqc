//
// Created by Chong Peng on 3/24/16.
//

#include "mpqc/util/external/madworld/parallel_file.h"

#include <fstream>

namespace mpqc {
namespace utility {

void parallel_read_file(madness::World &world,
                        const std::string &filename, std::unique_ptr<char[]>& buffer) {
  int size;
  std::string contents;
  if (world.rank() == 0) {
    std::ifstream input_file(filename, std::ifstream::in);

    if (input_file.fail()) {
      std::ostringstream oss;
      oss << "could not open file \"" << filename << "\"";
      throw std::invalid_argument(oss.str().c_str());
    }

    contents = std::string((std::istreambuf_iterator<char>(input_file)),
                           std::istreambuf_iterator<char>());

    input_file.close();
    size = contents.size() + 1;
  }

  world.gop.broadcast(size, 0);
  buffer = std::make_unique<char[]>(size);

  if (world.rank() == 0) {
    strcpy(buffer.get(), contents.c_str());
  }

  world.gop.broadcast(buffer.get(), size, 0);
}

void parallel_read_file(madness::World &world,
                        const std::string &filename,
                        std::stringstream &output) {
  std::unique_ptr<char[]> buffer;
  parallel_read_file(world, filename, buffer);
  output << buffer.get();
}

}  // namespace utility
}  // namespace mpqc

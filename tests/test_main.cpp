#define CATCH_CONFIG_RUNNER

#include <clocale>
#include <madness/world/world.h>
#include "catch.hpp"

int main( int argc, char* argv[] )
{
  Catch::Session session;
  // global setup...
  std::setlocale(LC_ALL,"en_US.UTF-8");

  auto &world = madness::initialize(argc, argv);

  int result = session.run( argc, argv );

  // global clean-up...

  return result;
}

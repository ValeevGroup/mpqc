#define CATCH_CONFIG_RUNNER

#include <clocale>
#include <madness/world/world.h>
#include "catch.hpp"

madness::World* default_world = nullptr;

int main( int argc, char* argv[] )
{
  Catch::Session session;
  // global setup...
  std::setlocale(LC_ALL,"en_US.UTF-8");

  default_world = &madness::initialize(argc, argv);

  int result = session.run( argc, argv );

  // global clean-up...

  return result;
}

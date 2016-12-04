#define CATCH_CONFIG_RUNNER

#include <clocale>
#include <tiledarray.h>
#include "catch.hpp"

int main( int argc, char* argv[] )
{
  Catch::Session session;
  // global setup...
  std::setlocale(LC_ALL,"en_US.UTF-8");

  auto& world = madness::initialize(argc, argv);
  TiledArray::set_default_world(world);

  int result = session.run( argc, argv );

  // global clean-up...
  madness::finalize();

  return result;
}

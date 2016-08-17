#define CATCH_CONFIG_RUNNER

#include <clocale>
#include "catch.hpp"

int main( int argc, char* const argv[] )
{
  Catch::Session session;
  // global setup...
  std::setlocale(LC_ALL,"en_US.UTF-8");

  int result = session.run( argc, argv );

  // global clean-up...

  return result;
}

#ifdef __GNUG__
#pragma implementation
#endif

#include <stdlib.h>
#include <iostream>

#include <util/misc/exenv.h>
#include <chemistry/qc/oint3/build.h>

using namespace std;

BuildIntV3::BuildIntV3()
{
}

BuildIntV3::~BuildIntV3()
{
}

int
BuildIntV3::impossible_integral()
{
  ExEnv::err()
      << "oint3/build.cc: tried to build a impossible integral" << endl;
  abort();
  return(0);
}

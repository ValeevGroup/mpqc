#include <stdlib.h>
#include <iostream>

#include <util/misc/exenv.h>
#include <chemistry/qc/oint3/build.h>

using namespace std;
using namespace sc;

BuildIntV3::BuildIntV3()
{
}

BuildIntV3::~BuildIntV3()
{
}

int
BuildIntV3::impossible_integral()
{
  ExEnv::errn()
      << "oint3/build.cc: tried to build a impossible integral" << endl;
  abort();
  return(0);
}

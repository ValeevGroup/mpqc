#ifdef __GNUG__
#pragma implementation
#endif

#include <chemistry/qc/oint3/build.h>

int
BuildIntV3::impossible_integral()
{
  fprintf(stderr,"oint3/build.cc: tried to build a impossible integral\n");
  abort();
  return(0);
}


#include <stdio.h>
#include <stdlib.h>

//
// this is called by some GNU generated files.  only make it if we're not
// using libg++
//

#ifndef __GNUC__
extern "C" void lib_error_handler(const char *kind, const char* msg)
{
  fprintf(stderr,"%s\n",kind);
  fprintf(stderr,"Error: %s\n",msg);
  exit(1);
}
#endif


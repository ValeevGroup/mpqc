
#ifndef _util_misc_string_h
#define _util_misc_string_h

#include <string.h>
#include <stdlib.h>

namespace sc {

inline char * 
strdup (const char *string)
{
  return string ? strcpy ((char *) malloc (strlen (string) + 1), string) : 0;
}

}

#endif


#ifndef _util_misc_newstring_h
#define _util_misc_newstring_h

#include <string.h>

#ifdef __cplusplus

inline static char *
new_string(const char* s)
{
  if (!s) return 0;

  char *ret = new char[strlen(s)+1];
  strcpy(ret,s);
  return ret;
}

#endif
#endif


#include <stdio.h>
#include <assert.h>

void check_alloc (addr, str)
void *addr;
char *str;
{
  if (addr != 0) return;
  printf ("Allocation for %s failed\n", str);
#ifndef I860
  assert (0);
#else
  exit(0);
#endif
  }


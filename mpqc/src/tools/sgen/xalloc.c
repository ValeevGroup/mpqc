
#include <stdio.h>
#include <stdlib.h>

static char rcsid[]="$Id$";

char *
xmalloc (size)
     unsigned int size;
{
  char *result = malloc (size);
  if (result == 0) {
    fprintf (stderr,"virtual memory exhausted");
    exit(1);
    }
  return result;
}


char *
xrealloc (ptr, size)
     char *ptr;
     unsigned int size;
{
  char *result = realloc (ptr, size);
  if (result == 0) {
    fprintf (stderr,"virtual memory exhausted");
    exit(1);
    }
  return result;
}


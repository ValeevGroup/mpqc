
static char rcsid[] = "$Id$";
#include <stdio.h>
#include "sgen.h"

void
sgen_chkmalloc(buf)
void *buf;
{
  if (!buf) {
    fprintf(stderr,"memory allocation failed in sgen\n");
    exit(1);
    }
  }

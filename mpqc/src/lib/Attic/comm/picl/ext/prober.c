
#include <comm/picl/picl.h>
#include "piclext.h"

/* This routines sees if any messages have been left lying around.
 * If anything is found, then a message is printed.
 */

void
picl_prober()
{
  int bytes,type,source;

  if (!probe0(-1)) return;

  recvinfo0(&bytes,&type,&source);

  printf("On node %3d found a message of type %3d, bytes %3d, from %3d\n",
         mynode0(),type,bytes,source);
  }

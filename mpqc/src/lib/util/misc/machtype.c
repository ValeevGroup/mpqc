
#include <util/misc/machtype.gbl>
#include <util/misc/machtype.lcl>

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

GLOBAL_FUNCTION char *
machine_type()
{
  char *r = "unknown";
#ifdef HOST_ARCH
  r = HOST_ARCH;
#endif
  return r;
  }

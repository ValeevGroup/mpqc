comm/picl/high/piclhigh.a
#if defined(USE_PICLPROC)
#  include <comm/picl/proc/LIBS.h>
#elif defined(USE_PICLMSG)
#  include <comm/picl/msg/LIBS.h>
#elif defined(USE_PICLSHM)
#  include <comm/picl/shm/LIBS.h>
#elif defined(USE_PICLPARAGON)
#  include <comm/picl/paragon/LIBS.h>
#elif defined(USE_PICLN6400)
#  include <comm/picl/n6400/LIBS.h>
#elif defined(USE_PICLPVM)
#  include <comm/picl/pvm/LIBS.h>
#elif defined(USE_PICLNODE)
#  include <comm/picl/piclnode/LIBS.h>
#elif defined(PICLLOW)
PICLLOW
#else
#  if defined(PARAGON)
#    include <comm/picl/paragon/LIBS.h>
#  elif defined(NCUBE_V2)
#    include <comm/picl/n6400/LIBS.h>
#  elif defined(NCUBE_V3)
#    include <comm/picl/n6400/LIBS.h>
#  else
#    include <comm/picl/proc/LIBS.h>
#  endif
#endif

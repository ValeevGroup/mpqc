
#include "machtype.gbl"
#include "machtype.lcl"

GLOBAL_FUNCTION char *
machine_type()
{
  char *r = "unknown";
#if defined(SGI)
  r = "SGI";
#elif defined(NCUBE_V2)
  r = "NCUBE_V2";
#elif defined(NCUBE_V3)
  r = "NCUBE_V3";
#elif defined(SUN4)
  r = "SUN sparc";
#elif defined(SUN)
  r = "SUN";
#elif defined(L486)
  r = "Intel x86/Linux";
#elif defined(PARAGON)
  r = "PARAGON (i860/XP)";
#elif defined(I860)
  r = "iPSC/860 (i860/XR)";
#elif defined(RS6000)
  r = "RS6000";
#endif
  return r;
  }

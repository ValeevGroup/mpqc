#ifndef _util_misc_libmisc_h
#define _util_misc_libmisc_h

#include <util/misc/memory.h>
#include <util/misc/timer.h>
#include <util/misc/ieee.h>

#ifdef __cplusplus
extern "C" {
#endif

void check_alloc(void*,char*);
int mtype_get();

#ifdef __cplusplus
}
#endif

#endif /* _util_misc_libmisc_h */

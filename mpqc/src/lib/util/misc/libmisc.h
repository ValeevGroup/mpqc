#ifndef _util_misc_libmisc_h
#define _util_misc_libmisc_h

#define UTIL_ASSERT

#ifdef DEBUG_ON_ASSERT_FAIL
#undef UTIL_ASSERT
#define OLD_ASSERT
#endif

#include <util/misc/assert.h>
#include <util/misc/memory.h>
#include <util/misc/timer.h>
#include <util/misc/machtype.gbl>
#include <util/misc/ieee.h>

#ifdef __cplusplus
extern "C" {
#endif

void debug_start(char*);
void check_alloc(void*,char*);
int mtype_get();

#ifndef __GNUC__
void lib_error_handler(const char*, const char *msg);
#endif

#ifdef __cplusplus
}
#endif

#endif /* _util_misc_libmisc_h */

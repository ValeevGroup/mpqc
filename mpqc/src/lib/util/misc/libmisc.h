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

void check_alloc(void*,char*);
int mtype_get();

#endif /* _util_misc_libmisc_h */

/*	assert.h 1.7 88/02/07 SMI; from UCB 4.2 85/01/21	*/
#if defined(UTIL_ASSERT)

#include <util/misc/assert.gbl>

#define _assert(ex) {if(!(ex)){(void)util_assert(__FILE__, __LINE__);}}
#define assert(ex) _assert(ex)


#elif defined(OLD_ASSERT)


#ifdef NDBX
# ifndef NDEBUG
# define _assert(ex)	{if (!(ex)){(void)fprintf\
(stdout,"Assertion failed: file \"%s\", line %d\n", __FILE__, __LINE__);exit(23);}}
# define assert(ex)	_assert(ex)
# else
# define _assert(ex)
# define assert(ex)
# endif

#else

# ifndef NDEBUG
# define _assert(ex)	{if (!(ex)){(void)fprintf\
(stderr,"Assertion failed: file \"%s\", line %d\n", __FILE__, __LINE__);debug_start("assertion failure");}}
# define assert(ex)	_assert(ex)
# else
# define _assert(ex)
# define assert(ex)
# endif
#endif


#else


#include "/usr/include/assert.h"


#endif


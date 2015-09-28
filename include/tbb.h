#pragma once
#ifndef TCC_INCLUDE_TBB_H
#define TCC_INCLUDE_TBB_H

// Keep this as a way to deal with thread local integral engines. 
// It is header only and doesn't depend on linking the library I think.
#include <tbb/enumerable_thread_specific.h>
#include <tbb/spin_mutex.h>

#endif // TCC_INCLUDE_TBB_H

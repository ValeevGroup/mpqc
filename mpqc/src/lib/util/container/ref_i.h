
#ifdef INLINE_FUNCTIONS
#define INLINE inline
#else
#define INLINE
#endif

INLINE RefCount::RefCount():_reference_count_(0) {}
INLINE RefCount::RefCount(const RefCount&):_reference_count_(0) {}
INLINE RefCount& RefCount::operator=(const RefCount&) { return *this; }
INLINE VRefCount::VRefCount():_reference_count_(0) {}
INLINE VRefCount::VRefCount(const VRefCount&):_reference_count_(0) {}
INLINE VRefCount& VRefCount::operator=(const VRefCount&) { return *this; }

#undef INLINE

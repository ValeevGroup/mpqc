
#ifdef __GNUC__
#pragma implementation
#endif

#include <voidptr.h>

VoidPtr::VoidPtr(): ptr_(0) {}
VoidPtr::VoidPtr(void* d): ptr_(d) {}
VoidPtr::VoidPtr(const VoidPtr& d): ptr_(d.ptr_) {}
VoidPtr::~VoidPtr() {}

void* VoidPtr::getptr() { return ptr_; }
int VoidPtr::operator==(VoidPtr d) { return ptr_ == d.ptr_; }
int VoidPtr::operator<=(VoidPtr d) { return ptr_ <= d.ptr_; }
VoidPtr& VoidPtr::operator=(void* d) { ptr_ = d; return *this; }
VoidPtr& VoidPtr::operator=(VoidPtr d) { ptr_ = d.ptr_; return *this; }
VoidPtr::operator void*() { return ptr_; }

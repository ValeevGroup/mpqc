
#ifndef _libQC_voidptr_h
#define _libQC_voidptr_h

class VoidPtr
{
  private:
    void* ptr_;
  public:
    inline VoidPtr(): ptr_(0) {}
    inline VoidPtr(void* d): ptr_(d) {}
    inline VoidPtr(VoidPtr& d): ptr_(d.ptr_) {}
    inline ~VoidPtr() {}

    inline void* getptr() { return ptr_; }
    inline int operator==(VoidPtr d) { return ptr_ == d.ptr_; }
    inline int operator<=(VoidPtr d) { return ptr_ <= d.ptr_; }
    inline VoidPtr& operator=(void* d) { ptr_ = d; return *this; }
    inline VoidPtr& operator=(VoidPtr d) { ptr_ = d.ptr_; return *this; }
    inline operator void*() { return ptr_; }
  };

#endif

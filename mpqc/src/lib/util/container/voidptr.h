
#ifdef __GNUC__
#pragma interface
#endif

#ifndef _libQC_voidptr_h
#define _libQC_voidptr_h

class VoidPtr
{
  private:
    void* ptr_;
  public:
    VoidPtr();
    VoidPtr(void* d);
    VoidPtr(VoidPtr& d);
    ~VoidPtr();
    
    void* getptr();
    int operator==(VoidPtr d);
    int operator<=(VoidPtr d);
    VoidPtr& operator=(void* d);
    VoidPtr& operator=(VoidPtr d);
    operator void*();
  };

#endif

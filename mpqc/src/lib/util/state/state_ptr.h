
#ifndef _libqc_state_ptr_h
#define _libqc_state_ptr_h

class StateDataNum {
  private:
    int num_;
    void* ptr_;
  public:
    StateDataNum(int num):num_(num),ptr_(0) {}
    StateDataNum(int num,void*p):num_(num),ptr_(p) {}
    int num() { return num_; }
    void* ptr() { return ptr_; }
    void assign_ptr(void*p) { ptr_ = p; }
    int operator==(StateDataNum&num) { return num_ == num.num_; }
    int compare(StateDataNum&num);
  };

inline int StateDataNum::compare(StateDataNum&p)
{
  return (p.num_ == num_)?0:((p.num_<num_)?1:-1);
}

class StateDataPtr {
  private:
    void* ptr_;
  public:
    StateDataPtr(void*p):ptr_(p) {}
    void* ptr() { return ptr_; }
    int operator==(StateDataPtr&p) { return ptr_ == p.ptr_; }
    int compare(StateDataPtr&p);
  };

inline int StateDataPtr::compare(StateDataPtr&p)
{
  return (p.ptr_ == ptr_)?0:((p.ptr_<ptr_)?1:-1);
}

#endif


#ifndef _libqc_state_ptr_h
#define _libqc_state_ptr_h

#ifdef __GNUC__
#pragma interface
#endif

class StateDataNum {
  private:
    int num_;
    void* ptr_;
  public:
    StateDataNum(int num):num_(num),ptr_(0) {}
    StateDataNum(int num,void*p):num_(num),ptr_(p) {}
    int num() const { return num_; }
    void* ptr() const { return ptr_; }
    void assign_ptr(void*p) { ptr_ = p; }
    int operator==(const StateDataNum&num) const { return num_ == num.num_; }
    int compare(const StateDataNum&num) const;
  };

inline int StateDataNum::compare(const StateDataNum&p) const
{
  return (p.num_ == num_)?0:((p.num_<num_)?1:-1);
}

class StateDataPtr {
  private:
    int num_;
    void* ptr_;
  public:
    StateDataPtr(void*p):ptr_(p) {}
    void* ptr() const { return ptr_; }
    int num() const { return num_; }
    void assign_num(int n) { num_ = n; }
    int operator==(const StateDataPtr&p) const { return ptr_ == p.ptr_; }
    int compare(const StateDataPtr&p) const;
  };

inline int StateDataPtr::compare(const StateDataPtr&p) const
{
  return (p.ptr_ == ptr_)?0:((p.ptr_<ptr_)?1:-1);
}

#endif

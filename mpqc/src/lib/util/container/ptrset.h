
#ifndef _ptrset_h
#define _ptrset_h

#include <util/container/voidptrSet.h>
#include <util/container/voidptrAVLSet.h>

template <class Type>
class PtrSet
{
 private:

  VoidPtrSet* setimpl;

 public:
                       PtrSet():setimpl(new VoidPtrAVLSet) {}
  virtual              ~PtrSet() { delete setimpl; }

  int                   length() { return setimpl->length(); }
  int                   empty() { return setimpl->empty(); }

  virtual Pix           add(Type* item) { return setimpl->add((void*)item); }
  virtual void          del(Type* item) { setimpl->del((void*)item); }
  virtual int           contains(Type*  item)
                            { return setimpl->contains((void*)item); }

  virtual void          clear() { setimpl->clear(); }

  virtual Pix           first() { return setimpl->first(); }
  virtual void          next(Pix& i) { return setimpl->next(i); }
  virtual Type*         operator () (Pix i) { return (*setimpl)(i); }

  virtual int           owns(Pix i) { return setimpl->owns(i); }
  virtual Pix           seek(Type*  item) { return setimpl->seek((void*)item); }

  void                  operator |= (PtrSet<Type>& b)
                            { setimpl->operator|=(*b.setimpl); }
  void                  operator -= (PtrSet<Type>& b)
                            { setimpl->operator-=(*b.setimpl); }
  void                  operator &= (PtrSet<Type>& b)
                            { setimpl->operator&=(*b.setimpl); }

  int                   operator == (PtrSet<Type>& b)
                            { return setimpl->operator==(*b.setimpl); }
  int                   operator != (PtrSet<Type>& b)
                            { return setimpl->operator!=(*b.setimpl); }
  int                   operator <= (PtrSet<Type>& b)
                            { return setimpl->operator<=(*b.setimpl); }

};

#endif


#ifndef _pixmap_h
#define _pixmap_h

#include <stdio.h>
#include <util/container/pixvpMap.h>
#include <util/container/pixvpRAVLMap.h>

template <class Type>
class PixMap
{
 private:
  // Since PixVoidPtrMap does use const members, trickery must be used.
  inline PixVoidPtrMap& cimpl() const { return *((PixVoidPtrMap*)&impl); }

  PixVoidPtrRAVLMap impl;
  Pix guess;
 public:
  inline PixMap():guess((Pix)10000),impl((VoidPtr)0) {}
  virtual ~PixMap() {
      for (Pix i=first(); i!=0; next(i)) { delete (Type*) impl.contents(i); }
    }
  inline int length() const { return cimpl().length(); }
  inline int empty() const { return cimpl().empty(); }
  inline int contains(Pix key) const { return cimpl().contains(key); }
  inline void clear() { impl.clear(); }
  inline Type& operator[](Pix key) { 
      if (impl.seek(key) == 0) {
	  impl[key] = (void*) new Type;
	}
      Type* tmp = (Type*)(impl[key].getptr());
      return *tmp;
    }
  inline const Type& operator[](Pix key) const { 
      if (cimpl().seek(key) == 0) {
	  cimpl()[key] = (void*) new Type;
	}
      Type* tmp = (Type*)(cimpl()[key].getptr());
      return *tmp;
    }
  inline void del(Pix key) { impl.del(key); }
  inline Pix first() const { return cimpl().first(); }
  inline void next(Pix&i) const { cimpl().next(i); }
  inline Type& contents(Pix i) { return *((Type*)impl.contents(i).getptr()); }
  inline int owns(Pix i) const { return cimpl().owns(i); }
  inline Pix seek(Pix key) const { return cimpl().seek(key); }
  inline Pix key(Pix pix) const { return cimpl().key(pix); }

  // this finds a key that isn't already in use
  inline Pix newkey() {
    while (impl.seek(guess)) ((int)guess)++;
    return guess;
    }
};

#endif

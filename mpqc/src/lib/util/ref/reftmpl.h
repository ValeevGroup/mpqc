//
// reftmpl.h --- template reference counting class
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@limitpt.com>
// Maintainer: LPS
//
// This file is part of the SC Toolkit.
//
// The SC Toolkit is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The SC Toolkit is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the SC Toolkit; see the file COPYING.LIB.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//
// The U.S. Government is granted a limited license as per AL 91-7.
//

#ifdef __GNUG__
#pragma interface
#endif

//. The \clsnm{Ref} template class maintains references counts.
//.
//. Several of these operations can cause a reference to an object to be
//. replaced by a reference to a different object.  If a reference to a
//. nonnull object is eliminated, the object's reference count is
//. decremented and the object is deleted if the reference count becomes
//. zero.
//.
//. There also may be a to \srccd{operator T*()}, some compilers have
//. bugs that prevent it's use.  The \srccd{pointer()} member
//. works everywhere.
template <class T>
class  Ref  : private RefBase {
  private:
    T* p;
  public:
    //. Create a reference to a null object.
    Ref(): p(0) {}
    //. Create a reference to the object \vrbl{a}.
    Ref(T*a): p(a)
    {
      if (p) {
          if (DO_REF_CHECK_STACK(p)) {
              DO_REF_UNMANAGE(p);
              warn_ref_to_stack();
            }
          p->reference();
        }
    }
    //. Create a reference to the object referred to by \vrbl{a}.
    Ref(const Ref<T> &a): p(a.p)
    {
      if (p) p->reference();
    }
    //. Delete this reference to the object.  Decrement the object's reference
    //. count and delete the object if the count is zero.
    ~ Ref ()
    {
      clear();
    }
    //. Returns the reference counted object.  The behaviour is undefined if
    //. the object is null.
    T* operator->() const { return p; }
    //. Returns a pointer the reference counted object.
    T* pointer() const { return p; }

    REF_TYPE_CAST_DEC(T);
    //. Returns a C++ reference to the reference counted object.
    //. The behaviour is undefined if the object is null.
    T& operator *() const { return *p; };
    //. Return 1 if this is a reference to a null object.  Otherwise
    //. return 0.
    int null() const { return p == 0; }
    //. Return \srccd{!null()}.
    int nonnull() const { return p != 0; }
    //. Ordering relations are provided using the \clsnmref{Identity}
    //. class.
    int operator!=(const  T * a) const { return ne(p,a); }
    int operator==(const  T * a) const { return eq(p,a); }
    int operator==(const Ref<T> &a) const { return eq(p,a.p); }
    int operator>=(const  Ref<T> &a) const { return ge(p,a.p); }
    int operator<=(const  Ref<T> &a) const { return le(p,a.p); }
    int operator>(const  Ref<T> &a) const { return gt(p,a.p); }
    int operator<(const  Ref<T> &a) const { return lt(p,a.p); }
    int operator!=(const Ref<T> &a) const { return ne(p,a.p); }
    int compare(const Ref<T> &a) const {
      return eq(p,a.p)?0:((lt(p,a.p)?-1:1));
    }
    //. Refer to the null object.
    void clear()
    {
      if (p && p->dereference()<=0) {
          delete p;
        }
      p = 0;
    }
    //. Assignment operators.
    Ref<T>& operator=(const Ref<T> & c)
    {
      if (c.p) c.p->reference();
      clear();
      p=c.p;
      return *this;
    }
    Ref<T>& operator=(T* cr)
    {
      assign_pointer(cr);
      return *this;
    }
    void assign_pointer(T* cr)
    {
      if (cr) {
          if (DO_REF_CHECK_STACK(cr)) {
              DO_REF_UNMANAGE(cr);
              warn_ref_to_stack();
            }
          cr->reference();
        }
      clear();
      p = cr;
    }
    //. Miscellaneous members.
    void check_pointer() const
    {
      if (p && p->nreference() <= 0) {
          warn_bad_ref_count();
        }
    }
    void ref_info(ostream& os) const
    {
      RefBase::ref_info(p,os);
    }
    void warn(const char*s) const { RefBase::warn(s); }
};

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:

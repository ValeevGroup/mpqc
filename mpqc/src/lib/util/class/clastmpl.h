//
// clastmpl.h
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

//. \clsnm{DCRef} is similar to the \clsnmref{Ref} class
//. template, however it provides new implementations of
//. constructor and assignment operators that take advantage
//. of the castdown ability of \clsnmref{DescribedClass}.  The
//. standard reference classes allows only assignment to
//. references or pointers to objects of the same exact type
//. as the contained object.  \clsnm{DCRef} has constructor
//. and assignment operators that take generic
//. \clsnmref{DescribedClass} references and pointers as
//. arguments.  They use the contained type's static
//. \srccd{castdown} operator to convert the target
//. \clsnmref{DescribedClass} object into an object of the
//. appropriate exact type.  If the cast fails a reference to
//. null is assigned.
template <class T>
class DCRef : public DCRefBase {
  protected:
    T* p;
  public:
    //. Implementation of \srccd{\clsnmref{DCRefBase}::parentpointer()}.
    DescribedClass* parentpointer() const { return p; }
    //. Create a reference to a null object.
    DCRef(): p(0) {}
    //. Create a reference to the object \vrbl{a}.
    DCRef(T*a): p(a)
    {
      reference(p);
    }
    //. Create a reference to the object \vrbl{a}.
    DCRef(const DCRef<T> &a): p(a.p)
    {
      reference(p);
    }
    //. Create a reference to the object \vrbl{a}.  Do a safe
    //. \srccd{castdown} to convert \vrbl{a} to the appropiate type.
    DCRef(const DCRefBase&a) {
        p = T::castdown(a.parentpointer());
        reference(p);
      }
    //. Delete this reference to the object.  Decrement the object's reference
    //. count and delete the object if the count is zero.
    ~ DCRef ()
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
    int compare(const DCRef<T> &a) const {
        return (eq(p,a.p)?0:(lt(p,a.p)?-1:1));
      }
    int operator==(const DCRef<T> &a) const { return eq(p,a.p); }
    int operator!=(const DCRef<T> &a) const { return ne(p,a.p); }
    int operator>=(const DCRef<T> &a) const { return ge(p,a.p); }
    int operator<=(const DCRef<T> &a) const { return le(p,a.p); }
    int operator> (const DCRef<T> &a) const { return gt(p,a.p); }
    int operator< (const DCRef<T> &a) const { return lt(p,a.p); }
    //. Delete this reference to the object.  Decrement the object's reference
    //. count and delete the object if the count is zero.
    void clear()
    {
      dereference(p);
      p = 0;
    }
    //. Assigment operators.
    DCRef<T>& operator=(const DCRef<T> & c)
    {
      reference(c.p);
      clear();
      p=c.p;
      return *this;
    }
    DCRef<T>& operator=(T* cr)
    {
      assign_pointer(cr);
      return *this;
    }
    DCRef<T>& operator=(const DCRefBase&a) {
        T* cr = T::castdown(a.parentpointer());
        reference(cr);
        clear();
        p = cr;
        return *this;
      }
    void assign_pointer(T* cr)
    {
      reference(cr);
      clear();
      p = cr;
    }
};

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:

//
// ref.h --- definitions of the reference counting classes
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

//   This is the main include file for the reference counting classes.
// This includes two other files: reftmpl.h and refmacr.h.  The
// former is a template declaration for the reference counted classes
// and the latter is generated from the former by a perl script and
// provides CPP macros that declare reference counting classes.
//
//   The behaviour of the package can be modified with the following five
// macros, each of which should be undefined, 0, or 1:
//
// REF_CHECK_STACK:  If this is 1 referenced objects are checked to see if they
// reside on the stack, in which case storage for the object is not managed,
// if management is enabled.  This feature can be confused by multiple threads
// and memory checking libraries.
//
// REF_MANAGE:  If this is 1 the manage and unmanage members are enabled.
//
// REF_CHECK_MAX_NREF:  If this is 1 the reference count is checked before
// it is incremented to make sure it isn't too big.
//
// REF_CHECK_MIN_NREF:  If this is 1 the reference count is checked before
// it is decremented to make sure it isn't already zero.
//
// REF_USE_LOCKS:  If this is 1 then critical regions are locked before they
// are entered.  This prevents erroneous behavior when multiple threads
// share reference counted objects.  This will slow down certain operations,
// so it should be set to 0 if your application does not need to be thread
// safe.
//
// If a macro is undefined, then the behaviour is architecture
// dependent--usually, the macro will be set to 1 in this case.
// For maximum efficiency and for normal operation after the program is
// debugged, compile with all of the above macros defined to zero.
// This can also be done with -DREF_OPTIMIZE.
//
//   An include file can be used to set these options as well.  This has
// the advantage that dependency checking will force an automatic
// recompile of all affected files if the options change.  The file
// <mpqc_config.h> will be include if -DHAVE_CONFIG_H is specified.
//
//   Note that all source code that uses references must be compiled with
// the same value REF_MANAGE.  Changing this can change the storage layout
// and the interpretation of the reference count data.


#ifndef _util_ref_ref_h
#define _util_ref_ref_h

#include <iostream>
#include <stdlib.h>
#include <limits.h>

#include <util/ref/identity.h>

#ifdef HAVE_CONFIG_H
#include <mpqc_config.h>
#endif

#ifdef REF_OPTIMIZE
#ifndef REF_CHECK_STACK
# define REF_CHECK_STACK   0
#endif
#ifndef REF_MANAGE
# define REF_MANAGE        0
#endif
#ifndef REF_CHECK_MAX_NREF
# define REF_CHECK_MAX_NREF 0
#endif
#ifndef REF_CHECK_MIN_NREF
# define REF_CHECK_MIN_NREF 0
#endif
#endif

#ifdef SUNMOS
#ifndef REF_CHECK_STACK
#define REF_CHECK_STACK 0
#endif
#else
#ifndef REF_CHECK_STACK
#define REF_CHECK_STACK 0
#endif
#endif

#ifndef REF_MANAGE
#define REF_MANAGE 1
#endif

#ifndef REF_CHECK_MAX_NREF
#define REF_CHECK_MAX_NREF 1
#endif

#ifndef REF_CHECK_MIN_NREF
#define REF_CHECK_MIN_NREF 1
#endif

#ifndef REF_USE_LOCKS
#if HAVE_STHREAD || HAVE_CREATETHREAD || HAVE_PTHREAD
#define REF_USE_LOCKS 1
#endif
#else // REF_USE_LOCKS
#warning "REF_USE_LOCKS not defined"
#endif

#ifndef REF_ALWAYS_USE_LOCKS
#define REF_ALWAYS_USE_LOCKS 1
#endif

#if REF_CHECK_STACK
#include <unistd.h>
#ifndef HAVE_SBRK_DEC
extern "C" void * sbrk(ssize_t);
#endif
#define DO_REF_CHECK_STACK(p) (((void*) (p) > sbrk(0)) && (p)->managed())
#else // REF_CHECK_STACK
#define DO_REF_CHECK_STACK(p) (0)
#endif // REF_CHECK_STACK

#if REF_MANAGE
#define DO_REF_UNMANAGE(p) ((p)->unmanage())
#else // REF_MANAGE
#define DO_REF_UNMANAGE(p)
#endif // REF_MANAGE

#if REF_USE_LOCKS
#define __REF_LOCK__(p) p->lock_ptr()
#define __REF_UNLOCK__(p) p->unlock_ptr()
#if REF_ALWAYS_USE_LOCKS
#define __REF_INITLOCK__() use_locks(true)
#else
#define __REF_INITLOCK__() ref_lock_ = 0xff
#endif
#else
#define __REF_LOCK__(p)
#define __REF_UNLOCK__(p)
#define __REF_INITLOCK__()
#endif

namespace sc {

typedef unsigned long refcount_t;

/** The base class for all reference counted objects.  If multiple
    inheritance is used, RefCount must be virtually inherited from,
    otherwise references to invalid memory will likely result.

    Reference counting information is usually maintained by smart
    pointer classes Ref, however this mechanism can be
    supplemented or replaced by directly using the public
    interface to RefCount.

    The unmanage() member is only needed for special cases where memory
    management must be turned off.  For example, if a reference counted
    object is created on the stack, memory management mechanisms based on
    reference counting must be prohibited from deleting it.  The unmanage()
    member accomplishes this, but a better solution would be to allocate
    the object on the heap with new and let a smart pointer manage the
    memory for the object.

    When using a debugger to look at reference counted objects the count is
    maintained in the _reference_count_ member.  However, this member is
    encoded so that memory overwrites can be sometimes detected.  Thus,
    interpretation of _reference_count_ is not always straightforward.

*/

class RefCount: public Identity {
  private:
#if REF_MANAGE
#  define REF_MAX_NREF (UINT_MAX - 1)
#  define REF_MANAGED_CODE UINT_MAX
#else
#  define REF_MAX_NREF UINT_MAX
#endif
    unsigned int _reference_count_;
#if REF_USE_LOCKS
    unsigned char ref_lock_;
#endif

    void error(const char*) const;
    void too_many_refs() const;
    void not_enough_refs() const;
  protected:
    RefCount(): _reference_count_(0) {
        __REF_INITLOCK__();
        //std::cout << "ref_lock_ = " << (int) ref_lock_ << std::endl;
      }
    RefCount(const RefCount&): _reference_count_(0) {
        __REF_INITLOCK__();
        //std::cout << "ref_lock_ = " << (int) ref_lock_ << std::endl;
      }

    // Assigment should not overwrite the reference count.
    RefCount& operator=(const RefCount&) { return *this; }
  public:
    virtual ~RefCount();

    /// Lock this object.
    int lock_ptr() const;
    /// Unlock this object.
    int unlock_ptr() const;

    /// start and stop using locks on this object
    void use_locks(bool inVal);

    /// Return the reference count.
    refcount_t nreference() const {
#       if REF_MANAGE
        if (!managed()) return 1;
#       endif
        return _reference_count_;
      }

    /// Increment the reference count and return the new count.
    refcount_t reference() {
#       if REF_MANAGE
        if (!managed()) return 1;
#       endif
        __REF_LOCK__(this);
#       if REF_CHECK_MAX_NREF
        if (_reference_count_ >= REF_MAX_NREF) too_many_refs();
#       endif
        _reference_count_++;
        refcount_t r = _reference_count_;
        __REF_UNLOCK__(this);
        return r;
      }

    /// Decrement the reference count and return the new count.
    refcount_t dereference() {
#       if REF_MANAGE
        if (!managed()) return 1;
#       endif
        __REF_LOCK__(this);
#       if REF_CHECK_MIN_NREF
        if (_reference_count_ == 0) not_enough_refs();
#       endif
        _reference_count_--;
        refcount_t r = _reference_count_;
        __REF_UNLOCK__(this);
        return r;
      }

#if REF_MANAGE
    int managed() const {
        return _reference_count_ != REF_MANAGED_CODE;
      }
    /** Turn off the reference counting mechanism for this object.
        The value returned by nreference() will always be
        1 after this is called.  The ability to unmanage()
        objects must be turned on at compile time by defining
        REF_MANAGE.  There is a slight performance penalty. */
    void unmanage() {
        _reference_count_ = REF_MANAGED_CODE;
      }
#else // REF_MANAGE
    /// Return 1 if the object is managed.  Otherwise return 0.
    int managed() const { return 1; }
#endif // REF_MANAGE
};

/** Provides a few utility routines common to all
    Ref template instantiations.
*/
class RefBase {
  protected:
    /// Print a warning message.
    void warn ( const char * msg) const;
    /// Called when stack data is referenced.
    void warn_ref_to_stack() const;
    /// Called when the deletion of stack data is skipped.
    void warn_skip_stack_delete() const;
    /// Called when the reference count is corrupted.
    void warn_bad_ref_count() const;
    /// Print information about the reference.
    void ref_info(RefCount*p,std::ostream& os) const;
    void ref_info(std::ostream& os) const;
    void check_pointer() const;
    void reference(RefCount *);
    int dereference(RefCount *);
  public:
    RefBase() {};
    virtual ~RefBase();
    /// Returns the DescribedClass pointer for the contained object.
    virtual RefCount* parentpointer() const = 0;
    /** Requires that a nonnull reference is held.  If not,
        the program will abort. */
    void require_nonnull() const;
};

/** A template class that maintains references counts.

    Several of these operations can cause a reference to an object to be
    replaced by a reference to a different object.  If a reference to a
    nonnull object is eliminated, the object's reference count is
    decremented and the object is deleted if the reference count becomes
    zero.

    There also may be a to convert to T*, where T is the type of the object
    which Ref references.  Some compilers have bugs that prevent the use of
    operator T*().  The pointer() member should be used instead.

*/
template <class T>
class  Ref  : public RefBase {
  public:
    typedef T element_type;
  private:
    T* p;
  public:
    /// Create a reference to a null object.
    Ref(): p(0) {}
    /// Create a reference to the object a.
    Ref(T*a) : p(0)
    {
      if (a) {
          p = a;
          reference(p);
        }
    }
    /// Create a reference to the object referred to by a.
    Ref(const Ref<T> &a) : p(0)
    {
      if (a.pointer()) {
          p = a.pointer();
          reference(p);
        }
    }
    /// Create a reference to the object referred to by a.
    template <class A> Ref(const Ref<A> &a): p(0)
    {
      if (a.pointer()) {
          p = a.pointer();
          reference(p);
        }
    }
//      /** Create a reference to the object a.  Do a
//          dynamic_cast to convert a to the appropiate type. */
//      Ref(const RefBase&a) {
//          p = dynamic_cast<T*>(a.parentpointer());
//          reference(p);
//        }
//      /** Create a reference to the object a.  Do a
//          dynamic_cast to convert a to the appropiate type. */
//      Ref(RefCount*a): p(0) {
//        operator<<(a);
//        }
    /** Delete this reference to the object.  Decrement the object's reference
        count and delete the object if the count is zero. */
    ~Ref()
    {
      clear();
    }
    /** Returns the reference counted object.  The behaviour is undefined if
        the object is null. */
    T* operator->() const { return p; }
    /// Returns a pointer the reference counted object.
    T* pointer() const { return p; }
    /// Implements the parentpointer pure virtual in the base class.
    RefCount *parentpointer() const { return p; }

    operator T*() const { return p; }
    /** Returns a C++ reference to the reference counted object.
        The behaviour is undefined if the object is null. */
    T& operator *() const { return *p; };
    /** Return 1 if this is a reference to a null object.  Otherwise
        return 0. */
    int null() const { return p == 0; }
    /// Return !null().
    int nonnull() const { return p != 0; }
    /** A variety of ordering and equivalence operators are provided using
        the Identity class. */
    template <class A> int operator==(const Ref<A>&a) const
        { return eq(p,a.pointer()); }
    template <class A> int operator>=(const Ref<A>&a) const
        { return ge(p,a.pointer()); }
    template <class A> int operator<=(const Ref<A>&a) const
        { return le(p,a.pointer()); }
    template <class A> int operator>(const Ref<A>&a) const
        { return gt(p,a.pointer()); }
    template <class A> int operator<(const Ref<A>&a) const
        { return lt(p,a.pointer()); }
    template <class A> int operator!=(const Ref<A>&a) const
        { return ne(p,a.pointer()); }
    /** Compare two objects returning -1, 0, or 1. Similar
        to the C library routine strcmp. */
    int compare(const Ref<T> &a) const {
      return eq(p,a.p)?0:((lt(p,a.p)?-1:1));
    }
    /// Refer to the null object.
    void clear()
    {
      if (p) {
          int ref = dereference(p);
          if (ref == 0)
              delete p;
          p = 0;
        }
    }
    /// Assignment to c.
    Ref<T>& operator=(const Ref<T> & c)
    {
      T *cp = c.pointer();
      if (cp) {
          cp->reference();
          clear();
          p=cp;
        }
      else {
          clear();
        }
      return *this;
    }
    /// Assignment to c.
    template <class A> Ref<T>& operator=(const Ref<A> & c)
    {
      A *cp = c.pointer();
      if (cp) {
          cp->reference();
          clear();
          p=cp;
        }
      else {
          clear();
        }
      return *this;
    }
    /// Assignment to the object that a references using dynamic_cast.
    Ref<T>& operator<<(const RefBase&a) {
        T* cr = dynamic_cast<T*>(a.parentpointer());
        if (cr) {
            reference(cr);
            clear();
          }
        p = cr;
        return *this;
      }
    /** Assigns to the given base class pointer using dynamic_cast.  If
        the dynamic_cast fails and the argument is nonnull and has a
        reference count of zero, then it is deleted. */
    Ref<T>& operator<<(RefCount *a) {
        T* cr = dynamic_cast<T*>(a);
        if (cr) assign_pointer(cr);
        else if (a && a->nreference() <= 0) delete a;
        return *this;
      }
    /// Assignment to cr.
    Ref<T>& operator=(T* cr)
    {
      assign_pointer(cr);
      return *this;
    }
    /// Assignment to cr.
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
    /// Check the validity of the pointer.
    void check_pointer() const
    {
      if (p && p->nreference() <= 0) {
          warn_bad_ref_count();
        }
    }
    /// Print information about the reference to os.
    void ref_info(std::ostream& os) const
    {
      RefBase::ref_info(p,os);
    }
    /// Print a warning concerning the reference.
    void warn(const char*s) const { RefBase::warn(s); }
};

/** this functor can be used as a binary predicate for standard algorithms. For example,
    it can be used as an argument to std::hash_map that uses keys of sc::Ref<T> type.
    Optional EqualTo argument is a default-constructible binary predicate
    that takes 2 values of type T. By default it calls operator==(const T&, const T&).
  */
template <typename T, typename EqualTo = std::equal_to<T> >
struct RefObjectEqual {
  bool operator()(const Ref<T>& obj1, const Ref<T>& obj2) {
    EqualTo equal_functor;
    return equal_functor(*obj1,*obj2);
  }
};

}

#endif

// ///////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:

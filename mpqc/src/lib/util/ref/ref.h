
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
// if management is enabled.
//
// REF_CHECKSUM:  If this is 1 checksums of the reference count are kept
// and checked to see if an object has been overwritten.
//
// REF_MANAGE:  If this is 1 the manage and unmanage members are enabled.
//
// REF_CHECK_MAX_NREF:  If this is 1 the reference count is checked before
// it is incremented to make sure it isn't too big.
//
// REF_CHECK_MIN_NREF:  If this is 1 the reference count is checked before
// it is decremented to make sure it isn't already zero.
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
// <config.h> will be include if -DHAVE_CONFIG_H is specified.
//
//   Note that all source code that uses references must be compiled with
// the same value for REF_CHECKSUM and REF_MANAGE.  Changing these can
// change the storage layout and the interpretation of the reference count
// data.


#ifdef __GNUC__
#pragma interface
#endif

#ifndef _util_ref_ref_h
#define _util_ref_ref_h

#include <stdio.h>
#include <stdlib.h>

#include <util/ref/identity.h>

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef REF_OPTIMIZE
#ifndef REF_CHECK_STACK
# define REF_CHECK_STACK   0
#endif
#ifndef REF_CHECKSUM
# define REF_CHECKSUM      0
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
#define REF_CHECK_STACK 1
#endif
#endif

#ifndef REF_CHECKSUM
#define REF_CHECKSUM 1
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

#if REF_CHECK_STACK
#include <unistd.h>
#ifndef HAVE_SBRK_DEC
extern "C" void * sbrk(int);
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

typedef unsigned long refcount_t;

//. \srccd{VRefCount} is the base class for all reference counted
//. objects.  If multiple inheritance is used, \srccd{VRefCount}
//. must be virtually inherited from, otherwise references to
//. invalid memory will likely result.
//.
//. Reference counting information is usually maintained by smart
//. pointer classes (\clsnm{Ref}), however this mechanism can be
//. supplemented or replaced by directly using the public
//. interface to \srccd{VRefCount}.
//.
//. The \srccd{unmanage()} member is only needed for special cases where
//. memory management must be turned off.  For example, if a reference
//. counted object is created on the stack, memory management
//. mechanisms based on reference counting must be prohibited from deleting it.
//. The \srccd{unmanage()} member accomplishes this, but a better solution
//. would be to allocate the object on the heap with \srccd{new} and let
//. a class declared using the \srccd{REF\_dec} macro manage the memory for
//. the object.
//. 
//. When using a debugger to look at reference counted objects the
//. count is maintained in the \srccd{\_reference\_count\_} member.
//. However, this member is encoded so that memory overwrites can be
//. sometimes detected.  Thus, interpretation of
//. \srccd{\_reference\_count\_} is not straightforward and will change
//. as the implementation of \srccd{VRefCount} becomes more
//. sophisticated.
class VRefCount: public Identity {
  private:
#if REF_CHECKSUM
#   if REF_MANAGE
#     define REF_MAX_NREF   0xfffffe
#     define REF_MANAGED_CODE 0xffffff
#   else
#     define REF_MAX_NREF 0xffffff
#   endif
    unsigned int _reference_count_:24;
    unsigned int _checksum_:8;
#else
#   if REF_MANAGE
#     define REF_MAX_NREF (UINT_MAX - 1)
#     define REF_MANAGED_CODE UINT_MAX
#   else
#     define REF_MAX_NREF UINT_MAX
#   endif
    unsigned int _reference_count_;
#endif

    void error(const char*) const;
    void too_many_refs() const;
    void not_enough_refs() const;
#if REF_CHECKSUM
    void bad_checksum() const;
    void check_checksum() const {
        if (_checksum_ != (~(( ( _reference_count_      & 0xff)
                               ^((_reference_count_)>>8  & 0xff)
                               ^((_reference_count_)>>16 & 0xff)))&0xff)) {
            bad_checksum();
          }
      }
    void update_checksum() {
        _checksum_ = ~( ( _reference_count_      & 0xff)
                        ^((_reference_count_)>>8  & 0xff)
                        ^((_reference_count_)>>16 & 0xff));
      }
#endif // REF_CHECKSUM
  protected:
    VRefCount(): _reference_count_(0) {
#       if REF_CHECKSUM
        update_checksum();
#       endif
      }
    VRefCount(const VRefCount&): _reference_count_(0) {
#       if REF_CHECKSUM
        update_checksum();
#       endif
      }

    // Assigment should not overwrite the reference count.
    VRefCount& operator=(const VRefCount&) { return *this; }
  public:
    virtual ~VRefCount();

    //. Return the reference count.
    refcount_t nreference() const {
#       if REF_MANAGE
        if (!managed()) return 1;
#       endif
#       if REF_CHECKSUM
        check_checksum();
#       endif
        return _reference_count_;
      }

    //. Increment the reference count and return the new count.
    refcount_t reference() {
#       if REF_MANAGE
        if (!managed()) return 1;
#       endif
#       if REF_CHECKSUM
        check_checksum();
#       endif
#       if REF_CHECK_MAX_NREF
        if (_reference_count_ >= REF_MAX_NREF) too_many_refs();
#       endif
        _reference_count_++;
#       if REF_CHECKSUM
        update_checksum();
#       endif
        return _reference_count_;
      }

    //. Decrement the reference count and return the new count.
    refcount_t dereference() {
#       if REF_MANAGE
        if (!managed()) return 1;
#       endif
#       if REF_CHECKSUM
        check_checksum();
#       endif
#       if REF_CHECK_MIN_NREF
        if (_reference_count_ == 0) not_enough_refs();
#       endif
        _reference_count_--;
#       if REF_CHECKSUM
        update_checksum();
#       endif
        return _reference_count_;
      }

#if REF_MANAGE
    int managed() const {
#       if REF_CHECKSUM
        check_checksum();
#       endif // REF_CHECKSUM
        return _reference_count_ != REF_MANAGED_CODE;
      }
    //. Turn off the reference counting mechanism for this object.
    //. The value returned by \srccd{nreference()} will always be
    //. 1 after this is called.  The ability to \srccd{unmanage()}
    //. objects must be turned on at compile time by defining
    //. \srccd{REF\_MANAGE}.  There is a slight performance penalty.
    void unmanage() {
        _reference_count_ = REF_MANAGED_CODE;
#       if REF_CHECKSUM
        update_checksum();
#       endif // REF_CHECKSUM
      }
#else // REF_MANAGE
    //. Return 1 if the object is managed.  Otherwise return 0.
    int managed() const { return 1; }
#endif // REF_MANAGE
};

//. \clsnm{RefBase} provide a few utility routines common to all
//. \clsnmref{Ref} template instantiations.
class RefBase {
  public:
    //. Reference utility routines.
    void warn ( const char * msg) const;
    void warn_ref_to_stack() const;
    void warn_skip_stack_delete() const;
    void warn_bad_ref_count() const;
    void ref_info(VRefCount*p,FILE*fp) const;
};

#ifdef TYPE_CONV_BUG
#  define REF_TYPE_CAST_DEC(T)
#else
#  define REF_TYPE_CAST_DEC(T) operator T*() const { return p; }
#endif

// The template reference declaration.
#include <util/ref/reftmpl.h>

// The macro reference declaration.
#include <util/ref/refmacr.h>
#define REF_dec(T) Ref_declare(T)
#define REF_def(T)

#endif

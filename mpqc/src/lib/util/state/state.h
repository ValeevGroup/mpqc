//
// state.h
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

#ifndef _util_state_state_h
#define _util_state_state_h

#ifdef __GNUC__
#pragma interface
#endif

#include <util/class/class.h>

/** This macro declares a smart pointer type with savable state
    capabilities.  If the class name is T, the smart pointer type will be
    RefT.  */

#define SavableState_REF_dec(T) SavableState_named_REF_dec(Ref ## T,T)

/** This macro gives the class definition a smart pointer type with savable
    state capabilities.  If the class name is T, the smart pointer type
    will be RefT.  */

#define SavableState_REF_def(T) SavableState_named_REF_def(Ref ## T,T)

#ifdef USE_REF_MACROS
#define SavableState_named_REF_dec(refname,T)				      \
   DCRef_declare(T); SSRef_declare(T); typedef class SSRef ## T refname;
#else
#define SSRef_declare(T) typedef class SSRef<T> SSRef ## T;
#define SavableState_named_REF_dec(refname,T) typedef class SSRef<T> refname;
#endif
#define SavableState_named_REF_def(refname,T)

// This does forward declarations of REF classes.
#ifdef USE_REF_MACROS
#define SavableState_REF_fwddec(T) class SSRef ## T; \
                                   typedef class SSRef ## T Ref ## T;
#else
#define SavableState_REF_fwddec(T) class T; \
                                   typedef class SSRef<T> Ref ## T;
#endif

class StateIn;
class StateOut;
class TranslateDataIn;
class TranslateDataOut;

// If multiple inheritance is used, SavableState should be a virtual
// parent.  This breaks some compilers so the virtual_base macro should
// be used instead of the virtual keyword (see identity.h).  Since for StateIn
// CTORs the SavableState base class must be initialized with the
// StateIn object, the following macro should be used to determine if
// nondirect descendants of SavableState need to call the SavableState
// CTOR.  s is the StateIn
#ifdef NO_VIRTUAL_BASES
#  define maybe_SavableState(s)
#  define maybe_SavableState_alone(s)
#else
#  define maybe_SavableState(s) SavableState(s),
#  define maybe_SavableState_alone(s) :SavableState(s)
#endif

/** Base class for objects that can save/restore state.
 */
class SavableState: public DescribedClass {
#   define CLASSNAME SavableState
#   include <util/class/classda.h>
  protected:
    SavableState();
    SavableState(const SavableState&);
    // helper for save_object_state overrides
    void save_object_state_(StateOut&, const ClassDesc *);
#ifndef __GNUC__
  public:
#endif
    SavableState& operator=(const SavableState&);
  public:
    virtual ~SavableState();

    // save functions

    /** Save the state of the object as specified by the
        StateOut object.  This routine saves the state
        of the object (which includes the nonvirtual bases),
        the virtual bases, and type information.  The default
        implementation should be adequate. */
    void save_state(StateOut&);

    // Like save_state(StateOut&), but will handle null pointers correctly.
    static void save_state(SavableState*s, StateOut&);

    /** This can be used for saving state when the exact type of
        the object is known for both the save and the restore.  To
        restore objects saved in this way the user must directly
        invoke the object's StateIn& constructor. */
    void save_object_state(StateOut&);

    /** Save the virtual bases for the object.
        This must be done in the same order that the ctor
        initializes the virtual bases.  This does not include
        the DescribedClass and SavableState virtual base classes.
        This must be implemented by the user if the class has other
        virtual bases.  (These virtual bases must come after
        SavableState, if SavableState is virtual.) */
    virtual void save_vbase_state(StateOut&);

    /** Save the base classes (with save_data_state) and the members in the
        same order that the StateIn CTOR initializes them.  This must be
        implemented by the derived class if the class has data. */
    virtual void save_data_state(StateOut&);

    // restore functions

    /** Restores objects saved with save_state.  The
        exact type of the next object in si can be any
        type publically derived from the SavableState.
        Derived classes implement a similar static function that
        returns a pointer to the derived class.  If the objectname is
        given the directory will be consulted to find and restore
        that object. */
    static SavableState* restore_state(StateIn& si);
    /** Like restore_state, but keyword is used to override
        values while restoring. */
    static SavableState* key_restore_state(StateIn& si,
                                           const char *keyword);
    static SavableState* dir_restore_state(StateIn& si,
                                           const char *objectname,
                                           const char *keyword = 0);

  protected:

    /** Each derived class StateIn CTOR handles the restore corresponding
        to calling save_object_state, save_vbase_state, and save_data_state
        listed above.  All derived class StateIn& constructors must invoke
        the SavableState(StateIn&) constructor. */
    SavableState(StateIn&);
  };

// just do a fwddec here, since StateIn and StateOut are not
// yet declared.
SavableState_REF_fwddec(SavableState);

// //////////////////////////////////////////////////////////////////

/** Provides a few utility routines common to all SSRef instantiations.
 */
class SSRefBase {
  protected:
    void check_castdown_result(void*, SavableState *, const ClassDesc *);
  public:
    virtual SavableState *sspointer() = 0;
    virtual void dir_restore_state(StateIn& si, const char *objectname,
                                   const char *keyword = 0) = 0;
    void key_restore_state(StateIn& si, const char *keyword) {
      dir_restore_state(si,0,keyword);
    }
    void restore_state(StateIn& si) {
      dir_restore_state(si,0,0);
    }
    void save_data_state(StateOut&);
    /// Save the state of the reference.
    void save_state(StateOut&);
};

// Include the smart pointer to SavableState templates and macros.
#include <util/state/stattmpl.h>
#ifdef USE_REF_MACROS
#include <util/state/statmacr.h>
#endif

SavableState_REF_dec(SavableState);

// //////////////////////////////////////////////////////////////////

#ifndef __GNUC__
static SavableState * att_hack_job(StateIn&si)
{
  return SavableState::restore_state(si);
}
#endif

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:

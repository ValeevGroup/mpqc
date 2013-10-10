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

#include <util/class/class.h>

namespace sc {

class StateIn;
class StateOut;
class TranslateDataIn;
class TranslateDataOut;

/// @addtogroup CoreState
/// @{

/** Base class for objects that can save/restore state.
 */
class SavableState: virtual public DescribedClass {
  protected:
    SavableState();
    SavableState(const SavableState&);
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

  /// useful as a dummy template argument
  struct DummySavableState : virtual public SavableState {
  public:
    typedef DummySavableState this_type;

    DummySavableState();
    DummySavableState(StateIn& si);
    void save_data_state(StateOut&);

  private:
    static ClassDesc class_desc_;
  };

/// @}
// end of addgroup State

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:

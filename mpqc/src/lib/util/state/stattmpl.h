//
// stattmpl.h
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

/** A smart reference to a SavableState object.

    SSRef is similar to the DCRef class
    template, however it provides members to save and restore
    the state of contained objects.  Null references are
    saved and restored correctly as well.
*/
template <class T>
class SSRef: public DCRef<T>, public SSRefBase {
  public:
    /// Make a reference to null.
    SSRef() {}
    /// Make a reference to the object held by o.
    SSRef (const SSRef<T> & o): DCRef<T> (o) {}
    /// Make a reference to o.
    SSRef (T * o): DCRef<T> (o) {}
    /// Make a reference to the object held by o.
    SSRef (const DCRefBase&o): DCRef<T> (o) {}
    /// Release the reference, possible freeing the object.
    ~SSRef () {}
    /// Refer to cr.
    SSRef<T>& operator=(T* cr) {
        DCRef<T>::operator=(cr); return *this; }
    /// Make a reference to the object held by c.
    SSRef<T>& operator=(const DCRefBase & c) {
        DCRef<T>::operator=(c); return *this; }
    /// Make a reference to the object held by c.
    SSRef<T>& operator=(const SSRef<T>&c) {
        DCRef<T>::operator=(c); return *this; }
    /// Restore the object held in s.
    SSRef (StateIn&s) { restore_state(s); }
    /// Return a SavableState pointer to the contained object.
    SavableState *sspointer() { return p; }
    /** Restore the state of the reference from a StateIn
        with an object directory. */
    void dir_restore_state(StateIn& si, const char *objectname,
                           const char *keyword = 0) {
      SavableState* ss
          = SavableState::dir_restore_state(si, objectname, keyword);
      T* t = T::castdown(ss);
      check_castdown_result((void*)t,ss,T::static_class_desc());
      assign_pointer(t);
    };
};

// ///////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:

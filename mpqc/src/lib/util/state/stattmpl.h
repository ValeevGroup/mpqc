//
// stattmpl.h
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@ca.sandia.gov>
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

//. \clsnm{SSRef} is similar to the \clsnmref{DCRef} class
//. template, however it provides members to save and restore
//. the state of contained objects.  Null references are
//. saved and restored correctly as well.
template <class T>
class SSRef: public DCRef<T>, public SSRefBase {
  public:
    SSRef() {}
    SSRef (const SSRef<T> & o): DCRef<T> (o) {}
    SSRef (T * o): DCRef<T> (o) {}
    SSRef (const DCRefBase&o): DCRef<T> (o) {}
    ~SSRef () {}
    SSRef<T>& operator=(T* cr) {
        DCRef<T>::operator=(cr); return *this; }
    SSRef<T>& operator=(const DCRefBase & c) {
        DCRef<T>::operator=(c); return *this; }
    SSRef<T>& operator=(const SSRef<T>&c) {
        DCRef<T>::operator=(c); return *this; }
    SSRef (StateIn&s) { restore_state(s); }
    SavableState *sspointer() { return p; }
    //. Restore the state of the reference.
    void restore_state(StateIn&si) {
        SavableState* ss = restore_ss(si);
        T* t = T::castdown(ss);
        check_castdown_result((void*)t,ss);
        assign_pointer(t);
      };
};

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// eval: (c-set-style "CLJ")
// End:

//
// identity.h --- definition of the Identity class
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

#ifndef _util_ref_identity_h
#define _util_ref_identity_h

#ifdef __GNUG__
#pragma interface
#endif


#ifdef NO_VIRTUAL_BASES
#  define virtual_base
#else
#  define virtual_base virtual
#endif

class Identity;

//. Identifier's are used to distinguish and order objects.  On
//. many architectures a pointer to the object will suffice, but
//. the C++ standard only guarantees that this works for two pointers
//. pointing within the same structure or array.  Classes need to
//. inherit from Identity to use this mechanism.  Identity,
//. Identifier, and the shorthand boolean operations may have to
//. be modified for some different architectures.
class Identifier {
  private:
    const void* id;
  public:
    //. Create an \clsnm{Identifier} for a null object.
    Identifier(): id(0) {}
    //. Create an \clsnm{Identifier} for the given object.
    Identifier(const Identity* i): id((void*)i) {}
    Identifier(const Identifier& i): id(i.id) {}
    //. The destructor does nothing.
    ~Identifier() {}

    //. Assign to the given \clsnm{Identifier}.
    void operator = (const Identifier& i) { id = i.id; }

    //. Ordering relationships for objects.
    int operator < (const Identifier&i) const { return id < i.id; }
    int operator > (const Identifier&i) const { return id > i.id; }
    int operator == (const Identifier&i) const { return id == i.id; }
    int operator <= (const Identifier&i) const { return id <= i.id; }
    int operator >= (const Identifier&i) const { return id >= i.id; }
    int operator != (const Identifier&i) const { return id != i.id; }
};

//. \clsnm{Identity} gives derivative objects the ability to have
//. a unique identity and ordering relationship to all other
//. objects.
//. Normally \clsnm{Identity} must be inherited from virtually if multiple
//. inheritance is to be used.  This breaks certain compilers so
//. \srccd{NO\_VIRTUAL\_BASES} must be defined in these cases.  Not
//. everything will work under these circumstances.
class Identity {
  public:
    virtual ~Identity();
    //. Return the \clsnmref{Identifier} for this argument.
    //. Usually this is just the pointer to the object.
    Identifier identifier() { return this; }
};
//. Shorthand boolean operation for pointer arguments
inline int lt(const Identity*i, const Identity*j) { return i < j; }
inline int gt(const Identity*i, const Identity*j) { return i > j; }
inline int le(const Identity*i, const Identity*j) { return i <= j; }
inline int ge(const Identity*i, const Identity*j) { return i >= j; }
inline int eq(const Identity*i, const Identity*j) { return i == j; }
inline int ne(const Identity*i, const Identity*j) { return i != j; }
inline int cmp(const Identity*i, const Identity*j)
{
  return (i==j)?0:((i<j)?-1:1);
}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:

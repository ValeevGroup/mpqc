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

#include <iostream>

#include <scconfig.h>

namespace sc {

class Identity;

/** Identifier's are used to distinguish and order objects.  On many
    architectures a pointer to the object will suffice, but the C++
    standard only guarantees that this works for two pointers pointing
    within the same structure or array.  Classes need to inherit from
    Identity to use this mechanism.  Identity, Identifier, and the
    shorthand boolean operations may have to be modified for certain
    architectures.  */
class Identifier {
  private:
    const void* id;
  public:
    /// Create an Identifier for a null object.
    Identifier(): id(0) {}
    /// Create an Identifier for the given object.
    Identifier(const Identity* i): id((void*)i) {}
    /// Create an Identifier for the given object.
    Identifier(const Identifier& i): id(i.id) {}
    /// The destructor does nothing.
    ~Identifier() {}

    /// Assign to the given Identifier.
    void operator = (const Identifier& i) { id = i.id; }

    /// Less than.
    int operator < (const Identifier&i) const { return id < i.id; }
    /// Greater than.
    int operator > (const Identifier&i) const { return id > i.id; }
    /// Equal.
    int operator == (const Identifier&i) const { return id == i.id; }
    /// Less than or equal.
    int operator <= (const Identifier&i) const { return id <= i.id; }
    /// Greater than or equal.
    int operator >= (const Identifier&i) const { return id >= i.id; }
    /// Not equal.
    int operator != (const Identifier&i) const { return id != i.id; }

    void print(std::ostream&) const;
};

std::ostream & operator << (std::ostream &o, const Identifier &i);

/** Identity gives objects a unique identity and ordering relationship
 relative to all other objects.

 Identity must be virtually inherited if multiple inheritance is
 to be used. */
class Identity {
  public:
    virtual ~Identity();
    /** Return the Identifier for this argument.
        Usually this is just the pointer to the object. */
    Identifier identifier() { return this; }
};
/// Less than for two Identity pointers.
inline int lt(const Identity*i, const Identity*j) { return i < j; }
/// Greater than for two Identity pointers.
inline int gt(const Identity*i, const Identity*j) { return i > j; }
/// Less than or equal for two Identity pointers.
inline int le(const Identity*i, const Identity*j) { return i <= j; }
/// Greater than or equal for two Identity pointers.
inline int ge(const Identity*i, const Identity*j) { return i >= j; }
/// Equal for two Identity pointers.
inline int eq(const Identity*i, const Identity*j) { return i == j; }
/// Not equal for two Identity pointers.
inline int ne(const Identity*i, const Identity*j) { return i != j; }
/** Compare for two Identity pointers.  Returns -1, 0, or 1, like
    the C library function strcmp. */
inline int cmp(const Identity*i, const Identity*j)
{
  return (i==j)?0:((i<j)?-1:1);
}

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:

//
// units.h
//
// Copyright (C) 1997 Limit Point Systems, Inc.
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

#ifndef _util_misc_units_h
#define _util_misc_units_h

#include <util/class/class.h>
#include <util/state/state.h>
#include <util/ref/ref.h>

//. The \clsnm{Units} class is used to perform unit converions.
class Units: public SavableState {
#define CLASSNAME Units
#define HAVE_STATEIN_CTOR
#include <util/state/stated.h>
#include <util/class/classd.h>
  protected:
    char *strrep_;
    double to_atomic_units_;

    void parse_unit();
  public:
    enum Storage { Steal, Copy };

    Units(const char *strrep);
    Units(char *strrep, Units::Storage = Units::Copy);
    Units(StateIn&);
    ~Units();

    double to_atomic_units() const;
    double from_atomic_units() const;

    const char *string_rep() const;

    void save_data_state(StateOut&);
};
SavableState_REF_dec(Units);

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:

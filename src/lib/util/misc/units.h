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

#ifndef _util_misc_units_h
#define _util_misc_units_h

#include <util/class/class.h>
#include <util/state/state.h>
#include <util/ref/ref.h>

namespace sc {

/// The Units class is used to perform unit conversions.
class Units: public SavableState {
  protected:
    char *strrep_;
    double to_atomic_units_;

    void parse_unit();
  public:
    enum Storage { Steal, Copy };

    /// Create using a string representation, like "kcal/mol".
    Units(const char *strrep);
    /** Create using a string representation, like "kcal/mol".
        if Units::Steal is given is the second argment, the new
        Units object will delete the strrep argument when it is
        destroyed. */
    Units(char *strrep, Units::Storage = Units::Copy);
    /// Restore the state of a Units object from s.
    Units(StateIn& s);
    ~Units();

    /// The conversion factor from this to u.
    double to(const Ref<Units> &u) const;
    /// The conversion factor from u to this.
    double from(const Ref<Units> &u) const;

    /// The conversion factor from this to atomic units.
    double to_atomic_units() const;
    /// The conversion factor from atom units to this.
    double from_atomic_units() const;

    /// The string representation of the units.
    const char *string_rep() const;

    /// Save the state of the Units object to s.
    void save_data_state(StateOut&s);
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:

//
// atominfo.h
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

#ifndef _chemistry_molecule_atominfo_h
#define _chemistry_molecule_atominfo_h

#include <util/class/class.h>
#include <util/keyval/keyval.h>

class Units;

/** The AtomInfo class provides information about atoms.  The information
    is kept in a file named atominfo.kv in the SC library directory.  That
    information can be overridden by the user. */
class AtomInfo: public SavableState {
  private:
    enum { MaxZ = 107+1 };

    struct atomname
    {
      char *name;
      char *symbol;
    };

    static struct atomname names_[MaxZ];
    double  mass_[MaxZ];
    double  atomic_radius_[MaxZ];
    double  vdw_radius_[MaxZ];
    double  bragg_radius_[MaxZ];
    double  maxprob_radius_[MaxZ];
    double  rgb_[MaxZ][3];
    double  ip_[MaxZ];
    double  atomic_radius_scale_;
    double  vdw_radius_scale_;
    double  bragg_radius_scale_;
    double  maxprob_radius_scale_;

    char *overridden_values_;

    void load_library_values();
    void override_library_values(const Ref<KeyVal> &keyval);
    void load_values(const Ref<KeyVal>& keyval, int override);
    void load_values(double *array, double *scale, const char *keyword,
                     const Ref<KeyVal> &keyval, int override,
                     const Ref<Units> &);
    void load_values(double array[][3], const char *keyword,
                     const Ref<KeyVal> &keyval, int override);
    void add_overridden_value(const char *assignment);
  public:
    AtomInfo();
    AtomInfo(const Ref<KeyVal>&);
    AtomInfo(StateIn&);
    ~AtomInfo();
    void save_data_state(StateOut& s);

    /// These return various measures of the atom's radius.
    double vdw_radius(int Z) const { return vdw_radius_[Z]*vdw_radius_scale_; }
    double bragg_radius(int Z) const { return bragg_radius_[Z]*bragg_radius_scale_; }
    double atomic_radius(int Z) const { return atomic_radius_[Z]*atomic_radius_scale_; }
    double maxprob_radius(int Z) const { return maxprob_radius_[Z]*maxprob_radius_scale_; }

    /// Returns the atomization potential for atomic number Z.
    double ip(int Z) const { return ip_[Z]; }

    /// Return the scale factor for the VdW radii.
    double vdw_radius_scale() const { return vdw_radius_scale_; }
    /// Return the scale factor for the Bragg radii.
    double bragg_radius_scale() const { return bragg_radius_scale_; }
    /// Return the scale factor for the atomic radii.
    double atomic_radius_scale() const { return atomic_radius_scale_; }
    /// Return the scale factor for the maximum probability radii.
    double maxprob_radius_scale() const { return maxprob_radius_scale_; }

    /** These return information about the color of the atom
        for visualization programs. */
    double rgb(int Z, int color) const { return rgb_[Z][color]; }
    double red(int Z) const { return rgb_[Z][0]; }
    double green(int Z) const { return rgb_[Z][1]; }
    double blue(int Z) const { return rgb_[Z][2]; }

    /// This returns the mass of the most abundant isotope.
    double mass(int Z) const { return mass_[Z]; }

    /// This returns the full name of the element.
    static const char *name(int Z) { return names_[Z].name; }
    /// This returns the symbol for the element.
    static const char *symbol(int Z) { return names_[Z].symbol; }

    /// This converts a name or symbol to the atomic number.
    static int string_to_Z(const char *, int allow_exceptions = 1);
};


#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:

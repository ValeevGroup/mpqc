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

#include <string>
#include <map>
#include <vector>

#include <util/class/class.h>
#include <util/keyval/keyval.h>

/// @defgroup Chemistry mpqc.Chemistry
/// Classes/functions related to the chemistry domain

/// @defgroup ChemistryMolecule mpqc.Chemistry.Molecule
/// Classes/functions related to Molecule class

namespace sc {

  class Units;

  /// @addtogroup ChemistryMolecule
  /// @{

/** The AtomInfo class provides information about atoms.  The information
    is kept in a file named atominfo.kv in the SC library directory.  That
    information can be overridden by the user. */
class AtomInfo: public SavableState {
  private:
    enum { Nelement = 118, DefaultZ = 0 };

    struct atom
    {
      int Z;
      const char *name;
      const char *symbol;
    };

    static struct atom elements_[Nelement];

    std::map<std::string,int> name_to_Z_;
    std::map<std::string,int> symbol_to_Z_;
    std::map<int,std::string> Z_to_names_;
    std::map<int,std::string> Z_to_symbols_;
    std::map<int,double> Z_to_mass_;
    std::map<int,double> Z_to_atomic_radius_;
    std::map<int,double> Z_to_vdw_radius_;
    std::map<int,double> Z_to_bragg_radius_;
    std::map<int,double> Z_to_maxprob_radius_;
    std::map<int,std::vector<double> > Z_to_rgb_;
    std::map<int,double> Z_to_ip_;
    double  atomic_radius_scale_;
    double  vdw_radius_scale_;
    double  bragg_radius_scale_;
    double  maxprob_radius_scale_;

    char *overridden_values_;

    void load_library_values();
    void override_library_values(const Ref<KeyVal> &keyval);
    void load_values(const Ref<KeyVal>& keyval, int override);
    void load_values(std::map<int,double>&,
                     double *scale, const char *keyword,
                     const Ref<KeyVal> &keyval, int override,
                     const Ref<Units> &);
    void load_values(std::map<int,std::vector<double> >&,
                     const char *keyword,
                     const Ref<KeyVal> &keyval, int override);
    void add_overridden_value(const char *assignment);
    void initialize_names();
    double lookup_value(const std::map<int,double>& values, int Z) const;
    double lookup_array_value(const std::map<int,std::vector<double> >& values,
                              int Z, int i) const;

    // set to true after load_library_values() has been called
    // to avoid repeatedly printing out "Reading file ..." message
    static bool has_announced_library_source_;

  public:
    AtomInfo();

    /** The AtomInfo KeyVal constructor is used to generate a AtomInfo
        object from the input.  Default values will be read in from the
        <tt>atominfo.kv</tt> file in library directory.  These can be
        overridden by specifying the keyword below.  The library file is
        also read using a KeyVal constructor syntax, so consult that file
        for an example.

        <table border="1">

        <tr><td>%Keyword<td>Type<td>Default<td>Description

        <tr><td><tt>mass:unit</tt><td>string<td><tt>amu</tt><td>The unit to
        be used for masses.  See the Units class for more information about
        units.

        <tr><td><tt>mass:</tt><em>symbol</em><td>double<td>library
        value<td>The mass associated with the given atomic symbol.

        <tr><td><tt>vdw_radius:unit</tt><td>string<td><tt>bohr</tt><td>The
        unit to be used for van der Waals radii.  See the Units class for
        more information about units.

        <tr><td><tt>vdw_radius:scaling_factor</tt><td>double<td>1.0<td>The
        scaling factor to be used for all van der Waals radii, including
        library values.

        <tr><td><tt>vdw_radius:</tt><em>symbol</em><td>double<td>library
        value <td>The van der Waals radius associated with the given atomic
        symbol.

        <tr><td><tt>atomic_radius:unit</tt><td>string<td><tt>bohr</tt><td>The
        unit to be used for atomic radii.  See the Units class for more
        information about units.

        <tr><td><tt>atomic_radius:scaling_factor</tt><td>double<td>1.0<td>The
        scaling factor to be used for all atomic radii, including library
        values.

        <tr><td><tt>atomic_radius:</tt><em>symbol</em><td>double<td>library
        value <td>The atomic radius associated with the given atomic
        symbol.

        <tr><td><tt>bragg_radius:unit</tt><td>string<td><tt>bohr</tt><td>The
        unit to be used for Bragg radii.  See the Units class for more
        information about units.

        <tr><td><tt>bragg_radius:scaling_factor</tt><td>double<td>1.0<td>The
        scaling factor to be used for all Bragg radii, including library
        values.

        <tr><td><tt>bragg_radius:</tt><em>symbol</em><td>double<td>library
        value <td>The Bragg radius associated with the given atomic symbol.

        <tr><td><tt>maxprob_radius:unit</tt><td>string<td><tt>bohr</tt><td>The
        unit to be used for maximum probability radii.  See the Units class
        for more information about units.

        <tr><td><tt>maxprob_radius:scaling_factor</tt><td>double<td>1.0<td>The
        scaling factor to be used for all maximum probability radii,
        including library values.

        <tr><td><tt>maxprob_radius:</tt><em>symbol</em><td>double<td>library
        value<td>The maximum probability radius associated with the given
        atomic symbol.

        <tr><td><tt>ip:unit</tt><td>string<td><tt>Hartree</tt><td>The unit
        to be used for ionization potentials.  See the Units class for more
        information about units.

        <tr><td><tt>ip:</tt><em>symbol</em><td>double<td>library
        value<td>The ionization potential for the given atom.

        <tr><td><tt>rgb:</tt><em>symbol</em><td>double[3]<td>library
        value<td>A vector with the red, green, and blue values used to
        color each atom.  Each element is between 0 (off) and 1 (on).

        </table>
    */

    AtomInfo(const Ref<KeyVal>&);
    AtomInfo(StateIn&);
    ~AtomInfo();
    void save_data_state(StateOut& s);

    /// These return various measures of the atom's radius.
    double vdw_radius(int Z) const;
    double bragg_radius(int Z) const;
    double atomic_radius(int Z) const;
    double maxprob_radius(int Z) const;

    /// Returns the atomization potential for atomic number Z.
    double ip(int Z) const;

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
    double rgb(int Z, int color) const;
    double red(int Z) const;
    double green(int Z) const;
    double blue(int Z) const;

    /// This returns the mass of the most abundant isotope.
    double mass(int Z) const;

    /// This returns the full name of the element.
    std::string name(int Z) const;
    /// This returns the symbol for the element.
    std::string symbol(int Z) const;

    /// This converts a name or symbol to the atomic number.
    int string_to_Z(const std::string &, int allow_exceptions = 1);

    /// prints out the contents of AtomInfo to ostream os
    void print(std::ostream& os = ExEnv::out0()) const;
};

/// @}
// end of addtogroup ChemistryMolecule

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:

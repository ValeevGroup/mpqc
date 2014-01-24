//
// atom.h
//
// Copyright (C) 2013 Drew Lewis
//
// Author: Drew Lewis <drew90@vt.edu>
// Maintainer: Drew Lewis and Edward Valeev
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

#ifndef _chemisty_molecule_atom_h
#define _chemisty_molecule_atom_h

#include <string>
#include <util/state/stateio.h>
#include <util/misc/xml.h>

using boost::property_tree::ptree;

namespace sc {

  /// @addtogroup ChemistryMolecule
  /// @{

    /**
     * Atom represents an atom in a Molecule. Its position is defined relative to the
     * coordinate system of the Molecule. It also has an atomic number. Optionally, it can
     * has a charge, mass, and it can belong to a fragment.
     */
    class Atom : public XMLWritable {

      private:
        /// Contains the vector to the atom in units determined by molecule.
        double r_[3];
        int Z_;

        bool have_charge_;
        bool have_fragment_;

        // optional prameters.
        double charge_;
        int fragment_;
        double mass_;
        std::string label_;

      public:
        /**
         * Creates an atom for use in the sc::Molecule class.
         *
         * @param Z Atomic charge of the nucleous.
         * @param x x Coordinate of position vector.
         * @param y y Coordinate of position vector.
         * @param z z Coordinate of position vector.
         * @param label A string which gives the atom a user defined label.
         * @param mass The mass of the nuclous.
         * @param have_charge A boolean which signals a charge different than Z
         *      has been specified by the user.
         * @param charge A user specified charge which is different than Z.
         * @param have_fragment_ A boolean which specifies that the atom belongs
         *      to a specific fragment.
         * @param fragment Specifies to which fragment the atom belongs if
         *      0 then the atom does not belong to a fragment.
         */
        Atom(int Z, double x, double y, double z, const std::string &label,
             double mass = 0, int have_charge = 0, double charge = 0,
             int have_fragment = 0, int fragment = 0)
             : Z_(Z), have_charge_(have_charge), have_fragment_(have_fragment),
               charge_(charge), fragment_(fragment), mass_(mass), label_(label)
        {
            r_[0] = x;
            r_[1] = y;
            r_[2] = z;
        }

        /**
         * Constructs an atom object that takes a c style char * array
         */
        Atom(int Z, double x, double y, double z, const char *label = 0,
             double mass = 0, int have_charge = 0, double charge = 0,
             int have_fragment = 0, int fragment = 0)
             : Z_(Z), have_charge_(have_charge), have_fragment_(have_fragment),
               charge_(charge), fragment_(fragment), mass_(mass), label_(label ? label : "")
        {
            r_[0] = x;
            r_[1] = y;
            r_[2] = z;
        }


        /**
         * Default constructor supplied so that Atom will work with
         * sc::SavableState.  The user should not use this.
         *
         * @Warning Do not use, meant for sc::SavableState
         */
        Atom() : Z_(-1), have_charge_(true), have_fragment_(true),
                 charge_(-1), fragment_(-1), mass_(-1), label_()
        { r_[0] = -1; r_[1] = -1; r_[2] = -1; }

        Atom(const Atom &other) :
            Z_(other.Z_),
            have_charge_(other.have_charge_),
            have_fragment_(other.have_fragment_),
            charge_(other.charge_),
            fragment_(other.fragment_),
            mass_(other.mass_),
            label_(other.label_)
        {
            r_[0] = other.r_[0];
            r_[1] = other.r_[1];
            r_[2] = other.r_[2];
        }

        /// Returns a reference to the x,y, or z coordinate.
        double& xyz(int xyz){return r_[xyz];}
        const double& xyz(int xyz) const {return r_[xyz];}

        /// Returns a pointer to the coordinate array
        double* r(){return r_;}
        const double* r() const {return r_;}

        /// Returns atomic number
        int Z() const {return Z_;}
        double mass() const {return mass_;}
        bool have_charge() const {return have_charge_;}
        double charge() const {return charge_;}
        bool have_fragment() const {return have_fragment_;}
        int fragment() const {return fragment_;}
        const std::string& label() const {return label_;}

        // Made friend for direct access for sc::SavableState
        friend void FromStateIn(Atom &a, StateIn &so, int &count);

        virtual ptree& write_xml(ptree& parent, const XMLWriter& writer);

   };

    /// writes Atom to sc::StateOut
    void ToStateOut(const Atom &a, StateOut &so, int &count);

    /// reads Atom from sc::StateIn
    void FromStateIn(Atom &a, StateIn &si, int &count);

    /// @}
    // end of addtogroup ChemistryMolecule

}


#endif /* _chemistry_molecule_atom_h */

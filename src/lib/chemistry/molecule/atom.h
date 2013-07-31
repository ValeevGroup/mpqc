//
// atom.h
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Drew Lewis <drew90@vt.edu>
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

#ifndef _chemisty_molecule_atom_h
#define _chemisty_molecule_atom_h

#include <string>
#include <util/state/stateio.h>

namespace sc {

    /**
     * The atom class constains information about atoms. Its construction
     * follows the form of add_atom in molecule.  Atom was intended to be
     * a replacement for standalone arrays holding information about
     * atoms.
     *
     */
    class Atom {

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
        Atom(int Z, double x, double y, double z, char *label = 0,
             double mass = 0, int have_charge = false, double charge = 0,
             int have_fragment = false, int fragment = 0)
             : Z_(Z), mass_(mass), label_(label ? label : ""), charge_(charge),
               have_charge_(have_charge), have_fragment_(have_fragment),
               fragment_(fragment)
        {
            r_[0] = x;
            r_[1] = y;
            r_[2] = z;
        }

        Atom(int Z, double x, double y, double z, const std::string &label,
             double mass = 0, int have_charge = false, double charge = 0,
             int have_fragment = false, int fragment = 0)
             : Z_(Z), mass_(mass), label_(label), charge_(charge),
               have_charge_(have_charge), have_fragment_(have_fragment),
               fragment_(fragment)
        {
            r_[0] = x;
            r_[1] = y;
            r_[2] = z;
        }

        //Don't use this guy.
        Atom() : Z_(-1), mass_(-1), label_(), have_charge_(false), charge_(-1),
                        have_fragment_(false), fragment_(-1)
        { r_[0] = -1; r_[1] = -1; r_[2] = -1; }
        /** Default constructor supplied so that Atom will work with
         * sc::SavableState.  The user should not use this.
         */

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

   };

    /// Interface to sc::StateOut
    void ToStateOut(const Atom &a, StateOut &so, int &count);

    /// Interface to sc::StateIn
    void FromStateIn(Atom &a, StateIn &so, int &count);
}


#endif /* _chemistry_molecule_atom_h */

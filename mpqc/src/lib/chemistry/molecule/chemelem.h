//
// chemelem.h
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

#ifndef _chemistry_molecule_chemelem_h
#define _chemistry_molecule_chemelem_h

#ifdef __GNUC__
#pragma interface
#endif

// This class is a repository of element information
#include <util/class/class.h>
#include <util/state/state.h>

// Macro to expand to query function
#define QUERY_FUNCTION_IMPL(type,property) \
  INLINE type ChemicalElement::property() const \
  { return atom_info[Z_].property; }
#define QUERY_FUNCTION_CONV_IMPL(type,property,unit_conversion) \
  INLINE type ChemicalElement::property() const \
  { return atom_info[Z_].property * unit_conversion; }
							      
#define ANGSTROMS_TO_AU 1.0/0.52917706;

//.  The \clsnm{ChemicalElement} class provides information about
//individual elements.  Internally, a \clsnmref{ChemicalElement} only stores
//the atomic number.  From this, a table lookup can be done to return
//various properties, such as atomic masses, boiling points, etc.
//
//  \clsnm{ChemicalElement} is a \clsnmref{SavableState} so it has a
//\clsnmref{StateIn} constructor.  There is no \clsnmref{KeyVal}
//constructor.
class ChemicalElement: public SavableState
{
#   define CLASSNAME ChemicalElement
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  public:
    struct atom_info_type
    {
      char   *name;
      char   *symbol;
      int     number;
      double     mass;
      int     family;
      int     row;
      int     valence;
      int     melting_pt;
      int     boiling_pt;
      int     first_ion_pt;
      double  bragg_radius;
      double  electronegativity;
      double  specific_heat;
      double  density;
      double  atomic_radius;
      double  vdw_radius;
    };

  private:
    int Z_;
  
    static int max_atomic_number;
    static atom_info_type atom_info[];

  public:
    //. This constructor takes the atomic number as its argument.
    ChemicalElement(int Z=1);
    //. This constructor can take a string containing the atomic number
    // of the element, a string containing the full name of the element (can
    // be upper or lower case), or a string containing the atomic symbol (again
    // can be upper or lower case).
    ChemicalElement(const char* name);
    
    ChemicalElement(StateIn&);
    ChemicalElement(const ChemicalElement&);

    ~ChemicalElement();

    int operator == (const ChemicalElement &e) const { return Z_ == e.Z_; }
    
    void save_data_state(StateOut&);

    // Here are all of the query functions for the properties
    //. Returns full name of the element (all lower case).
    const char * name() const;
    //. Returns the atomic symbol for the element (mixed upper and lower case).
    const char * symbol() const;
    //. Returns the atomic number for the element.
    int number() const;
    //. Returns the mass of the most abundant isotope in AMU.
    double mass() const;
    //. Returns a number indicating which column of the periodic table the
    // element belongs to.  This will be 1 for the group IA or alkali metals,
    // up to 18 for Noble gases.
    int family() const;
    //. Returns the row of the periodic table the element belongs to.
    int row() const;
    //. Returns the number of electrons in the valence orbitals.
    int valence() const;
    //. Returns the melting point in degrees Celsius.
    int melting_pt() const;
    //. Returns the boiling point in degrees Celsius.
    int boiling_pt() const;
    //. Returns the first ionization potential in kcal/mol.
    int first_ion_pt() const;
    //. Returns some number Curt needs.
    double bragg_radius() const;
    //. Returns the electronegativity value for the element.
    double electronegativity() const;
    //. Returns the specific head of the element.
    double specific_heat() const;
    //. No clue.
    double density() const;
    //. Returns the atomic radius in bohr.
    double atomic_radius() const;
    //. Returns the van der Waals radius in bohr.
    double vdw_radius() const;
    //. Returns the nuclear charge.
    double charge() const;
};
DescribedClass_REF_dec(ChemicalElement);

#ifdef INLINE_FUNCTIONS
#include <chemistry/molecule/chmelm_i.h>
#endif

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:

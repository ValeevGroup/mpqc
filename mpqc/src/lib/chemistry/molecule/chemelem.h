
#ifndef _chemistry_molecule_chemelem_h
#define _chemistry_molecule_chemelem_h

#ifdef __GNUC__
#pragma interface
#endif

// This class is a repository of element information
#include <stdio.h>
#include <string.h>
#include <ctype.h>

#include <util/class/class.h>
#include <util/state/state.h>

// Macro to expand to query function
#define QUERY_FUNCTION_IMPL(type,property) \
  INLINE const type ChemicalElement::property() const \
  { return atom_info[Z_].property; }
#define QUERY_FUNCTION_CONV_IMPL(type,property,unit_conversion) \
  INLINE const type ChemicalElement::property() const \
  { return atom_info[Z_].property * unit_conversion; }
							      
#define ANGSTROMS_TO_AU 1.0/0.52917706;

//texi
// The @code{ChemicalElement} class provides information about individual
// elements.  Internally, a @code{ChemicalElement} only stores the atomic
// number.  From this, a table lookup can be done to return various properties,
// such as atomic masses, boiling points, etc.
//
// @code{ChemicalElement} is a @code{SavableState} so it has a @code{StateIn}
// constructor.  There is no @code{KeyVal} constructor.
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
    //texi This constructor takes the atomic number as its argument.
    ChemicalElement(int Z=1);
    //texi This constructor can take a string containing the atomic number
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
    //texi Returns full name of the element (all lower case).
    const char * name() const;
    //texi Returns the atomic symbol for the element (mixed upper and lower case).
    const char * symbol() const;
    //texi Returns the atomic number for the element.
    const int number() const;
    //texi Returns the mass of the most abundant isotope in AMU.
    const double mass() const;
    //texi Returns a number indicating which column of the periodic table the
    // element belongs to.  This will be 1 for the group IA or alkali metals,
    // up to 18 for Noble gases.
    const int family() const;
    //texi Returns the row of the periodic table the element belongs to.
    const int row() const;
    //texi Returns the number of electrons in the valence orbitals.
    const int valence() const;
    //texi Returns the melting point in degrees Celsius.
    const int melting_pt() const;
    //texi Returns the boiling point in degrees Celsius.
    const int boiling_pt() const;
    //texi Returns the first ionization potential in kcal/mol.
    const int first_ion_pt() const;
    //texi Returns some number Curt needs.
    const double bragg_radius() const;
    //texi Returns the electronegativity value for the element.
    const double electronegativity() const;
    //texi Returns the specific head of the element.
    const double specific_heat() const;
    //texi No clue.
    const double density() const;
    //texi Returns the atomic radius in bohr.
    const double atomic_radius() const;
    //texi Returns the van der Waals radius in bohr.
    const double vdw_radius() const;
    //texi Returns the nuclear charge.
    double charge() const;
};
DescribedClass_REF_dec(ChemicalElement);

#ifdef INLINE_FUNCTIONS
#include <chemistry/molecule/chmelm_i.h>
#endif

#endif

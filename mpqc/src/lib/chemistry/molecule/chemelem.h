
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
#define QUERY_FUNCTION_PROTO(type,property) \
  const type property() const
#define QUERY_FUNCTION_IMPL(type,property) \
  INLINE const type ChemicalElement::property() const \
  { return atom_info[Z_].property; }
#define QUERY_FUNCTION_CONV_IMPL(type,property,unit_conversion) \
  INLINE const type ChemicalElement::property() const \
  { return atom_info[Z_].property * unit_conversion; }
							      
#define ANGSTROMS_TO_AU 1.0/0.52917706;

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
    ChemicalElement(int Z=1);
    ChemicalElement(const ChemicalElement&);
    ChemicalElement(const char* name);
    ChemicalElement(StateIn&);

    ~ChemicalElement();
    
    void save_data_state(StateOut&);

    // Here are all of the query functions for the properties
    QUERY_FUNCTION_PROTO(char* ,name);
    QUERY_FUNCTION_PROTO(char* ,symbol);
    QUERY_FUNCTION_PROTO(int   ,number);
    QUERY_FUNCTION_PROTO(double,mass);
    QUERY_FUNCTION_PROTO(int   ,family);
    QUERY_FUNCTION_PROTO(int   ,row);
    QUERY_FUNCTION_PROTO(int   ,valence);
    QUERY_FUNCTION_PROTO(int   ,melting_pt);
    QUERY_FUNCTION_PROTO(int   ,boiling_pt);
    QUERY_FUNCTION_PROTO(int   ,first_ion_pt);
    QUERY_FUNCTION_PROTO(double,bragg_radius);
    QUERY_FUNCTION_PROTO(double,electronegativity);
    QUERY_FUNCTION_PROTO(double,specific_heat);
    QUERY_FUNCTION_PROTO(double,density);
    QUERY_FUNCTION_PROTO(double,atomic_radius);
    QUERY_FUNCTION_PROTO(double,vdw_radius);
    double charge() const;
};
DescribedClass_REF_dec(ChemicalElement);

#ifdef INLINE_FUNCTIONS
#include <chemistry/molecule/chmelm_i.h>
#endif

#endif

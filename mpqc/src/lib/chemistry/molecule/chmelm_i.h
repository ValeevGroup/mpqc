
#ifdef INLINE_FUNCTIONS
#define INLINE inline
#else
#define INLINE
#endif

QUERY_FUNCTION_IMPL(char* ,name);
QUERY_FUNCTION_IMPL(char* ,symbol);
QUERY_FUNCTION_IMPL(int   ,number);
QUERY_FUNCTION_IMPL(double,mass);
QUERY_FUNCTION_IMPL(int   ,family);
QUERY_FUNCTION_IMPL(int   ,row);
QUERY_FUNCTION_IMPL(int   ,valence);
QUERY_FUNCTION_IMPL(int   ,melting_pt);
QUERY_FUNCTION_IMPL(int   ,boiling_pt);
QUERY_FUNCTION_IMPL(int   ,first_ion_pt);
QUERY_FUNCTION_IMPL(double,bragg_radius);
QUERY_FUNCTION_IMPL(double,electronegativity);
QUERY_FUNCTION_IMPL(double,specific_heat);
QUERY_FUNCTION_IMPL(double,density);
QUERY_FUNCTION_IMPL(double,atomic_radius);
QUERY_FUNCTION_IMPL(double,vdw_radius);

INLINE double ChemicalElement::atomic_radius_au() const
{
  return atom_info[Z_].atomic_radius * ang_to_au;
}

INLINE double ChemicalElement::charge() const
{
  return (double) number();
}

#undef INLINE


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
QUERY_FUNCTION_CONV_IMPL(double,bragg_radius,ang_to_au);
QUERY_FUNCTION_IMPL(double,electronegativity);
QUERY_FUNCTION_IMPL(double,specific_heat);
QUERY_FUNCTION_IMPL(double,density);
QUERY_FUNCTION_CONV_IMPL(double,atomic_radius,ang_to_au);
QUERY_FUNCTION_CONV_IMPL(double,vdw_radius,ang_to_au);

INLINE double
ChemicalElement::charge() const
{
  return (double) number();
}

#undef INLINE


char *
compile_date()
{
#ifdef NCUBE_V2
  static char *date;
  date = "No date available on NCUBE_V2";
#else
  static char *date = COMPILE_DATE;
#endif
  return date;
  }

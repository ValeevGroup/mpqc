
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <strings.h>
#include <errno.h>

#include <fstream.h>

#include "mkclasses.h"

#ifdef __GNUG__
// forgive me
#define INSTANTIATE_VECTOR_OPT(theclass)		\
  template class vector<theclass>;		\
  template class vector_iterator<theclass>;	\
  template class vector_const_iterator<theclass>; \
  template void distance(const theclass *, const theclass *, unsigned int &); \
  template void __distance(const theclass *, const theclass *, \
                           unsigned int &, random_access_iterator_tag); \
  template void destroy(theclass *, theclass *); \
  template void destroy(theclass *);		\
  template void deallocate(theclass *);		\
  template void construct(theclass *, const theclass &); \
  template void iterator_category(const theclass *); \
  template void uninitialized_copy(theclass *, theclass *, theclass *); \
  template void uninitialized_copy(const theclass *, \
                                   const theclass *, theclass *); \
  template void copy(const theclass *, const theclass *, theclass *); \
  template void copy(theclass *, theclass *, theclass *); \
  template void copy_backward(theclass *, theclass *, theclass *); \
  template void swap(theclass *&, theclass *&);	\
  template void uninitialized_fill_n(theclass *, \
                                     unsigned int, const theclass &); \
  template void fill(theclass *, theclass *, const theclass &); \
theclass* allocate(int size, theclass*) {	\
    set_new_handler(0);				\
    theclass* tmp = (theclass*)			\
                    (::operator new((unsigned int)(size * sizeof(theclass)))); \
    if (tmp == 0) {				\
        cerr << "out of memory" << endl;	\
        exit(1);				\
    }						\
    return tmp;					\
}

#define INSTANTIATE_VECTOR_DEBUG(theclass) \
  template class allocator<theclass>; \
  template class reverse_iterator<const theclass *,theclass, \
                                  theclass&,int>; \
  template class reverse_iterator<theclass *,theclass, \
                                  theclass&,int>; \
  template class reverse_iterator<const theclass *,theclass, \
                                  const theclass&,int>;

#define INSTANTIATE_VECTOR(theclass) \
  INSTANTIATE_VECTOR_OPT(theclass) \
  INSTANTIATE_VECTOR_DEBUG(theclass)

template const unsigned int &max(unsigned int const &, unsigned int const &);

//  template theclass * allocate(long, theclass *);

INSTANTIATE_VECTOR(string)
INSTANTIATE_VECTOR(Parent)
INSTANTIATE_VECTOR(DataMember)
INSTANTIATE_VECTOR(Class)

#endif

int
main(int argc, char** argv)
{
  MkClassesOptions options;

  // process the command line
  int i;
  char *errormsg = 0;
  for (i=1; i<argc; i++) {
      if (!strcmp(argv[i], "--dbname")) {
          if (++i == argc) { errormsg = "--dbname requires argument"; break; }
          options.set_dbname(argv[i]);
        }
      else if (!strcmp(argv[i], "--libname")) {
          if (++i == argc) { errormsg = "--libname requires argument"; break; }
          options.set_libname(argv[i]);
        }
      else if (!strcmp(argv[i], "--gensrc")) {
          options.set_gensrc(true);
        }
      else if (!strcmp(argv[i], "--geninc")) {
          options.set_geninc(true);
        }
      else {
          ifstream input(argv[i]);
          options.set_hdrname(to_outfilename(argv[i], ".h"));
          options.set_srcname(to_outfilename(argv[i], ".cc"));
          MkClasses mk(options);
          mk.process(input);
        }
    }

  if (errormsg) {
      fprintf(stderr, "%s: %s\n", argv[0], errormsg);
      exit(1);
    }

  return 0;
}

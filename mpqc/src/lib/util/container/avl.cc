
#include <util/container/avl.h>

#ifdef __GNUG__
#define INST_COMP(T) \
  template int compare(const T &, const T &)

INST_COMP(int);
INST_COMP(long);
INST_COMP(double);
INST_COMP(char);
INST_COMP(unsigned char);
#endif

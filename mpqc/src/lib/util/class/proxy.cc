
#include <util/class/proxy.h>

#define CLASSNAME DescribedClassProxy
#define PARENTS public DescribedClass
#include <util/class/classia.h>
void *
DescribedClassProxy::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = DescribedClass::_castdown(cd);
  return do_castdowns(casts,cd);
}

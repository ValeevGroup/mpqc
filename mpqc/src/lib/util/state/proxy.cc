
#ifdef __GNUG__
#pragma implementation
#endif

#include <util/state/state_bin.h>
#include <util/state/proxy.h>
#include <util/keyval/keyval.h>

#define CLASSNAME SavableStateProxy
#define PARENTS public DescribedClassProxy
#define HAVE_KEYVAL_CTOR
#include <util/class/classi.h>
void *
SavableStateProxy::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] =  DescribedClassProxy::_castdown(cd) ;
  return do_castdowns(casts,cd);
}

SavableStateProxy::SavableStateProxy(const RefKeyVal &keyval)
{
  RefStateIn statein = keyval->describedclassvalue("statein");
  if (statein.nonnull()) {
      char *objectname = keyval->pcharvalue("object");
      StateIn &si = *(statein.pointer());
      if (keyval->exists("override")) {
          si.set_override(new PrefixKeyVal(keyval,"override"));
        }
      if (objectname) {
          object_.dir_restore_state(si, objectname);
          delete[] objectname;
        }
      else {
          object_.restore_state(si);
        }
    }
}

RefDescribedClass
SavableStateProxy::object()
{
  return object_;
}


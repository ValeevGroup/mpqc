
#include <Python.h>
#include <extern/python/util.h>
#include <boost/python.hpp>
#include <boost/python/implicit.hpp>
#include <util/state/state.h>

using namespace boost::python;
using namespace sc;

namespace sc {
  void
  export_util_state()
  {
    class_<SavableState, Ref<SavableState>, bases<DescribedClass>, boost::noncopyable >
      ("SavableState", no_init)
      ;
    implicitly_convertible<Ref<SavableState>, Ref<DescribedClass> >();
  }
}

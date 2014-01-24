

#include <Python.h>
#include <boost/python.hpp>
#include <util/ref/ref.h>
#include <extern/python/util.h>

using namespace boost::python;
using namespace sc;

namespace sc {
  void
  export_util_ref()
  {
    class_<RefCount, Ref<RefCount>, boost::noncopyable >
      ("RefCount", no_init)
      .def("nreference", &RefCount::nreference)
    ;
  }
}



#include <Python.h>
#include <boost/python.hpp>
#include <boost/python/implicit.hpp>
#include <util/class/class.h>
#include <extern/python/util.h>

using namespace boost::python;
using namespace sc;

namespace sc {

  void
  export_util_class()
  {
    void (DescribedClass::*p)(std::ostream&) const = &DescribedClass::print;
    class_<DescribedClass, bases<RefCount>, Ref<DescribedClass>, boost::noncopyable >
      ("DescribedClass", no_init)
      .def("class_name", &DescribedClass::class_name)
      .def("class_version", &DescribedClass::class_version)
      .def("ref", &DescribedClass::ref)
      .def("dcprint", p)
    ;
    implicitly_convertible<Ref<DescribedClass>, Ref<RefCount> >();
  }
}


#include <Python.h>
#include <boost/python.hpp>
#include <boost/python/implicit.hpp>
#include <util/misc/exenv.h>
#include <util/misc/compute.h>
#include <util/misc/units.h>
#include <extern/python/util.h>

using namespace boost::python;
using namespace sc;

namespace sc {
  void
  export_util_misc()
  {
    class_<std::ostream, boost::noncopyable> ("ostream", no_init);

    def("out0",&ExEnv::out0,
        return_value_policy<reference_existing_object>() )
      ;

    class_<Compute, boost::noncopyable >
      ("Compute", no_init)
      ;

    class_<Units, bases<SavableState>, Ref<Units>, boost::noncopyable >
      ("Units", init<const char*>())
      .def("to_atomic_units",&Units::to_atomic_units)
      .def("from_atomic_units",&Units::from_atomic_units)
      ;
    implicitly_convertible<Ref<Units>, Ref<SavableState> >();
  }
}

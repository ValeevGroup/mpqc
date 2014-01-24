

#include <Python.h>
#include <boost/python.hpp>
#include <boost/python/implicit.hpp>
#include <chemistry/qc/intv3/intv3.h>
#include <extern/python/util.h>

using namespace boost::python;
using namespace sc;

namespace sc {

  void
  export_chemistry_qc_intv3()
  {
    class_<IntegralV3, Ref<IntegralV3>, bases<Integral>, boost::noncopyable>
      ("IntegralV3", init<const Ref<KeyVal>&>())
      .def(init<const Ref<GaussianBasisSet> &,
           const Ref<GaussianBasisSet> &,
           const Ref<GaussianBasisSet> &,
           const Ref<GaussianBasisSet> &>())
      ;
    implicitly_convertible<Ref<IntegralV3>, Ref<Integral> >();
  }
}

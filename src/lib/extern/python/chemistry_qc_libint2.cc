

#include <Python.h>
#include <boost/python.hpp>
#include <boost/python/implicit.hpp>
#include <mpqc_config.h>
#include <chemistry/qc/libint2/libint2.h>
#include <extern/python/util.h>

using namespace boost::python;
using namespace sc;

namespace sc {

  void
  export_chemistry_qc_libint2()
  {
#ifdef HAVE_LIBINT2
    class_<IntegralLibint2, Ref<IntegralLibint2>, bases<Integral>, boost::noncopyable>
      ("IntegralLibint2", init<const Ref<KeyVal>&>())
      .def(init<const Ref<GaussianBasisSet> &,
           const Ref<GaussianBasisSet> &,
           const Ref<GaussianBasisSet> &,
           const Ref<GaussianBasisSet> &>())
      ;
    implicitly_convertible<Ref<IntegralLibint2>, Ref<Integral> >();
#endif
  }
}

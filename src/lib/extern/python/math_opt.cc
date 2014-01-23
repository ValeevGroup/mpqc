
#include <Python.h>
#include <boost/python.hpp>
#include <boost/python/implicit.hpp>
#include <math/optimize/opt.h>
#include <math/optimize/scextrap.h>
#include <math/optimize/diis.h>
#include <math/optimize/scextrapmat.h>
#include <extern/python/util.h>

using namespace boost::python;
using namespace sc;

namespace sc {

  void
  export_math_opt()
  {

    class_<Function, Ref<Function>, bases<SavableState,Compute>, boost::noncopyable >
      ("Function", no_init)
      .def("matrixkit",&Function::matrixkit)
      .def("dimension",&Function::dimension)
      .def("set_x",&Function::set_x)
      .def("get_x",&Function::get_x)
      .def("value",&Function::value)
      .def("gradient",&Function::gradient)
      .def("hessian",&Function::hessian)
    ;
    implicitly_convertible<Ref<Function>, Ref<SavableState> >();

    class_<SCExtrapData, Ref<SCExtrapData>, bases<SavableState>, boost::noncopyable>
      ("SCExtrapData", no_init)
      .def("zero",&SCExtrapData::zero)
      .def("accumulate_scaled",&SCExtrapData::accumulate_scaled)
      ;
    implicitly_convertible<Ref<SCExtrapData>, Ref<SavableState> >();

    class_<SCExtrapError, Ref<SCExtrapError>, bases<SavableState>, boost::noncopyable>
      ("SCExtrapError", no_init)
      .def("error",&SCExtrapError::error)
      .def("scalar_product",&SCExtrapError::scalar_product)
      ;
    implicitly_convertible<Ref<SCExtrapError>, Ref<SavableState> >();

    class_<SymmSCMatrixSCExtrapData, Ref<SymmSCMatrixSCExtrapData>, bases<SCExtrapData>, boost::noncopyable>
      ("SymmSCMatrixSCExtrapData", init<const RefSymmSCMatrix&>())
      ;
    implicitly_convertible<Ref<SymmSCMatrixSCExtrapData>, Ref<SCExtrapData> >();

    class_<SymmSCMatrixSCExtrapError, Ref<SymmSCMatrixSCExtrapError>, bases<SCExtrapError>, boost::noncopyable>
      ("SymmSCMatrixSCExtrapError", init<const RefSymmSCMatrix&>())
      ;
    implicitly_convertible<Ref<SymmSCMatrixSCExtrapError>, Ref<SCExtrapError> >();

    class_<SelfConsistentExtrapolation, Ref<SelfConsistentExtrapolation>, bases<SavableState>, boost::noncopyable>
      ("SelfConsistentExtrapolation", no_init)
      .def("set_tolerance",&SelfConsistentExtrapolation::set_tolerance)
      .def("tolerance",&SelfConsistentExtrapolation::tolerance)
      .def("error",&SelfConsistentExtrapolation::error)
      .def("converged",&SelfConsistentExtrapolation::converged)
      .def("extrapolate",&SelfConsistentExtrapolation::extrapolate)
      .def("start_extrapolation",&SelfConsistentExtrapolation::start_extrapolation)
      .def("reinitialize",&SelfConsistentExtrapolation::reinitialize)
      ;
    implicitly_convertible<Ref<SelfConsistentExtrapolation>, Ref<SavableState> >();

    class_<DIIS, Ref<DIIS>, bases<SelfConsistentExtrapolation>, boost::noncopyable>
      ("DIIS", init<const Ref<KeyVal>&>())
      ;
    implicitly_convertible<Ref<DIIS>, Ref<SelfConsistentExtrapolation> >();

  }
}

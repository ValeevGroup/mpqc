
#include <Python.h>
#include <boost/python.hpp>
#include <boost/python/implicit.hpp>
#include <chemistry/qc/wfn/obwfn.h>
#include <extern/python/util.h>

using namespace boost::python;
using namespace sc;

namespace sc {
  void
  export_chemistry_qc_wfn()
  {
    class_<Wavefunction, Ref<Wavefunction>, bases<MolecularEnergy>, boost::noncopyable >
      ("Wavefunction", no_init)
      .def("nelectron", &Wavefunction::nelectron)
      .def("magnetic_moment", &Wavefunction::magnetic_moment)
      ;
    implicitly_convertible<Ref<Wavefunction>, Ref<MolecularEnergy> >();

    class_<OneBodyWavefunction, Ref<OneBodyWavefunction>, bases<Wavefunction>, boost::noncopyable >
      ("OneBodyWavefunction", no_init)
      ;
    implicitly_convertible<Ref<OneBodyWavefunction>, Ref<Wavefunction> >();
  }
}

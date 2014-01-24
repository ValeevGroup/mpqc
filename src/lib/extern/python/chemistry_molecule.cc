

#include <Python.h>
#include <boost/python.hpp>
#include <boost/python/implicit.hpp>
#include <chemistry/molecule/energy.h>
#include <extern/python/util.h>

using namespace boost::python;
using namespace sc;

namespace sc {
  void
  export_chemistry_molecule()
  {
    void (Molecule::*aa)(int,double,double,double,
                         const std::string &, double mass,
                         int have_charge, double charge,
                         int have_fragment, int fragment) = &Molecule::add_atom;

    class_<Molecule, bases<SavableState>, Ref<Molecule>, boost::noncopyable >
      ("Molecule", init<>())
      .def(init<Ref<KeyVal> >())
      .def("add_atom",aa,(arg("Z"),arg("x"),arg("y"),arg("z"),arg("label")="",
                          arg("mass")=0.0,
                          arg("have_charge")=0, arg("charge")=0.0,
                          arg("have_fragment")=0, arg("fragment") = 0))
      .def("natom",&Molecule::natom)
      .def("Z",&Molecule::Z)
      .def("mass",&Molecule::mass)
      .def("label",&Molecule::label)
      .def("atom_at_position",&Molecule::atom_at_position)
      .def("atom_label_to_index",&Molecule::atom_label_to_index)
      .def("charge",&Molecule::charge)
      .def("total_charge",&Molecule::total_charge)
      .def("total_Z",&Molecule::total_Z)
      .def("nuclear_repulsion_energy",&Molecule::nuclear_repulsion_energy)
      .def("geometry_units",&Molecule::geometry_units)
      ;
    implicitly_convertible<Ref<Molecule>, Ref<SavableState> >();
  
    class_<MolecularEnergy, Ref<MolecularEnergy>, bases<Function>, boost::noncopyable >
      ("MolecularEnergy", no_init)
      .def("molecule", &MolecularEnergy::molecule)
      .def("energy", &MolecularEnergy::energy)
      ;
    implicitly_convertible<Ref<MolecularEnergy>, Ref<Function> >();

  }
}

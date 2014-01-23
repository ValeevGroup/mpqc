
#include <Python.h>
#include <boost/python.hpp>
#include <boost/python/implicit.hpp>
#include <chemistry/qc/scf/clhf.h>
#include <chemistry/qc/lcao/fockbuild.h>
#include <chemistry/qc/lcao/clhfcontrib.h>
#include <extern/python/util.h>

using namespace boost::python;
using namespace sc;

namespace sc {

  void
  export_chemistry_qc_scf()
  {
    class_<SCF, Ref<SCF>, bases<OneBodyWavefunction>, boost::noncopyable >
      ("SCF", no_init)
      ;
    implicitly_convertible<Ref<SCF>, Ref<OneBodyWavefunction> >();

    class_<CLSCF, Ref<CLSCF>, bases<SCF>, boost::noncopyable >
      ("CLSCF", no_init)
      ;
    implicitly_convertible<Ref<CLSCF>, Ref<SCF> >();

    class_<CLHF, Ref<CLHF>, bases<CLSCF>, boost::noncopyable >
      ("CLHF", init<Ref<KeyVal> >())
      ;
    implicitly_convertible<Ref<CLHF>, Ref<CLSCF> >();

    class_<FockBuild, Ref<FockBuild>, bases<RefCount>, boost::noncopyable >
      ("FockBuild", init<const Ref<FockDistribution> &, const Ref<FockContribution> &, bool, const Ref<GaussianBasisSet> &>())
      .def(init<const Ref<FockDistribution> &,
           const Ref<FockContribution> &,
           bool,
           const Ref<GaussianBasisSet> &,
           const Ref<GaussianBasisSet> &,
           const Ref<GaussianBasisSet> &>())
      .def(init<const Ref<FockDistribution> &,
           const Ref<FockContribution> &,
           bool,
           const Ref<GaussianBasisSet> &,
           const Ref<GaussianBasisSet> &,
           const Ref<GaussianBasisSet> &,
           const Ref<MessageGrp> &,
           const Ref<ThreadGrp> &,
           const Ref<Integral> &>())
      .def("build",&FockBuild::build)
      // this returns a const Ref<...>& and a policy must be specified:
      //.def("contrib",&FockBuild::contrib)
      .def("set_accuracy",&FockBuild::set_accuracy)
      .def("set_compute_J",&FockBuild::set_compute_J)
      .def("set_compute_K",&FockBuild::set_compute_K)
      .def("set_coef_K",&FockBuild::set_coef_K)
      .def("compute_J",&FockBuild::compute_J)
      .def("compute_K",&FockBuild::compute_K)
      .def("coef_K",&FockBuild::coef_K)
      ;
    implicitly_convertible<Ref<FockBuild>, Ref<RefCount> >();

    class_<FockContribution, Ref<FockContribution>, bases<RefCount>, boost::noncopyable >
      ("FockContribution", no_init)
      ;

    {
      void (GenericFockContribution::*set_fmat_1)(int i, const RefSCMatrix &)
        = &GenericFockContribution::set_fmat;
      void (GenericFockContribution::*set_fmat_2)(int i, const RefSymmSCMatrix &)
        = &GenericFockContribution::set_fmat;
      void (GenericFockContribution::*set_jmat_1)(int i, const RefSCMatrix &)
        = &GenericFockContribution::set_jmat;
      void (GenericFockContribution::*set_jmat_2)(int i, const RefSymmSCMatrix &)
        = &GenericFockContribution::set_jmat;
      void (GenericFockContribution::*set_kmat_1)(int i, const RefSCMatrix &)
        = &GenericFockContribution::set_kmat;
      void (GenericFockContribution::*set_kmat_2)(int i, const RefSymmSCMatrix &)
        = &GenericFockContribution::set_kmat;

      class_<GenericFockContribution, Ref<GenericFockContribution>, bases<FockContribution>, boost::noncopyable >
        ("GenericFockContribution", no_init)
        .def("set_fmat",set_fmat_1)
        .def("set_fmat",set_fmat_2)
        .def("set_jmat",set_jmat_1)
        .def("set_jmat",set_jmat_2)
        .def("set_kmat",set_kmat_1)
        .def("set_kmat",set_kmat_2)
        .def("set_pmat",&GenericFockContribution::set_pmat)
        ;
      implicitly_convertible<Ref<GenericFockContribution>, Ref<FockContribution> >();
    }

    class_<CLHFContribution, Ref<CLHFContribution>, bases<GenericFockContribution>, boost::noncopyable >
      ("CLHFContribution", init< Ref<GaussianBasisSet> &, Ref<GaussianBasisSet> &, Ref<GaussianBasisSet> &, const std::string & >())
      ;
    implicitly_convertible<Ref<CLHFContribution>, Ref<GenericFockContribution> >();

    class_<FockDistribution, Ref<FockDistribution>, bases<SavableState>, boost::noncopyable >
      ("FockDistribution", init<>())
      ;
    implicitly_convertible<Ref<FockDistribution>, Ref<SavableState> >();

  }
}

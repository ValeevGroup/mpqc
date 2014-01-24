
#include <Python.h>
#include <numpy/arrayobject.h>
#include <boost/python/dict.hpp>
#include <boost/python.hpp>
#include <boost/python/implicit.hpp>
#include <chemistry/qc/basis/basis.h>
#include <chemistry/qc/basis/orthog.h>
#include <chemistry/qc/basis/integral.h>
#include <chemistry/qc/basis/operator.h>
#include <chemistry/qc/basis/petite.h>
#include <chemistry/qc/basis/obint.h>
#include <chemistry/qc/basis/tbint.h>
#include <chemistry/qc/basis/symmint.h>
#include <extern/python/util.h>

using namespace boost::python;
using namespace sc;

namespace sc {

  struct convert_tbarr4_to_numpy
  {
    static PyObject* convert(const std::pair<std::map<TwoBodyOper::type,const double*>,std::array<unsigned long, 4> >& x)
    {
      npy_intp dim[4];
      PyObject *pydict = PyDict_New();
      for (int i=0; i<4; i++) dim[i] = x.second[i];
      for (std::map<TwoBodyOper::type,const double*>::const_iterator iter = x.first.begin();
           iter != x.first.end();
           iter++) {
        PyObject *arr = PyArray_New(&PyArray_Type, 4, dim, NPY_DOUBLE, NULL,
                                    (char*)iter->second, 0, NPY_CARRAY, NULL);
        object inttype(iter->first);
        PyDict_SetItem(pydict, inttype.ptr(), arr);
      }
      return pydict;
    }
  };

  void
  export_chemistry_qc_basis()
  {
    import_array();

    {
      int (OperatorDescr::*perm_symm_i)(unsigned int i) const
        = &OperatorDescr::perm_symm;
      int (OperatorDescr::*perm_symm_ij)(unsigned int i,
                                         unsigned int j) const
        = &OperatorDescr::perm_symm;
      class_<OperatorDescr, Ref<OperatorDescr>, bases<RefCount>,
        boost::noncopyable >
        ("OperatorDescr", no_init)
        .def("num_particles", &OperatorDescr::num_particles)
        .def("perm_symm_i", perm_symm_i)
        .def("perm_symm_ij", perm_symm_ij);
      implicitly_convertible<Ref<OperatorDescr>, Ref<RefCount> >();
    }

    {
      int (OneBodyOperDescr::*perm_symm_i)(unsigned int i) const
        = &OneBodyOperDescr::perm_symm;
      int (OneBodyOperDescr::*perm_symm_ij)(unsigned int i,
                                            unsigned int j) const
        = &OneBodyOperDescr::perm_symm;
      class_<OneBodyOperDescr, Ref<OneBodyOperDescr>, bases<OperatorDescr>,
        boost::noncopyable >
        ("OneBodyOperDescr", no_init)
        .def("num_particles", &OneBodyOperDescr::num_particles)
        .def("perm_symm_i", perm_symm_i)
        .def("perm_symm_ij", perm_symm_ij);
      implicitly_convertible<Ref<OneBodyOperDescr>, Ref<OperatorDescr> >();
    }

    {
      int (TwoBodyOperDescr::*perm_symm_i)(unsigned int i) const
        = &TwoBodyOperDescr::perm_symm;
      int (TwoBodyOperDescr::*perm_symm_ij)(unsigned int i,
                                            unsigned int j) const
        = &TwoBodyOperDescr::perm_symm;
      class_<TwoBodyOperDescr, Ref<TwoBodyOperDescr>, bases<OperatorDescr>,
        boost::noncopyable >
        ("TwoBodyOperDescr", no_init)
        .def("num_particles", &TwoBodyOperDescr::num_particles)
        .def("perm_symm_i", perm_symm_i)
        .def("perm_symm_ij", perm_symm_ij);
      implicitly_convertible<Ref<TwoBodyOperDescr>, Ref<OperatorDescr> >();
    }

    {
      scope inner = class_<OneBodyOper>("OneBodyOper", no_init)
        .def_readonly("max_ntypes", &OneBodyOper::max_ntypes)
        .def("descr", &OneBodyOper::descr)
        .staticmethod("descr")
        .def("to_string", &OneBodyOper::to_string)
        .staticmethod("to_string")
        ;

      enum_<OneBodyOper::type>("type")
        .value("T",OneBodyOper::T)
        .value("V",OneBodyOper::V)
        .value("h",OneBodyOper::h)
        ;
    }

    {
      scope inner = class_<TwoBodyOper>("TwoBodyOper", no_init)
        .def_readonly("max_ntypes", &TwoBodyOper::max_ntypes)
        .def("descr", &TwoBodyOper::descr)
        .staticmethod("descr")
        .def("to_string", &TwoBodyOper::to_string)
        .staticmethod("to_string")
        ;

      enum_<TwoBodyOper::type>("type")
        .value("eri",TwoBodyOper::eri)
        .value("r12",TwoBodyOper::r12)
        .value("r12t1",TwoBodyOper::r12t1)
        .value("r12t2",TwoBodyOper::r12t2)
        .value("r12_0_g12",TwoBodyOper::r12_0_g12)
        .value("r12_m1_g12",TwoBodyOper::r12_m1_g12)
        .value("t1g12",TwoBodyOper::t1g12)
        .value("t2g12",TwoBodyOper::t2g12)
        .value("g12t1g12",TwoBodyOper::g12t1g12)
        .value("g12p4g12_m_g12t1g12t1",TwoBodyOper::g12p4g12_m_g12t1g12t1)
        .value("anti_g12g12",TwoBodyOper::anti_g12g12)
        .value("delta",TwoBodyOper::delta)
        ;
    }
        
    {
      scope inner = class_<TwoBodyOperSet>("TwoBodyOperSet", no_init);
      enum_<TwoBodyOperSet::type>("type")
        .value("ERI", TwoBodyOperSet::ERI)
        .value("R12", TwoBodyOperSet::R12)
        .value("G12", TwoBodyOperSet::G12)
        .value("G12NC", TwoBodyOperSet::G12NC)
        .value("G12DKH", TwoBodyOperSet::G12DKH);
    }

    {
      class_<GaussianBasisSet, Ref<GaussianBasisSet>, bases<SavableState>,
        boost::noncopyable >
        ("GaussianBasisSet", init<Ref<KeyVal> >())
        .def("ncenter",&GaussianBasisSet::ncenter)
        .def("nshell",&GaussianBasisSet::nshell)
        .def("nshell_on_center",&GaussianBasisSet::nshell_on_center)
        .def("shell_on_center",&GaussianBasisSet::shell_on_center)
        .def("shell_to_center",&GaussianBasisSet::shell_to_center)
        .def("nbasis",&GaussianBasisSet::nbasis)
        .def("nbasis_on_center",&GaussianBasisSet::nbasis_on_center)
        .def("nprimitive",&GaussianBasisSet::nprimitive)
        .def("has_pure",&GaussianBasisSet::has_pure)
        .def("max_ncartesian_in_shell",&GaussianBasisSet::max_ncartesian_in_shell,
             (arg("aminc")=0))
        .def("max_nfunction_in_shell",&GaussianBasisSet::max_nfunction_in_shell)
        .def("max_nprimitive_in_shell",&GaussianBasisSet::max_nprimitive_in_shell)
        .def("max_angular_momentum",&GaussianBasisSet::max_angular_momentum)
        .def("max_ncontraction",&GaussianBasisSet::max_ncontraction)
        .def("max_am_for_contraction",&GaussianBasisSet::max_am_for_contraction)
        .def("max_cartesian",&GaussianBasisSet::max_cartesian)
        .def("shell_to_function",&GaussianBasisSet::shell_to_function)
        .def("function_to_shell",&GaussianBasisSet::function_to_shell)
        .def("basisdim",&GaussianBasisSet::basisdim)
        .def("matrixkit",&GaussianBasisSet::matrixkit)
        .def("so_matrixkit",&GaussianBasisSet::so_matrixkit)
        .def("molecule",&GaussianBasisSet::molecule)
        ;
    }

    {
      void (OverlapOrthog::*reinit)(OverlapOrthog::OrthogMethod,
                                    const RefSymmSCMatrix &,
                                    const Ref<SCMatrixKit> &,
                                    double, int)
        = &OverlapOrthog::reinit;
      scope inner =
        class_<OverlapOrthog, Ref<OverlapOrthog>, bases<SavableState>, boost::noncopyable>
        ("OverlapOrthog",init<OverlapOrthog::OrthogMethod,
         const RefSymmSCMatrix &, const Ref<SCMatrixKit> &, double, int>
         ((arg("method"),"overlap","result_kit","lindep_tolerance",arg("debug")=0)))
        .def("reinit",reinit)
        .def("min_orthog_res",&OverlapOrthog::min_orthog_res)
        .def("max_orthog_res",&OverlapOrthog::max_orthog_res)
        .def("copy",&OverlapOrthog::copy)
        .def("orthog_method",&OverlapOrthog::orthog_method)
        .def("lindep_tol",&OverlapOrthog::lindep_tol)
        .def("basis_to_orthog_basis",&OverlapOrthog::basis_to_orthog_basis)
        .def("basis_to_orthog_basis_inverse",&OverlapOrthog::basis_to_orthog_basis_inverse)
        .def("overlap_inverse",&OverlapOrthog::overlap_inverse)
        .def("dim",&OverlapOrthog::dim)
        .def("orthog_dim",&OverlapOrthog::orthog_dim)
        .def("nlindep",&OverlapOrthog::nlindep)
        ;

      enum_<OverlapOrthog::OrthogMethod>("OrthogMethod")
        .value("Symmetric",OverlapOrthog::Symmetric)
        .value("Canonical",OverlapOrthog::Canonical)
        .value("GramSchmidt",OverlapOrthog::GramSchmidt)
        ;

      implicitly_convertible<Ref<OverlapOrthog>, Ref<SavableState> >();

    }

    {
      int (PetiteList::*in_p2_1)(int) const = &PetiteList::in_p2;
      int (PetiteList::*in_p2_2)(int,int) const = &PetiteList::in_p2;
      int (PetiteList::*in_p4_4)(int,int,int,int) const = &PetiteList::in_p4;
      int (PetiteList::*in_p4_6)(int,int,int,int,int,int) const = &PetiteList::in_p4;
      class_<PetiteList, Ref<PetiteList>, bases<RefCount>, boost::noncopyable>
        ("PetiteList", init<const Ref<GaussianBasisSet>&, const Ref<Integral>&>())
        .def("basis",&PetiteList::basis)
        .def("integral",&PetiteList::integral)
        .def("clone",&PetiteList::clone)
        .def("nirrep",&PetiteList::nirrep)
        .def("order",&PetiteList::order)
        .def("AO_basisdim",&PetiteList::AO_basisdim)
        .def("SO_basisdim",&PetiteList::SO_basisdim)
        .def("aotoso",&PetiteList::aotoso)
        .def("sotoao",&PetiteList::sotoao)
        .def("symmetrize",&PetiteList::symmetrize)
        .def("to_SO_basis",&PetiteList::to_SO_basis)
        .def("to_AO_basis",&PetiteList::to_AO_basis)
        .def("evecs_to_SO_basis",&PetiteList::evecs_to_SO_basis)
        .def("evecs_to_AO_basis",&PetiteList::evecs_to_AO_basis)
        .def("in_p1",&PetiteList::in_p1)
        .def("in_p2_1",in_p2_1)
        .def("in_p2_2",in_p2_2)
        .def("in_p4_4",in_p4_4)
        .def("in_p4_6",in_p4_6)
        .def("nfunction",&PetiteList::nfunction)
        .def("nblocks",&PetiteList::nblocks)
        .def("r",&PetiteList::r)
        ;
      implicitly_convertible<Ref<PetiteList>, Ref<RefCount> >();
    }

    Ref<PetiteList> (Integral::*pl0)() = &Integral::petite_list;
    Ref<PetiteList> (Integral::*pl1)(const Ref<GaussianBasisSet>&)
      = &Integral::petite_list;
    void (Integral::*sb4)(const Ref<GaussianBasisSet>&,
                          const Ref<GaussianBasisSet>&,
                          const Ref<GaussianBasisSet>&,
                          const Ref<GaussianBasisSet>&)
      = &Integral::set_basis;
    class_<Integral, Ref<Integral>, bases<SavableState>, boost::noncopyable>
        ("Integral",no_init)
      .def("petite_list",pl0)
      .def("petite_list",pl1)
      .def("set_basis",sb4)
      .def("overlap",&Integral::overlap)
      .def("kinetic",&Integral::kinetic)
      .def("nuclear",&Integral::nuclear)
      .def("hcore",&Integral::hcore)
      .def("electron_repulsion",&Integral::electron_repulsion)
      .def("set_storage",&Integral::set_storage)
      .def("set_default_integral",&Integral::set_default_integral)
      .staticmethod("set_default_integral")
      ;
    implicitly_convertible<Ref<Integral>, Ref<SavableState> >();

    class_<OneBodyInt, Ref<OneBodyInt>, bases<RefCount>, boost::noncopyable>
      ("OneBodyInt", no_init)
      .def("nbasis",&OneBodyInt::nbasis)
      .def("nbasis1",&OneBodyInt::nbasis1)
      .def("nbasis2",&OneBodyInt::nbasis2)
      .def("nshell",&OneBodyInt::nshell)
      .def("nshell1",&OneBodyInt::nshell1)
      .def("nshell2",&OneBodyInt::nshell2)
      .def("basis",&OneBodyInt::basis)
      .def("basis1",&OneBodyInt::basis1)
      .def("basis2",&OneBodyInt::basis2)
      //.def("buffer",&OneBodyInt::buffer,return_internal_reference<>())
      .def("compute_shell",&OneBodyInt::compute_shell)
      .def("compute_shell_array",&OneBodyInt::compute_shell_array)
      .def("reinitialize",&OneBodyInt::reinitialize)
      ;
    implicitly_convertible<Ref<OneBodyInt>, Ref<RefCount> >();

    {
      scope inner
        = class_<TwoBodyInt, Ref<TwoBodyInt>, bases<RefCount>, boost::noncopyable>
        ("TwoBodyInt", no_init)
        .def("nbasis",&TwoBodyInt::nbasis)
        .def("nbasis1",&TwoBodyInt::nbasis1)
        .def("nbasis2",&TwoBodyInt::nbasis2)
        .def("nbasis3",&TwoBodyInt::nbasis3)
        .def("nbasis4",&TwoBodyInt::nbasis4)
        .def("nshell",&TwoBodyInt::nshell)
        .def("nshell1",&TwoBodyInt::nshell1)
        .def("nshell2",&TwoBodyInt::nshell2)
        .def("nshell3",&TwoBodyInt::nshell3)
        .def("nshell4",&TwoBodyInt::nshell4)
        .def("basis",&TwoBodyInt::basis)
        .def("basis1",&TwoBodyInt::basis1)
        .def("basis2",&TwoBodyInt::basis2)
        .def("basis3",&TwoBodyInt::basis3)
        .def("basis4",&TwoBodyInt::basis4)
        //.def("buffer",&TwoBodyInt::buffer,return_internal_reference<>())
        .def("compute_shell",&TwoBodyInt::compute_shell)
        .def("compute_shell_arrays",&TwoBodyInt::compute_shell_arrays)
        .def("type",&TwoBodyInt::type)
        ;
      implicitly_convertible<Ref<TwoBodyInt>, Ref<RefCount> >();
    }

    to_python_converter<std::pair<std::map<TwoBodyOper::type,const double*>,
      std::array<unsigned long,4> >,
      convert_tbarr4_to_numpy>();

    class_<OneBodyIntIter, Ref<OneBodyIntIter>, bases<RefCount>, boost::noncopyable >
      ("OneBodyIntIter",init<const Ref<OneBodyInt>&>())
       ;
    implicitly_convertible<Ref<OneBodyIntIter>, Ref<RefCount> >();

    class_<SymmOneBodyIntIter, Ref<SymmOneBodyIntIter>, bases<OneBodyIntIter>, boost::noncopyable>
       ("SymmOneBodyIntIter",init<const Ref<OneBodyInt>&,const Ref<PetiteList>&>())
        ;
    implicitly_convertible<Ref<SymmOneBodyIntIter>, Ref<OneBodyIntIter> >();

    class_<OneBodyIntOp, Ref<OneBodyIntOp>, bases<SCElementOp>, boost::noncopyable >
      ("OneBodyIntOp",init<const Ref<OneBodyInt>&>())
       .def(init<const Ref<OneBodyIntIter>&>())
       ;
    implicitly_convertible<Ref<OneBodyIntOp>, Ref<SCElementOp> >();

  }
}

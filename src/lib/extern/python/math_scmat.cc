

#include <Python.h>
#include <boost/python.hpp>
#include <boost/python/implicit.hpp>
#include <math/scmat/matrix.h>
#include <math/scmat/elemop.h>
#include <extern/python/util.h>

using namespace boost::python;
using namespace sc;

namespace sc {

  void
  export_math_scmat()
  {
    class_<RefSCDimension,bases<SavableState> >
      ("SCDimension", no_init)
      .def("n",&RefSCDimension::n)
      ;

    class_<SCMatrixKit,Ref<SCMatrixKit>,bases<DescribedClass>, boost::noncopyable >
      ("SCMatrixKit", no_init)
      ;
    implicitly_convertible<Ref<SCMatrixKit>, Ref<DescribedClass> >();

    class_<SCElementOp,Ref<SCElementOp>,bases<SavableState>(),boost::noncopyable>
      ("ElementOp", no_init)
      ;
    implicitly_convertible<Ref<SCElementOp>, Ref<SavableState> >();

    class_<SCElementOp2,Ref<SCElementOp2>,bases<SavableState>(),boost::noncopyable>
      ("ElementOp2", no_init)
      ;
    implicitly_convertible<Ref<SCElementOp2>, Ref<SavableState> >();

    class_<SCElementOp3,Ref<SCElementOp3>,bases<SavableState>(),boost::noncopyable>
      ("ElementOp3", no_init)
      ;
    implicitly_convertible<Ref<SCElementOp3>, Ref<SavableState> >();

    void (RefSCVector::*vecp)(std::ostream&) const = &RefSCVector::print;
    void (RefSCVector::*vecelemop1)(const Ref<SCElementOp>&) const
      = &RefSCVector::element_op;
    void (RefSCVector::*vecassd)(double) const
      = &RefSCVector::assign;

    class_<RefSCVector> ("SCVector",
                         init<const RefSCDimension &,
                         const Ref<SCMatrixKit>&>())
      .def("n",&RefSCVector::n)
      .def("dim",&RefSCVector::dim)
      .def("get_element",&RefSCVector::get_element)
      .def("set_element",&RefSCVector::set_element)
      .def("kit",&RefSCVector::kit)
      .def("clone",&RefSCVector::clone)
      .def("copy",&RefSCVector::copy)
      .def("element_op",vecelemop1)
      .def("assign",vecassd)
      .def("operator_add",&RefSCVector::operator+)
      .def("operator_mul",&RefSCVector::operator*)
      .def("scale",&RefSCVector::scale)
      .def("vprint",vecp);
    ;

    void (RefSCMatrix::*matp)(std::ostream&) const = &RefSCMatrix::print;
    void (RefSCMatrix::*matelemop1)(const Ref<SCElementOp>&) const
      = &RefSCMatrix::element_op;
    void (RefSCMatrix::*matassd)(double) const
      = &RefSCMatrix::assign;

    {
      RefSCMatrix (RefSCMatrix::*op_times1)(const RefSCMatrix&) const
        = &RefSCMatrix::operator*;
      RefSCMatrix (RefSCMatrix::*op_times2)(const RefSymmSCMatrix&) const
        = &RefSCMatrix::operator*;
      RefSCMatrix (RefSCMatrix::*op_times3)(const RefDiagSCMatrix&) const
        = &RefSCMatrix::operator*;
      RefSCVector (RefSCMatrix::*op_times4)(const RefSCVector&) const
        = &RefSCMatrix::operator*;
      RefSCMatrix (RefSCMatrix::*op_times5)(double) const
        = &RefSCMatrix::operator*;
      scope inner = class_<RefSCMatrix> ("SCMatrix",
                                         init<const RefSCDimension &,
                                         const RefSCDimension &,
                                         const Ref<SCMatrixKit>&>())
        .def("nrow",&RefSCMatrix::nrow)
        .def("ncol",&RefSCMatrix::ncol)
        .def("rowdim",&RefSCMatrix::rowdim)
        .def("coldim",&RefSCMatrix::coldim)
        .def("get_element",&RefSCMatrix::get_element)
        .def("set_element",&RefSCMatrix::set_element)
        .def("kit",&RefSCMatrix::kit)
        .def("clone",&RefSCMatrix::clone)
        .def("copy",&RefSCMatrix::copy)
        .def("element_op",matelemop1)
        .def("assign",matassd)
        .def("scale",&RefSCMatrix::scale)
        .def("t",&RefSCMatrix::t)
        .def("trace",&RefSCMatrix::trace)
        .def("operator_add",&RefSCMatrix::operator+)
        .def("operator_mul",op_times1)
        .def("operator_mul",op_times2)
        .def("operator_mul",op_times3)
        .def("operator_mul",op_times4)
        .def("operator_mul",op_times5)
        .def("mprint",matp)
        ;

      enum_<SCMatrix::Transform>("Transform")
        .value("NormalTransform",SCMatrix::NormalTransform)
        .value("TransposeTransform",SCMatrix::TransposeTransform)
        ;
    }

    {
      RefSCMatrix (RefSymmSCMatrix::*op_times1)(const RefSCMatrix&) const
        = &RefSymmSCMatrix::operator*;
      RefSCMatrix (RefSymmSCMatrix::*op_times2)(const RefSymmSCMatrix&) const
        = &RefSymmSCMatrix::operator*;
      RefSCMatrix (RefSymmSCMatrix::*op_times3)(const RefDiagSCMatrix&) const
        = &RefSymmSCMatrix::operator*;
      RefSCVector (RefSymmSCMatrix::*op_times4)(const RefSCVector&) const
        = &RefSymmSCMatrix::operator*;
      RefSymmSCMatrix (RefSymmSCMatrix::*op_times5)(double) const
        = &RefSymmSCMatrix::operator*;

      void (RefSymmSCMatrix::*smatp)(std::ostream&) const = &RefSymmSCMatrix::print;
      void (RefSymmSCMatrix::*smatelemop1)(const Ref<SCElementOp>&) const
        = &RefSymmSCMatrix::element_op;
      void (RefSymmSCMatrix::*smatassd)(double) const
        = &RefSymmSCMatrix::assign;
      void (RefSymmSCMatrix::*smatat1)(const RefSCMatrix&a,const RefSymmSCMatrix&b,
                                       SCMatrix::Transform) const
        = &RefSymmSCMatrix::accumulate_transform;
      void (RefSymmSCMatrix::*smatat2)(const RefSCMatrix&a,const RefDiagSCMatrix&b,
                                       SCMatrix::Transform) const
        = &RefSymmSCMatrix::accumulate_transform;
      void (RefSymmSCMatrix::*smatat3)(const RefSymmSCMatrix&a,
                                       const RefSymmSCMatrix&b) const
        = &RefSymmSCMatrix::accumulate_transform;

      class_<RefSymmSCMatrix> ("SymmSCMatrix",
                               init<const RefSCDimension&,
                               const Ref<SCMatrixKit>&>())
        .def("n",&RefSymmSCMatrix::n)
        .def("dim",&RefSymmSCMatrix::dim)
        .def("get_element",&RefSymmSCMatrix::get_element)
        .def("set_element",&RefSymmSCMatrix::set_element)
        .def("kit",&RefSymmSCMatrix::kit)
        .def("clone",&RefSymmSCMatrix::clone)
        .def("copy",&RefSymmSCMatrix::copy)
        .def("accumulate_transform",smatat1)
        .def("accumulate_transform",smatat2)
        .def("accumulate_transform",smatat3)
        .def("diagonalize",&RefSymmSCMatrix::diagonalize)
        .def("element_op",smatelemop1)
        .def("assign",smatassd)
        .def("scale",&RefSymmSCMatrix::scale)
        .def("scale_diagonal",&RefSymmSCMatrix::scale_diagonal)
        .def("operator_add",&RefSymmSCMatrix::operator+)
        .def("trace",&RefSymmSCMatrix::trace)
        .def("operator_mul",op_times1)
        .def("operator_mul",op_times2)
        .def("operator_mul",op_times3)
        .def("operator_mul",op_times4)
        .def("operator_mul",op_times5)
        .def("mprint",smatp)
        ;
    }

    {
      RefSCMatrix (RefDiagSCMatrix::*op_times1)(const RefSCMatrix&) const
        = &RefDiagSCMatrix::operator*;
      RefSCMatrix (RefDiagSCMatrix::*op_times2)(const RefSymmSCMatrix&) const
        = &RefDiagSCMatrix::operator*;
      RefDiagSCMatrix (RefDiagSCMatrix::*op_times3)(const RefDiagSCMatrix&) const
        = &RefDiagSCMatrix::operator*;
      RefDiagSCMatrix (RefDiagSCMatrix::*op_times4)(double) const
        = &RefDiagSCMatrix::operator*;

      void (RefDiagSCMatrix::*dmatp)(std::ostream&) const = &RefDiagSCMatrix::print;
      void (RefDiagSCMatrix::*dmatelemop1)(const Ref<SCElementOp>&) const
        = &RefDiagSCMatrix::element_op;
      void (RefDiagSCMatrix::*dmatassd)(double) const
        = &RefDiagSCMatrix::assign;

      class_<RefDiagSCMatrix> ("DiagSCMatrix",
                               init<const RefSCDimension &,
                               const Ref<SCMatrixKit>&>())
        .def("n",&RefDiagSCMatrix::n)
        .def("dim",&RefDiagSCMatrix::dim)
        .def("get_element",&RefDiagSCMatrix::get_element)
        .def("set_element",&RefDiagSCMatrix::set_element)
        .def("kit",&RefDiagSCMatrix::kit)
        .def("clone",&RefDiagSCMatrix::clone)
        .def("copy",&RefDiagSCMatrix::copy)
        .def("element_op",dmatelemop1)
        .def("assign",dmatassd)
        .def("scale",&RefDiagSCMatrix::scale)
        .def("trace",&RefDiagSCMatrix::trace)
        .def("operator_add",&RefSCMatrix::operator+)
        .def("operator_mul",op_times1)
        .def("operator_mul",op_times2)
        .def("operator_mul",op_times3)
        .def("operator_mul",op_times4)
        .def("mprint",matp)
        ;
    }
  }
}

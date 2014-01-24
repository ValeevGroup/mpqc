
#include <Python.h>
#include <extern/python/util.h>
#include <boost/python.hpp>
#include <boost/python/implicit.hpp>
#include <util/keyval/keyval.h>

using namespace boost::python;
using namespace sc;

namespace sc {
  void
  export_util_keyval()
  {
    class_<KeyValValue, boost::noncopyable >
      ("KeyValValue", no_init)
    ;

    class_<KeyValValueint, bases<KeyValValue>, boost::noncopyable >
      ("KeyValValueint", init<>())
      .def(init<int>())
    ;

    class_<KeyValValuedouble, bases<KeyValValue>, boost::noncopyable >
      ("KeyValValuedouble", init<>())
      .def(init<double>())
    ;

    class_<KeyValValuefloat, bases<KeyValValue>, boost::noncopyable >
      ("KeyValValuefloat", init<>())
      .def(init<float>())
    ;

    class_<KeyValValuechar, bases<KeyValValue>, boost::noncopyable >
      ("KeyValValuechar", init<>())
      .def(init<char>())
    ;

    class_<KeyValValuestring, bases<KeyValValue>, boost::noncopyable >
      ("KeyValValuestring", init<>())
      .def(init<std::string>())
    ;

    int (KeyVal::*v10)(const char*,const KeyValValue&) = &KeyVal::intvalue;
    int (KeyVal::*v11)(const char*,int,const KeyValValue&) = &KeyVal::intvalue;
    int (KeyVal::*v12)(const char*,int,int,const KeyValValue&) = &KeyVal::intvalue;
    double (KeyVal::*v20)(const char*,const KeyValValue&) = &KeyVal::doublevalue;
    double (KeyVal::*v21)(const char*,int,const KeyValValue&) = &KeyVal::doublevalue;
    double (KeyVal::*v22)(const char*,int,int,const KeyValValue&) = &KeyVal::doublevalue;
    float (KeyVal::*v30)(const char*,const KeyValValue&) = &KeyVal::floatvalue;
    float (KeyVal::*v31)(const char*,int,const KeyValValue&) = &KeyVal::floatvalue;
    float (KeyVal::*v32)(const char*,int,int,const KeyValValue&) = &KeyVal::floatvalue;
    char (KeyVal::*v40)(const char*,const KeyValValue&) = &KeyVal::charvalue;
    char (KeyVal::*v41)(const char*,int,const KeyValValue&) = &KeyVal::charvalue;
    char (KeyVal::*v42)(const char*,int,int,const KeyValValue&) = &KeyVal::charvalue;
    std::string (KeyVal::*v50)(const char*,const KeyValValue&) = &KeyVal::stringvalue;
    std::string (KeyVal::*v51)(const char*,int,const KeyValValue&) = &KeyVal::stringvalue;
    std::string (KeyVal::*v52)(const char*,int,int,const KeyValValue&) = &KeyVal::stringvalue;
    Ref<DescribedClass> (KeyVal::*v80)(const char*,const KeyValValue&) = &KeyVal::describedclassvalue;
    Ref<DescribedClass> (KeyVal::*v81)(const char*,int,const KeyValValue&) = &KeyVal::describedclassvalue;
    Ref<DescribedClass> (KeyVal::*v82)(const char*,int,int,const KeyValValue&) = &KeyVal::describedclassvalue;
    int (KeyVal::*v60)(const char*) = &KeyVal::exists;
    int (KeyVal::*v61)(const char*,int) = &KeyVal::exists;
    int (KeyVal::*v62)(const char*,int,int) = &KeyVal::exists;
    int (KeyVal::*v70)(const char*) = &KeyVal::count;
    int (KeyVal::*v71)(const char*,int) = &KeyVal::count;
    int (KeyVal::*v72)(const char*,int,int) = &KeyVal::count;

    class_<KeyVal, Ref<KeyVal>, boost::noncopyable >
      ("KeyVal", no_init)
      .def("intvalue",v10)
      .def("intvalue",v11)
      .def("intvalue",v12)
      .def("doublevalue",v20)
      .def("doublevalue",v21)
      .def("doublevalue",v22)
      .def("floatvalue",v30)
      .def("floatvalue",v31)
      .def("floatvalue",v32)
      .def("charvalue",v40)
      .def("charvalue",v41)
      .def("charvalue",v42)
      .def("stringvalue",v50)
      .def("stringvalue",v51)
      .def("stringvalue",v52)
      .def("exists",v60)
      .def("exists",v61)
      .def("exists",v62)
      .def("count",v70)
      .def("count",v71)
      .def("count",v72)
      .def("describedclassvalue",v80)
      .def("describedclassvalue",v81)
      .def("describedclassvalue",v82)
    ;

    void (AssignedKeyVal::*a1)(const char*,int) = &AssignedKeyVal::assign;
    void (AssignedKeyVal::*a2)(const char*,double) = &AssignedKeyVal::assign;
    void (AssignedKeyVal::*a3)(const char*,float) = &AssignedKeyVal::assign;
    void (AssignedKeyVal::*a4)(const char*,char) = &AssignedKeyVal::assign;
    void (AssignedKeyVal::*a5)(const char*,const char *) = &AssignedKeyVal::assign;
    void (AssignedKeyVal::*a6)(const char*,const Ref<DescribedClass>&) = &AssignedKeyVal::assign;

    class_<AssignedKeyVal, Ref<AssignedKeyVal>, bases<KeyVal>, boost::noncopyable >
      ("AssignedKeyVal", init<>())
      // Python overloading can unexpected results, so each member is
      // given a type-specific name.
      .def("assignboolean",&AssignedKeyVal::assignboolean)
      .def("assignint",a1)
      .def("assigndouble",a2)
      .def("assignfloat",a3)
      .def("assignchar",a4)
      .def("assignstring",a5)
      .def("assigndescribedclass",a6)
      .def("clear",&AssignedKeyVal::clear)
    ;
    implicitly_convertible<Ref<AssignedKeyVal>, Ref<KeyVal> >();

    class_<StringKeyVal, Ref<StringKeyVal>, bases<KeyVal>, boost::noncopyable >
      ("StringKeyVal", no_init)
    ;
    implicitly_convertible<Ref<StringKeyVal>, Ref<KeyVal> >();

    class_<ParsedKeyVal, Ref<ParsedKeyVal>, bases<StringKeyVal>, boost::noncopyable >
      ("ParsedKeyVal", init<const char *>())
    ;
    implicitly_convertible<Ref<ParsedKeyVal>, Ref<KeyVal> >();
  }
}

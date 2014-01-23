
#include <Python.h>
#include <numpy/arrayobject.h>
#include <extern/python/util.h>
#include <boost/python.hpp>
#include <chemistry/qc/basis/tbint.h>

using namespace boost::python;
using namespace sc;

namespace sc {
  void export_util_ref();
  void export_util_misc();
  void export_util_class();
  void export_util_keyval();
  void export_util_state();
  void export_math_scmat();
  void export_math_opt();
  void export_chemistry_molecule();
  void export_chemistry_qc_basis();
  void export_chemistry_qc_intv3();
  void export_chemistry_qc_libint2();
  void export_chemistry_qc_wfn();
  void export_chemistry_qc_scf();

  struct convert_arr1_to_numpy
  {
    static PyObject* convert(std::pair<const double *,unsigned long> const& x)
    {
      npy_intp dim = x.second;
      PyObject *r = PyArray_SimpleNewFromData(1,&dim,NPY_DOUBLE,(char*)x.first);
      return r;
    }
  };

  struct convert_arr2_to_numpy
  {
    static PyObject* convert(std::pair<const double *, std::array<unsigned long, 2> > const& x)
    {
      npy_intp dim[2];
      dim[0] = x.second[0];
      dim[1] = x.second[1];
      PyObject *r = PyArray_SimpleNewFromData(2,dim,NPY_DOUBLE,(char*)x.first);
      return r;
    }
  };

  struct convert_arr4_to_numpy
  {
    static PyObject* convert(std::pair<const double *, std::array<unsigned long, 4> > const& x)
    {
      npy_intp dim[4];
      dim[0] = x.second[0];
      dim[1] = x.second[1];
      dim[2] = x.second[2];
      dim[3] = x.second[3];
      PyObject *r = PyArray_SimpleNewFromData(4,dim,NPY_DOUBLE,(char*)x.first);
      return r;
    }
  };

  class TestArray1 {
    int n;
    double *d;
  public:
    ~TestArray1() { delete[] d; }
    TestArray1() { d = new double[1]; d[0] = 1.2345; n = 1; }
    TestArray1(int dim): n(dim) {
      d = new double[dim];
      for (int i=0; i<dim; i++) d[i] = i + 0.9876;
    }
    std::pair<const double*,unsigned long> result() { return std::make_pair<const double*,unsigned long>(d,n); }
  };

  class TestArray2 {
    int n1,n2;
    double *d;
  public:
    ~TestArray2() { delete[] d; }
    TestArray2() { d = new double[1]; d[0] = 1.2345; n1 = 1; n2 = 1; }
    TestArray2(int dim1, int dim2): n1(dim1),n2(dim2) {
      d = new double[dim1*dim2];
      for (int i=0; i<dim1; i++)
        for (int j=0; j<dim2; j++)
          d[i*dim2+j] = i + 0.001*j;
    }
    std::pair<const double*,std::array<unsigned long, 2> > result() {
      std::pair<const double*,std::array<unsigned long, 2> > r;
      r.first = d;
      r.second[0] = n1;
      r.second[1] = n2;
      return r;
    }
  };

  class TestArray4 {
    int n1,n2,n3,n4;
    double *d;
  public:
    ~TestArray4() { delete[] d; }
    TestArray4() { d = new double[1]; d[0] = 1.2345; n1 = 1; n2 = 1; n3 = 1; n4 = 1; }
    TestArray4(int dim1, int dim2, int dim3, int dim4): n1(dim1),n2(dim2),n3(dim3),n4(dim4) {
      d = new double[dim1*dim2*dim3*dim4];
      for (int i=0; i<dim1; i++)
        for (int j=0; j<dim2; j++)
          for (int k=0; k<dim3; k++)
            for (int l=0; l<dim4; l++)
              d[((i*dim2+j)*dim3+k)*dim4+l] = i + 0.001*j + 0.000001*k + 0.0000000001*l;
    }
    std::pair<const double*,std::array<unsigned long, 4> > result() {
      std::pair<const double*,std::array<unsigned long, 4> > r;
      r.first = d;
      r.second[0] = n1;
      r.second[1] = n2;
      r.second[2] = n3;
      r.second[3] = n4;
      return r;
    }
  };

}

BOOST_PYTHON_MODULE(libmpqc_py)
{
  import_array();

  to_python_converter<std::pair<const double*,unsigned long>,
                      convert_arr1_to_numpy>();

  to_python_converter<std::pair<const double*,std::array<unsigned long, 2> >,
                      convert_arr2_to_numpy>();

  to_python_converter<std::pair<const double*,std::array<unsigned long, 4> >,
                      convert_arr4_to_numpy>();

  class_<TestArray1, boost::noncopyable>("TestArray1",init<>())
    .def(init<int>())
    .def("result",&TestArray1::result)
    ;

  class_<TestArray2, boost::noncopyable>("TestArray2",init<>())
    .def(init<int,int>())
    .def("result",&TestArray2::result)
    ;

  class_<TestArray4, boost::noncopyable>("TestArray4",init<>())
    .def(init<int,int,int,int>())
    .def("result",&TestArray4::result)
    ;

  ///////////////////////////////////////////////////////////
  // Utility classes

  export_util_ref();
  export_util_class();
  export_util_keyval();
  export_util_state();
  export_util_misc();

  ///////////////////////////////////////////////////////////
  // Math classes

  export_math_scmat();
  export_math_opt();

  ///////////////////////////////////////////////////////////
  // Chemistry classes

  export_chemistry_molecule();
  export_chemistry_qc_basis();
  export_chemistry_qc_intv3();
  export_chemistry_qc_libint2();
  export_chemistry_qc_wfn();
  export_chemistry_qc_scf();

}


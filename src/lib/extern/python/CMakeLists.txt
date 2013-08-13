# mpqc.Python library

if (HAVE_PYTHON AND HAVE_PYTHON_NUMPY)

  include_directories(${PYTHON_INCLUDE_DIRS})
  include_directories(${PY_NUMPY}/core/include)

  add_library(mpqc.Python
    chemistry_molecule.cc
    chemistry_qc_basis.cc
    chemistry_qc_intv3.cc
    chemistry_qc_libint2.cc
    chemistry_qc_scf.cc
    chemistry_qc_wfn.cc
    math_opt.cc
    math_scmat.cc
    mpqc.cc
    util_class.cc
    util_keyval.cc
    util_misc.cc
    util_ref.cc
    util_state.cc
  )

  target_link_libraries(mpqc.Python chemistry ${PYTHON_LIBRARIES} ${Boost_LIBRARIES})
endif()

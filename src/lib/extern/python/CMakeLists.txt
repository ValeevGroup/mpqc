# mpqc.Python library

if (HAVE_PYTHON AND HAVE_PYTHON_NUMPY)

  MESSAGE("====== To use the MPQC python packages, must install first =====")

  include_directories(${PYTHON_INCLUDE_DIRS})
  include_directories(${PY_NUMPY}/core/include)

  add_library(mpqc_py MODULE
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

  target_link_libraries(mpqc_py chemistry ${PYTHON_LIBRARIES} ${Boost_LIBRARIES})
  
  INSTALL( TARGETS mpqc_py DESTINATION lib/pympqc )
  INSTALL( FILES __init__.py DESTINATION lib/pympqc )
  INSTALL( FILES ase_mpqc.py DESTINATION lib/pympqc )
  
  # announce to the user how to use pympqc
  GET_FILENAME_COMPONENT(PY_NUMPY_PATH ${PY_NUMPY}/.. ABSOLUTE)
  INSTALL( CODE "MESSAGE(\"To use pympqc package, add these directories to PYTHONPATH: ${DESTDIR}${CMAKE_INSTALL_PREFIX}/lib and ${PY_NUMPY_PATH}\")" )

  configure_file(testbasic.py.in testbasic.py @ONLY IMMEDIATE)
  configure_file(testhf.py.in testhf.py @ONLY IMMEDIATE)
  configure_file(testhf2.py.in testhf2.py @ONLY IMMEDIATE)
  configure_file(testintegral.py.in testintegral.py @ONLY IMMEDIATE)  
  configure_file(testase.py.in testase.py @ONLY IMMEDIATE)  
  
endif()

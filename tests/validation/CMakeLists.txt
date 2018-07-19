if (${PYTHONINTERP_FOUND})
  enable_testing(true)
  add_custom_target(check_validation COMMAND ${CMAKE_CTEST_COMMAND} -V)

  set(REF_INPUT_PATH "${CMAKE_CURRENT_SOURCE_DIR}/reference/inputs")

  # Find all the input files under ${REF_INPUT_PATH}
  file(GLOB inFiles RELATIVE "${REF_INPUT_PATH}"
      "${REF_INPUT_PATH}/*.json" "${REF_INPUT_PATH}/*.xml" "${REF_INPUT_PATH}/*.info")

  foreach(infileName ${inFiles})

    string(REGEX REPLACE ".(json|xml|info)\$" "" baseName "${infileName}")

    # filter out tests based on the available MPQC features
    # TODO also filter based on the default TA policy
    # skip CP-CCSD tests until fixed
    if ((${baseName} MATCHES "(cp-ccsd)"))
        continue()
    endif()
    if ((${baseName} MATCHES "(ccsd|cc3)") AND NOT "${MPQC_FEATURES}" MATCHES "lcao_cc")
        continue()
    endif()
    if ((${baseName} MATCHES "eom" OR ${baseName} MATCHES "cis") AND NOT "${MPQC_FEATURES}" MATCHES "lcao_ci")
        continue()
    endif()
    if ((${baseName} MATCHES "f12") AND NOT "${MPQC_FEATURES}" MATCHES "lcao_f12")
        continue()
    endif()
    if ((${baseName} MATCHES "mp2" OR ${baseName} MATCHES "gf2") AND NOT "${MPQC_FEATURES}" MATCHES "lcao_mbpt")
        continue()
    endif()
    if ((${baseName} MATCHES "hf") AND NOT "${MPQC_FEATURES}" MATCHES "lcao_scf")
        continue()
    endif()

    add_test(NAME validation/${baseName}-np1
            COMMAND ${CMAKE_COMMAND}
            -DtestName=${baseName}
            -DprintOutput=${MPQC_VALIDATION_TEST_PRINT}
            -DsrcDir=${CMAKE_CURRENT_SOURCE_DIR}
            -DpythonExec=${PYTHON_EXECUTABLE}
            -DmpiExec=${MPIEXEC}
            -DmpiNPFlags=${MPIEXEC_NUMPROC_FLAG}
            -DmpiNProc=1
            -DmpiPre=${MPIEXEC_PREFLAGS}
            -DmpiPost=${MPIEXEC_POSTFLAGS}
            -P ${CMAKE_CURRENT_SOURCE_DIR}/test.cmake)
    set_tests_properties(validation/${baseName}-np1
            PROPERTIES DEPENDS ../../src/bin/mpqc/mpqc)

    add_test(NAME validation/${baseName}-np2
             COMMAND ${CMAKE_COMMAND}
              -DtestName=${baseName}
              -DprintOutput=${MPQC_VALIDATION_TEST_PRINT}
              -DsrcDir=${CMAKE_CURRENT_SOURCE_DIR}
              -DpythonExec=${PYTHON_EXECUTABLE}
              -DmpiExec=${MPIEXEC}
              -DmpiNPFlags=${MPIEXEC_NUMPROC_FLAG}
              -DmpiNProc=2
              -DmpiPre=${MPIEXEC_PREFLAGS}
              -DmpiPost=${MPIEXEC_POSTFLAGS}
              -P ${CMAKE_CURRENT_SOURCE_DIR}/test.cmake)
    set_tests_properties(validation/${baseName}-np2
                         PROPERTIES DEPENDS ../../src/bin/mpqc/mpqc)
    
  endforeach(infileName)
else (${PYTHONINTERP_FOUND})
  add_custom_target(check_validation COMMENT "Python interpreter is not found, will skip validation tests")
endif (${PYTHONINTERP_FOUND})
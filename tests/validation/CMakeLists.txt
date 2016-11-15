if (${PYTHONINTERP_FOUND})
  enable_testing(true)
  add_custom_target(check_validation COMMAND ${CMAKE_CTEST_COMMAND} -V)

  set(REF_INPUT_PATH "${CMAKE_CURRENT_SOURCE_DIR}/reference/inputs")

  # Find all the input files under ${REF_INPUT_PATH}
  file(GLOB inFiles RELATIVE "${REF_INPUT_PATH}"
      "${REF_INPUT_PATH}/*.json")

  foreach(infileName ${inFiles})

    string(REGEX REPLACE ".json\$" "" baseName "${infileName}")
                         
    add_test(NAME validation/${baseName}
             COMMAND ${CMAKE_COMMAND}
              -DtestName=${baseName}
              -DsrcDir=${CMAKE_CURRENT_SOURCE_DIR}
              -DpythonExec=${PYTHON_EXECUTABLE}
              -DmpiExec=${MPIEXEC}
              -DmpiNPFlags=${MPIEXEC_NUMPROC_FLAG}
              -DmpiNProc=2
              -DmpiPre=${MPIEXEC_PREFLAGS}
              -DmpiPost=${MPIEXEC_POSTFLAGS}
              -P ${CMAKE_CURRENT_SOURCE_DIR}/test.cmake)
    set_tests_properties(validation/${baseName}
                         PROPERTIES DEPENDS ../../src/bin/mpqc/mpqc)
    
  endforeach(infileName)
else (${PYTHONINTERP_FOUND})
  add_custom_target(check_validation COMMENT "Python interpreter is not found, will skip validation tests")
endif (${PYTHONINTERP_FOUND})
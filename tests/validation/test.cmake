
macro(runtest)

set(OUTPUT_FILE_NAME "${CMAKE_BINARY_DIR}/${testName}.out")

set(CHECK_CMD "${pythonExec}")
set(CHECK_ARGS "${srcDir}/check.py"
"${OUTPUT_FILE_NAME}"
"${srcDir}/reference/outputs/${testName}.out")

set(ENV{MAD_NUM_THREADS} 2)
#if (NOT EXISTS "${OUTPUT_FILE_NAME}")
if (1)

  set(MPQC_CMD "${CMAKE_BINARY_DIR}/../../src/bin/mpqc/mpqc")

  # Retrieve the list of DescribedClass classes registered with this current executable
  execute_process(COMMAND
          ${MPQC_CMD} "-k"
          RESULT_VARIABLE MPQC_DC_RESULT
          OUTPUT_VARIABLE MPQC_DC_OUTPUT)

  string(REPLACE "\n" " " MPQC_DC_OUTPUT "${MPQC_DC_OUTPUT}")

  # filter out tests based on the registered classes
  # parse the wfn type from the input file, make sure it has a match in the registered class list
  if (${MPQC_DC_RESULT} EQUAL 0)
      file(READ "${srcDir}/reference/inputs/${testName}.json" infileContents)
      string(REGEX REPLACE ".*\"wfn\"[\r\n\t ]*:[\r\n\t ]*{[\r\n\t ]*\"type\"[\r\n\t ]*:[\r\n\t ]*\"\([-a-zA-Z0-9 _]+\)\".*"
              "\\1" wfnType "${infileContents}")
      if (NOT "${MPQC_DC_OUTPUT}" MATCHES "${wfnType}")
          message(STATUS "skipped test ${testName}")
          return()
      endif()
  endif()

  set(MPQC_ARGS "-D" "-p" "${srcDir}/reference/inputs"
  "-i" "${srcDir}/reference/inputs/${testName}.json")
  execute_process(COMMAND
                  ${mpiExec}
                  ${mpiNPFlags}
                  ${mpiNProc}
                  ${mpiPre}
                  $ENV{MPQC_PRE_CMD}
                  ${MPQC_CMD} ${MPQC_ARGS}
                  ${mpiPost}
                  OUTPUT_FILE "${OUTPUT_FILE_NAME}"
                  RESULT_VARIABLE MPQC_RESULT)

  if(MPQC_RESULT)
    if(printOutput)
      message(STATUS "\nOUTPUT of " ${testName})
      execute_process(COMMAND
              cat
              ${OUTPUT_FILE_NAME}
              RESULT_VARIABLE
              CAT_RESULT
              )
    endif()
    message(FATAL_ERROR "Error running ${MPQC_CMD}")
  endif()
endif()

execute_process(COMMAND
                ${CHECK_CMD} ${CHECK_ARGS}
                RESULT_VARIABLE CHECK_RESULT)
if(CHECK_RESULT)
    if(printOutput)
      message(STATUS "\nOUTPUT of " ${testName})
      execute_process(COMMAND
              cat
              ${OUTPUT_FILE_NAME}
              RESULT_VARIABLE
              CAT_RESULT
              )
    endif()

    message(FATAL_ERROR "Error running ${CHECK_CMD}")
endif(CHECK_RESULT)

endmacro(runtest)

runtest()

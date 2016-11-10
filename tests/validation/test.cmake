
macro(runtest srcDir baseName)

set(OUTPUT_FILE_NAME "${CMAKE_BINARY_DIR}/${baseName}.out")

set(CHECK_CMD "${pythonExec}")
set(CHECK_ARGS "${srcDir}/check.py"
"${OUTPUT_FILE_NAME}"
"${srcDir}/reference/outputs/${baseName}.out")

if (NOT EXISTS "${OUTPUT_FILE_NAME}")

  set(MPQC_CMD "${CMAKE_BINARY_DIR}/../../src/bin/mpqc/mpqc")
  set(MPQC_ARGS "-p ${srcDir}/reference/inputs"
      "${srcDir}/reference/inputs/${baseName}.json")

  execute_process(COMMAND
                  ${MPQC_CMD} ${MPQC_ARGS}
                  OUTPUT_FILE "${OUTPUT_FILE_NAME}"
                  RESULT_VARIABLE MPQC_RESULT)

  if(MPQC_RESULT)
    message(STATUS "output: ${OUTPUT_VARIABLE}")
    message(FATAL_ERROR "Error running ${MPQC_CMD}")
  endif()
endif()

execute_process(COMMAND
                ${CHECK_CMD} ${CHECK_ARGS}
                RESULT_VARIABLE CHECK_RESULT)
if(CHECK_RESULT)
    message(FATAL_ERROR "Error running ${CHECK_CMD}")
endif(CHECK_RESULT)

endmacro(runtest)

runtest(${srcDir} ${testName} ${pythonExec})

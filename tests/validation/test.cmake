
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
  set(MPQC_ARGS "-p" "${srcDir}/reference/inputs"
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

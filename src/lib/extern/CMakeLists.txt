
add_subdirectory(qc)
add_subdirectory(moinfo)

if (HAVE_PYTHON AND HAVE_PYTHON_NUMPY)
  add_subdirectory(python)
endif()

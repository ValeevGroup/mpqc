add_subdirectory(qc)

if (HAVE_PYTHON AND HAVE_PYTHON_NUMPY)
  add_subdirectory(python)
endif()

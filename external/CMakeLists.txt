add_custom_target(External)

include(external/Lapack)

if (MPQC_NEW_FEATURES)
  include(external/Eigen)
  include(external/Boost)
  include(external/HDF5)
  include(external/ARMCI)
  include(external/Libint)
endif()



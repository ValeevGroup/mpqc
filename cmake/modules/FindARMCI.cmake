# Try to find ARMCI headers and libraries.
#
# Usage of this module as follows:
#
#     find_package(ARMCI)
#
# Variables used by this module, they can change the default behaviour and need
# to be set before calling find_package:
#
#  ARMCI_PREFIX         Set this variable to the root installation of
#                      libpapi if the module has problems finding the
#                      proper installation path.
#
# Variables defined by this module:
#
#  ARMCI_FOUND              System has ARMCI libraries and headers
#  ARMCI_LIBRARIES          The ARMCI library
#  ARMCI_INCLUDE_DIRS       The location of ARMCI headers

find_path(ARMCI_PREFIX
    NAMES include/armci.h
)

find_library(ARMCI_LIBRARIES
    # Pick the static library first for easier run-time linking.
    NAMES libarmci.a armci
    HINTS ${ARMCI_PREFIX}/lib ${HILTIDEPS}/lib
)

find_path(ARMCI_INCLUDE_DIRS
    NAMES armci.h
    HINTS ${ARMCI_PREFIX}/include ${HILTIDEPS}/include
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(ARMCI DEFAULT_MSG
    ARMCI_LIBRARIES
    ARMCI_INCLUDE_DIRS
)

mark_as_advanced(
    ARMCI_PREFIX_DIRS
    ARMCI_LIBRARIES
    ARMCI_INCLUDE_DIRS
)

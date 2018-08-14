#
#  This file is a part of TiledArray.
#  Copyright (C) 2013  Virginia Tech
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#  Justus Calvin
#  Department of Chemistry, Virginia Tech
#
#  CMakeLists.txt
#  Jul 19, 2013
#
#  Copied from TiledArray on June 21 2015, by Drew Lewis

set (CMAKE_DOCINSTALL_PATH ${CMAKE_INSTALL_PREFIX}/share/doc/mpqc-${MPQC_VERSION})

find_package(Doxygen)
if(DOXYGEN_FOUND)
  configure_file(Doxyfile.in Doxyfile @ONLY IMMEDIATE)
  add_custom_target (doc
    COMMAND ${DOXYGEN_EXECUTABLE} ${PROJECT_BINARY_DIR}/doc/Doxyfile
    SOURCES ${PROJECT_BINARY_DIR}/doc/Doxyfile)
#  INSTALL(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html DESTINATION ${CMAKE_DOCINSTALL_PATH}
#          COMPONENT mpqc-docs)
endif()

INSTALL(FILES ${PROJECT_SOURCE_DIR}/doc/devsamp/mp2/mp2.cpp
              ${PROJECT_SOURCE_DIR}/doc/devsamp/mp2/mp2.json
              ${PROJECT_BINARY_DIR}/doc/devsamp/mp2/CMakeLists.txt
        DESTINATION ${CMAKE_DOCINSTALL_PATH}/examples/mp2 COMPONENT mpqc-docs)
INSTALL(FILES ${PROJECT_SOURCE_DIR}/doc/devsamp/aoints/aoints.cpp
        ${PROJECT_SOURCE_DIR}/doc/devsamp/aoints/aoints.json
        ${PROJECT_BINARY_DIR}/doc/devsamp/aoints/CMakeLists.txt
        DESTINATION ${CMAKE_DOCINSTALL_PATH}/examples/aoints COMPONENT mpqc-docs)

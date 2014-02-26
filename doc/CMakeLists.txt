set (CMAKE_DOCINSTALL_PATH ${CMAKE_INSTALL_PREFIX}/share/doc/mpqc-${MPQC_VERSION})
INSTALL(DIRECTORY ./devsamp/hf DESTINATION ${CMAKE_DOCINSTALL_PATH}/examples)
INSTALL(DIRECTORY ./devsamp/mp2 DESTINATION ${CMAKE_DOCINSTALL_PATH}/examples)
INSTALL(DIRECTORY ./devsamp/aoints DESTINATION ${CMAKE_DOCINSTALL_PATH}/examples)

#-- Add an Option to toggle the generation of the API documentation

find_package(Doxygen)
if(DOXYGEN_FOUND)
  configure_file(Doxyfile.in Doxyfile @ONLY IMMEDIATE)
  
  add_custom_command( OUTPUT ${PROJECT_BINARY_DIR}/doc/html/index.html
                      COMMAND ${DOXYGEN_EXECUTABLE} ${PROJECT_BINARY_DIR}/doc/Doxyfile
                      MAIN_DEPENDENCY ${PROJECT_BINARY_DIR}/doc/Doxyfile ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in
                      DEPENDS
                        ${CMAKE_CURRENT_SOURCE_DIR}/DoxygenLayout.xml
                        ${CMAKE_CURRENT_SOURCE_DIR}/compile.dox
                        ${CMAKE_CURRENT_SOURCE_DIR}/develop.dox
                        ${CMAKE_CURRENT_SOURCE_DIR}/footer.html
                        ${CMAKE_CURRENT_SOURCE_DIR}/overview.dox
                        ${CMAKE_CURRENT_SOURCE_DIR}/sf_footer.html
                        ${CMAKE_CURRENT_SOURCE_DIR}/srctree.dox
                        ${CMAKE_CURRENT_SOURCE_DIR}/modules.dox
                        ${CMAKE_CURRENT_SOURCE_DIR}/top.dox
                        ${CMAKE_SOURCE_DIR}/bin/mpqcrun.dox
                        ${CMAKE_SOURCE_DIR}/bin/mpqcrunning.dox
                        ${CMAKE_SOURCE_DIR}/src/bin/molrender/molrender.dox
                        ${CMAKE_SOURCE_DIR}/src/bin/mpqc/mpqc.dox
                        ${CMAKE_SOURCE_DIR}/src/bin/mpqc/mpqccomp.dox
                        ${CMAKE_SOURCE_DIR}/src/bin/mpqc/mpqcinp.dox
                        ${CMAKE_SOURCE_DIR}/src/bin/mpqc/mpqclic.dox
                        ${CMAKE_SOURCE_DIR}/src/bin/mpqc/mpqcoo.dox
                        ${CMAKE_SOURCE_DIR}/src/bin/mpqc/mpqcover.dox
                        ${CMAKE_SOURCE_DIR}/src/bin/mpqc/mpqcsimp.dox
                        ${CMAKE_SOURCE_DIR}/src/bin/mpqc/mpqcval.dox
                        ${CMAKE_SOURCE_DIR}/src/bin/mpqc/mpqcwar.dox
                        ${CMAKE_SOURCE_DIR}/src/bin/scls/scls.dox
                        ${CMAKE_SOURCE_DIR}/src/bin/scpr/scpr.dox
                        ${CMAKE_SOURCE_DIR}/src/lib/chemistry/qc/psi/mpqcpsi.dox
                        ${CMAKE_SOURCE_DIR}/src/lib/math/scmat/scmat.dox
                        ${CMAKE_SOURCE_DIR}/src/lib/math/symmetry/symmetry.dox
                        ${CMAKE_SOURCE_DIR}/src/lib/util/class/class.dox
                        ${CMAKE_SOURCE_DIR}/src/lib/util/group/group.dox
                        ${CMAKE_SOURCE_DIR}/src/lib/util/keyval/keyval.dox
                        ${CMAKE_SOURCE_DIR}/src/lib/util/ref/ref.dox
                        ${CMAKE_SOURCE_DIR}/src/lib/util/state/state.dox
                      IMPLICIT_DEPENDS
                        CXX ${CMAKE_SOURCE_DIR}/src/bin/mpqc/mpqc.cc
                        CXX ${CMAKE_SOURCE_DIR}/src/bin/pt2r12/pt2r12.cc
                      COMMENT "Generating MPQC API documentation with Doxygen")
  
  INSTALL( DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html DESTINATION ${CMAKE_DOCINSTALL_PATH})

  add_custom_target (html
    DEPENDS ${PROJECT_BINARY_DIR}/doc/html/index.html)

endif()

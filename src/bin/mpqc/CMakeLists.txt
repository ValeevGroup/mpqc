add_library(
  MPQCinit
  mpqcinit.cpp
)

target_link_libraries(
  MPQCinit
  MPQCutil_options
  MPQCutil_keyval
  MPQCutil_misc
)

add_executable(mpqc EXCLUDE_FROM_ALL mpqc.cpp)
# add the libraries you need and force their linkage in mpqc.cpp
# developers will want to trim this list down to the minimum needed
target_link_libraries(mpqc
                      MPQCf12 MPQCcc MPQCmbpt MPQCscf
                      MPQCutil_options MPQCutil_misc
                      MPQCinit)

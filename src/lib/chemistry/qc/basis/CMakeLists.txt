# basis library

add_library(basis OBJECT
  aotoso.cc
  cart.cc
  cartiter.cc
  dercent.cc
  distsh.cc
  distshpair.cc
  extent.cc
  files.cc
  fjt.cc
  gaussbas.cc
  gaussbaval.cc
  gaussshell.cc
  gaussshval.cc
  gpetite.cc
  integral.cc
  intparams.cc
  inttraits.cc
  lselect.cc
  obint.cc
  operator.cc
  orthog.cc
  petite.cc
  shellrot.cc
  sobasis.cc
  sointegral.cc
  split.cc
  symmint.cc
  tbint.cc
  transform.cc
  uncontract.cc
  union.cc
)

# tests

#set_property(SOURCE btest.cc PROPERTY COMPILE_DEFINITIONS
#    -DSRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
#add_executable(btest EXCLUDE_FROM_ALL btest.cc $<TARGET_OBJECTS:basis>)
#add_test(btest btest)

#set_property(SOURCE tbint_batch_test.cc PROPERTY COMPILE_DEFINITIONS
#    -DSRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
#add_executable(tbint_batch_test EXCLUDE_FROM_ALL tbint_batch_test.cc $<TARGET_OBJECTS:basis>)
#add_test(tbint_batch_test tbint_batch_test)


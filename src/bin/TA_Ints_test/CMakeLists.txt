# Test of Drew's New Stuff


set(CMAKE_BUILD_TYPE DEBUG)
set(CMAKE_CXX_FLAGS_DEBUG "-g -O2 -std=c++11 -stdlib=libc++")
add_executable(IntegralEnginePool_Test
   EXCLUDE_FROM_ALL IntegralEnginePool_Test.cc )

target_link_libraries(IntegralEnginePool_Test chemistry)
target_link_libraries(IntegralEnginePool_Test math)
target_link_libraries(IntegralEnginePool_Test util)
target_link_libraries(IntegralEnginePool_Test chemistry)
target_link_libraries(IntegralEnginePool_Test math)
target_link_libraries(IntegralEnginePool_Test util)
target_link_libraries(IntegralEnginePool_Test ${LINK_LIBRARIES})

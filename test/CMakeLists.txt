add_library(gtest_mpi_parallel_main SHARED)
target_sources(gtest_mpi_parallel_main PRIVATE gtest_mpi_parallel_main.cpp)
target_link_libraries(gtest_mpi_parallel_main PUBLIC GTest::gtest_mpi)
add_library(GTest::gtest_mpi_main ALIAS gtest_mpi_parallel_main)

# implementation of the test infrastructure
include(GoogleTest)

add_executable(FSClippingTest FSClipping/FSClippingTestGeo.cpp 
                              FSClipping/FSClippingTestFaceExtractor.cpp 
                              FSClipping/FSClippingTestMatches.cpp 
                              FSClipping/FSClippingTestCell2NodeBuilder.cpp 
                              FSClipping/FSClippingTestPolyFaceBuilder.cpp
                              FSClipping/FSClippingTestHelper.cpp)
target_link_libraries(FSClippingTest PUBLIC GTest::gtest_mpi_main fsclipping::fsclipping)
set_property(TARGET FSClippingTest PROPERTY CROSSCOMPILING_EMULATOR '${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1')
gtest_discover_tests(FSClippingTest DISCOVERY_TIMEOUT 20 PROPERTIES LABELS "CppUnitTests")

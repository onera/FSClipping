add_library(gtest_mpi_parallel_main SHARED)
target_sources(gtest_mpi_parallel_main PRIVATE gtest_mpi_parallel_main.cpp)
target_link_libraries(gtest_mpi_parallel_main PUBLIC GTest::gtest_mpi)
add_library(GTest::gtest_mpi_main ALIAS gtest_mpi_parallel_main)

add_compile_definitions(FS_TEST_DATA_DIR="${CMAKE_SOURCE_DIR}/test/Mesh")
include(GoogleTest)

# Sequential tests — 1 proc
add_executable(FSClippingTest FSClipping/FSClippingTestGeo.cpp
                              FSClipping/FSClippingTestBoundaryExtractor.cpp
                              FSClipping/FSClippingTestMatches.cpp
                              FSClipping/FSClippingTestCell2NodeBuilder.cpp
                              FSClipping/FSClippingTestPolyFaceBuilder.cpp
                              FSClipping/FSClippingTestInterface.cpp)
target_link_libraries(FSClippingTest PUBLIC GTest::gtest_mpi_main fsclipping::fsclipping)
set_property(TARGET FSClippingTest PROPERTY CROSSCOMPILING_EMULATOR '${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1')
gtest_discover_tests(FSClippingTest DISCOVERY_TIMEOUT 20 PROPERTIES LABELS "CppUnitTests")

# Parallel tests — 3 procs (clipper proc receives both face sets and computes matches)
add_executable(FSClippingParallelTest FSClipping/Parallel/FSClippingTestInterfaceParallel3.cpp)
target_link_libraries(FSClippingParallelTest PUBLIC GTest::gtest_mpi_main fsclipping::fsclipping)
target_include_directories(FSClippingParallelTest PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/FSClipping)
set_property(TARGET FSClippingParallelTest PROPERTY CROSSCOMPILING_EMULATOR '${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 3')
gtest_discover_tests(FSClippingParallelTest DISCOVERY_TIMEOUT 20 PROPERTIES LABELS "CppUnitTests")

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
                              FSClipping/FSClippingTestGlobalNumbering.cpp
                              FSClipping/FSClippingTestInterface.cpp)
target_link_libraries(FSClippingTest PUBLIC GTest::gtest_mpi_main fsclipping::fsclipping)
set_property(TARGET FSClippingTest PROPERTY CROSSCOMPILING_EMULATOR '${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1')
gtest_discover_tests(FSClippingTest DISCOVERY_TIMEOUT 20 PROPERTIES LABELS "CppUnitTests")

# Parallel tests — one executable per MPI proc count, each launched at exactly
# that proc count. gtest_discover_tests reads the proc count from the target's
# CROSSCOMPILING_EMULATOR property, which can hold only one value, so each proc
# count must be a separate executable. Every test is compiled into the single
# executable matching its proc layout, so no test skips itself and the CTest
# history has no SKIP noise:
#   - np3 : send/receive match tests (3 procs: 1 subject + 1 clipped + 1 matcher)
#   - np5 : symmetric layout 1 clipper + 2 (mesh 1) + 2 (mesh 2)
#           cube surface/volume + 2-cylinder volume
#   - np6 : asymmetric layout 1 clipper + 2 (mesh 1) + 3 (mesh 2)
#           cube surface/volume + 2-cylinder surface/volume
#
# add_parallel_test(<nproc> <source>...) builds one such executable and registers
# it with CTest under the np<nproc>. prefix.
function(add_parallel_test nproc)
  set(target FSClippingParallelTest_np${nproc})
  add_executable(${target} ${ARGN})
  target_link_libraries(${target} PUBLIC GTest::gtest_mpi_main fsclipping::fsclipping)
  target_include_directories(${target} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/FSClipping
                                              ${CMAKE_CURRENT_SOURCE_DIR}/FSClipping/Parallel)
  set_property(TARGET ${target}
               PROPERTY CROSSCOMPILING_EMULATOR '${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${nproc}')
  gtest_discover_tests(${target}
    TEST_PREFIX "np${nproc}."
    DISCOVERY_TIMEOUT 40
    PROPERTIES LABELS "CppUnitTests")
endfunction()

add_parallel_test(3 FSClipping/Parallel/FSClippingTestMatchesSendReceive.cpp)
add_parallel_test(5 FSClipping/Parallel/FSClippingTestInterfaceSym.cpp
                    FSClipping/Parallel/FSClippingTestMatchRouting.cpp)
add_parallel_test(6 FSClipping/Parallel/FSClippingTestInterfaceAsym.cpp)

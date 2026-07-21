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

# Parallel tests — the same test binary is run in CTest at several proc counts.
# Each TEST body skips itself (SKIP_UNLESS_*_PROCS) when the current world size
# is not the one its proc layout needs, so the whole suite can be launched at 3,
# 5 and 7 procs and each test runs only where it is meant to (others SKIPPED):
#   - 3 procs : send/receive match tests, cube surface/volume interface
#   - 5 procs : symmetric 2-cylinder volume reconstruction
#   - 7 procs : asymmetric 4+2+1 cylinder surface/volume reconstruction
#
# gtest_discover_tests reads the proc count from the target's CROSSCOMPILING_EMULATOR
# property, which can hold only one value — so we build one executable per proc
# count. The sources are compiled once into an OBJECT library and shared, so the
# extra executables cost only a link. TEST_PREFIX (np3./np5./np7.) keeps the CTest
# names unique across the three registrations.
add_library(FSClippingParallelTestObj OBJECT
            FSClipping/Parallel/FSClippingTestInterfaceParallel.cpp
            FSClipping/Parallel/FSClippingTestMatchesSendReceive.cpp)
target_link_libraries(FSClippingParallelTestObj PUBLIC GTest::gtest_mpi_main fsclipping::fsclipping)
target_include_directories(FSClippingParallelTestObj PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/FSClipping
                                                            ${CMAKE_CURRENT_SOURCE_DIR}/FSClipping/Parallel)

foreach(nproc 3 5 7)
  add_executable(FSClippingParallelTest${nproc} $<TARGET_OBJECTS:FSClippingParallelTestObj>)
  target_link_libraries(FSClippingParallelTest${nproc} PUBLIC FSClippingParallelTestObj)
  set_property(TARGET FSClippingParallelTest${nproc}
               PROPERTY CROSSCOMPILING_EMULATOR '${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${nproc}')
  gtest_discover_tests(FSClippingParallelTest${nproc}
    TEST_PREFIX "np${nproc}."
    DISCOVERY_TIMEOUT 40
    PROPERTIES LABELS "CppUnitTests")
endforeach()

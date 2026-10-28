add_library(gtest_mpi_parallel_main SHARED)
target_sources(gtest_mpi_parallel_main PRIVATE gtest_mpi_parallel_main.cpp)
target_link_libraries(gtest_mpi_parallel_main PUBLIC GTest::gtest_mpi)
add_library(GTest::gtest_mpi_main ALIAS gtest_mpi_parallel_main)

add_compile_definitions(FS_TEST_DATA_DIR="${CMAKE_SOURCE_DIR}/test/Mesh")
include(GoogleTest)

# Detect available logical processors at configure time so that
# parallel tests can be auto-disabled when the hardware cannot run them.
cmake_host_system_information(RESULT AVAILABLE_CORES QUERY NUMBER_OF_LOGICAL_CORES)

# Sequential tests — 1 proc
add_executable(FSClippingTest FSClipping/FSClippingTestGeo.cpp
                              FSClipping/FSClippingTestBoundaryExtractor.cpp
                              FSClipping/FSClippingTestMatches.cpp
                              FSClipping/FSClippingTestGlobalNumbering.cpp
                              FSClipping/FSClippingTestInterface.cpp
                              FSClipping/FSClippingTestCellClipped.cpp
                              FSClipping/FSClippingTestCellSource.cpp
                              FSClipping/FSClippingTestNodeRegistry.cpp
                              FSClipping/FSClippingTestPolyCell.cpp
                              FSClipping/FSClippingTestTopologySink.cpp
                              FSClipping/FSClippingTestCellAttributes.cpp)
target_link_libraries(FSClippingTest PUBLIC GTest::gtest_mpi_main fsclipping::fsclipping)
set_property(TARGET FSClippingTest PROPERTY CROSSCOMPILING_EMULATOR '${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1')
gtest_discover_tests(FSClippingTest DISCOVERY_TIMEOUT 20 PROPERTIES LABELS "CppUnitTests")

# Death tests end in MPI_Abort, which kills the whole job under mpiexec: run them as an MPI singleton
add_executable(FSClippingDeathTest FSClipping/FSClippingTestCellClipped.cpp
                                   FSClipping/FSClippingTestNodeRegistry.cpp)
target_link_libraries(FSClippingDeathTest PUBLIC GTest::gtest_mpi_main fsclipping::fsclipping)
gtest_discover_tests(FSClippingDeathTest DISCOVERY_TIMEOUT 20 TEST_FILTER "*Death*" PROPERTIES LABELS "CppUnitTests")

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
  if(AVAILABLE_CORES AND AVAILABLE_CORES LESS nproc)
    # Doc on DIASBLED https://cmake.org/cmake/help/latest/prop_test/DISABLED.html
    message(STATUS "Disabling ${target}: requires ${nproc} cores but only ${AVAILABLE_CORES} available")
    set(_discover_props "LABELS;CppUnitTests;DISABLED;TRUE")
  else()
    set(_discover_props "LABELS;CppUnitTests")
  endif()
  # DISCOVERY_MODE PRE_TEST seems to solve issue related to "mpirun -np 1 cmd"
  gtest_discover_tests(${target}
    TEST_PREFIX "np${nproc}."
    DISCOVERY_TIMEOUT 40
    PROPERTIES ${_discover_props})
endfunction()

add_parallel_test(3 FSClipping/Parallel/FSClippingTestMatchesSendReceive.cpp)
add_parallel_test(4 FSClipping/Parallel/FSClippingTestBoundaryExtractionPar.cpp)
add_parallel_test(5 FSClipping/Parallel/FSClippingTestInterfaceSym.cpp
                    FSClipping/Parallel/FSClippingTestMatchRouting.cpp)
add_parallel_test(9 FSClipping/Parallel/FSClippingTestPoly2DAttributesPar.cpp)
add_parallel_test(6 FSClipping/Parallel/FSClippingTestInterfaceAsym.cpp)

# Python examples, run against the module of the build tree (build/py) rather than the installed one
set(example_dir ${CMAKE_SOURCE_DIR}/example/meshOp)
set(example_env "PYTHONPATH=${CMAKE_BINARY_DIR}/py:$ENV{PYTHONPATH}")

add_test(NAME example.checkForDataManagerOp
         COMMAND ${Python_EXECUTABLE} ${example_dir}/checkForDataManagerOp.py)
set_tests_properties(example.checkForDataManagerOp PROPERTIES
                     ENVIRONMENT "${example_env}"
                     PASS_REGULAR_EXPRESSION "FSClippedMesh is registered : *True"
                     LABELS "PythonExamples")

add_test(NAME example.howToComputeClippedMesh
         COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 6 ${Python_EXECUTABLE}
                 ${example_dir}/howToComputeClippedMesh.py)
set_tests_properties(example.howToComputeClippedMesh PROPERTIES
                     ENVIRONMENT "${example_env}"
                     FAIL_REGULAR_EXPRESSION "\\*\\*\\* ERROR \\*\\*\\*"
                     LABELS "PythonExamples")

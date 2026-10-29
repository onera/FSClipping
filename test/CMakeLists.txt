add_library(gtest_mpi_parallel_main SHARED)
target_sources(gtest_mpi_parallel_main PRIVATE gtest_mpi_parallel_main.cpp)
target_link_libraries(gtest_mpi_parallel_main PUBLIC GTest::gtest_mpi)
fs_enable_warnings(gtest_mpi_parallel_main)
add_library(GTest::gtest_mpi_main ALIAS gtest_mpi_parallel_main)

include(GoogleTest)

# Detect available logical processors at configure time so that
# parallel tests can be auto-disabled when the hardware cannot run them.
cmake_host_system_information(RESULT AVAILABLE_CORES QUERY NUMBER_OF_LOGICAL_CORES)
set(test_data_dir "${CMAKE_SOURCE_DIR}/test/Mesh")

# Tests must load the build-tree library, even when an installed one is on the loader path
if(APPLE)
  set(loader_path_var DYLD_LIBRARY_PATH)
else()
  set(loader_path_var LD_LIBRARY_PATH)
endif()
set(build_lib_env "${loader_path_var}=path_list_prepend:${PROJECT_BINARY_DIR}/src")

# Sequential tests, 1 proc
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
target_compile_definitions(FSClippingTest PRIVATE FS_TEST_DATA_DIR="${test_data_dir}")
fs_enable_warnings(FSClippingTest)
set_property(TARGET FSClippingTest PROPERTY CROSSCOMPILING_EMULATOR '${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1')
gtest_discover_tests(FSClippingTest DISCOVERY_TIMEOUT 20 PROPERTIES LABELS "CppUnitTests" ENVIRONMENT_MODIFICATION "${build_lib_env}")

# Parallel tests: CROSSCOMPILING_EMULATOR holds a single proc count, so one executable per count
function(add_parallel_test nproc)
  set(target FSClippingParallelTest_np${nproc})
  add_executable(${target} ${ARGN})
  target_link_libraries(${target} PUBLIC GTest::gtest_mpi_main fsclipping::fsclipping)
  target_compile_definitions(${target} PRIVATE FS_TEST_DATA_DIR="${test_data_dir}")
  fs_enable_warnings(${target})
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
    PROPERTIES LABELS "CppUnitTests" ENVIRONMENT_MODIFICATION "${build_lib_env}")
endfunction()

add_parallel_test(3 FSClipping/Parallel/FSClippingTestMatchesSendReceive.cpp)
add_parallel_test(4 FSClipping/Parallel/FSClippingTestBoundaryExtractionPar.cpp)
add_parallel_test(5 FSClipping/Parallel/FSClippingTestInterfaceSym.cpp
                    FSClipping/Parallel/FSClippingTestMatchRouting.cpp)
add_parallel_test(9 FSClipping/Parallel/FSClippingTestPoly2DAttributesPar.cpp)
add_parallel_test(6 FSClipping/Parallel/FSClippingTestInterfaceAsym.cpp)

# Python examples, run against the build-tree module rather than the installed one
set(example_dir ${CMAKE_SOURCE_DIR}/example/meshOp)
set(example_env "PYTHONPATH=path_list_prepend:${PROJECT_BINARY_DIR}/py" "${build_lib_env}")

add_test(NAME example.checkForDataManagerOp
         COMMAND ${Python_EXECUTABLE} ${example_dir}/checkForDataManagerOp.py)
set_tests_properties(example.checkForDataManagerOp PROPERTIES
                     ENVIRONMENT_MODIFICATION "${example_env}"
                     PASS_REGULAR_EXPRESSION "FSClippedMesh is registered : *True"
                     LABELS "PythonExamples")

set(nproc 6)
add_test(NAME example.howToComputeClippedMesh
         COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${nproc} ${Python_EXECUTABLE}
                 ${example_dir}/howToComputeClippedMesh.py)
set_tests_properties(example.howToComputeClippedMesh PROPERTIES
                     ENVIRONMENT_MODIFICATION "${example_env}"
                     FAIL_REGULAR_EXPRESSION "\\*\\*\\* ERROR \\*\\*\\*"
                     LABELS "PythonExamples")
if(AVAILABLE_CORES AND AVAILABLE_CORES LESS nproc)
  # Doc on DIASBLED https://cmake.org/cmake/help/latest/prop_test/DISABLED.html
  message(STATUS "Disabling test example.howToComputeClippedMesh: requires ${nproc} cores but only ${AVAILABLE_CORES} available")
  set_tests_properties(example.howToComputeClippedMesh PROPERTIES DISABLED TRUE)
endif()

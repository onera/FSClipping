cmake_minimum_required(VERSION 3.22)  # ENVIRONMENT_MODIFICATION for tests

# Default to the MPI wrapper so FSDM and MPI headers are found without extra flags
if(NOT CMAKE_CXX_COMPILER)
    find_program(MPICXX mpicxx)
    if(MPICXX)
        set(CMAKE_CXX_COMPILER "${MPICXX}" CACHE FILEPATH "MPI C++ wrapper")
    endif()
endif()

project(FSClipping
        VERSION 0.1.0
        LANGUAGES CXX
        DESCRIPTION "Sutherland‑Hodgman clipping for two FSMesh")

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Optimized with debug symbols unless the user chooses otherwise
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Build type" FORCE)
endif()

find_package(PkgConfig REQUIRED)
set(FSDM_MIN_VERSION 25.08)
pkg_check_modules(FSDM REQUIRED IMPORTED_TARGET "fsdm >= ${FSDM_MIN_VERSION}")
message(STATUS "Found FSDM v${FSDM_VERSION}")

find_package(MPI REQUIRED COMPONENTS CXX)
message(STATUS "MPI implementation: ${MPI_CXX_LIBRARIES}")

# SWIG module and Python examples
find_package(Python 3.10 REQUIRED COMPONENTS Interpreter Development)

if (NOT APPLE)
    set(CMAKE_SKIP_BUILD_RPATH OFF)
    # CMAKE_INSTALL_RPATH_USE_LINK_PATH misses the FSDM dirs coming from $LIBRARY_PATH
    set(CMAKE_INSTALL_RPATH "${FSDM_LIBRARY_DIRS}")
endif ()

include(GNUInstallDirs)

option(FSCLIPPING_WARNINGS_AS_ERRORS "Treat compiler warnings as errors" OFF)
# Warnings on our own targets only, not on gtest nor the SWIG-generated code
function(fs_enable_warnings tgt)
    target_compile_options(${tgt} PRIVATE -Wall -Wextra -Wpedantic
                           $<$<BOOL:${FSCLIPPING_WARNINGS_AS_ERRORS}>:-Werror>)
endfunction()

add_subdirectory(src)
add_subdirectory(py)

include(CTest)
if (BUILD_TESTING)
    message(STATUS "Enabling building of C++ tests")
    add_subdirectory(thirdparty/gtest-gmock-1.14.0)
    add_subdirectory(test)
endif ()

# pkg-config file for non-CMake consumers
set(LINK_FLAGS "-lfsclipping")
configure_file(tools/fsclipping.pc.in
               ${CMAKE_CURRENT_BINARY_DIR}/fsclipping.pc @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/fsclipping.pc
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)

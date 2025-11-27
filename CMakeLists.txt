# --------------------------------------------------------------
# Root CMakeLists.txt – project entry point
# --------------------------------------------------------------
cmake_minimum_required(VERSION 3.20)

project(FSClipping
        VERSION 0.1.0
        LANGUAGES CXX
        DESCRIPTION "Sutherland‑Hodgman clipping for two FSMesh")

# --------------------------------------------------------------
# Options 
# --------------------------------------------------------------
option(FSCLIPPING_BUILD_TESTS "Compile unit tests" OFF)

# --------------------------------------------------------------
# Compiler settings
# --------------------------------------------------------------
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)          # no compiler‑specific extensions
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)   # for clang‑tidy, VSCode, etc.

# --------------------------------------------------------------
# External dependencies
# --------------------------------------------------------------
# 1) FSDM ≥ 25.08 (pkg‑config)
find_package(PkgConfig REQUIRED)
pkg_check_modules(FSDM REQUIRED IMPORTED_TARGET "fsdm >= 25.08")
message(STATUS "Found FSDM v${FSDM_VERSION}")

# 2) MPI (OpenMPI 4.0+)
find_package(MPI REQUIRED COMPONENTS CXX)
message(STATUS "MPI implementation: ${MPI_CXX_LIBRARIES}")

# 3) Python ≥ 3.10 (interpreter + development headers) – only needed for includes
find_package(Python 3.10 REQUIRED COMPONENTS Interpreter Development)

# --------------------------------------------------------------
# RPATH handling (non‑Apple platforms)
# --------------------------------------------------------------
if (NOT APPLE)
    set(CMAKE_SKIP_BUILD_RPATH OFF)
    # CMAKE_INSTALL_RPATH_USE_LINK_PATH does not work for library paths contained in $LIBRARY_PATH,
    # so we have to set the RPATH manually.
    set(CMAKE_INSTALL_RPATH "${FSDM_LIBRARY_DIRS}")
endif ()

# --------------------------------------------------------------
# Warning handling (disabled for now – uncomment if you need it)
# --------------------------------------------------------------
#function(enable_warnings tgt)
#    if (MSVC)
#        target_compile_options(${tgt} PRIVATE /W4 /WX)
#    else()
#        target_compile_options(${tgt} PRIVATE -Wall -Wextra -Wpedantic -Werror)
#    endif()
#endfunction()

# --------------------------------------------------------------
# Standard installation directories (GNUInstallDirs defines the variables)
# --------------------------------------------------------------
include(GNUInstallDirs)   # defines CMAKE_INSTALL_LIBDIR, CMAKE_INSTALL_INCLUDEDIR, ...

# --------------------------------------------------------------
# Add the src sub‑directory (contains the simplified CMakeLists.txt)
# --------------------------------------------------------------
add_subdirectory(src)

# --------------------------------------------------------------
# Add the test sub‑directory (contains the simplified CMakeLists.txt)
# --------------------------------------------------------------
include(CTest)
if (BUILD_TESTING)
    add_subdirectory(thirdparty/gtest-gmock-1.14.0)
    add_subdirectory(test)
endif ()

# --------------------------------------------------------------
# Install a pkg‑config file so that other projects can link against us
# --------------------------------------------------------------
configure_file(tools/fsclipping.pc.in
               ${CMAKE_CURRENT_BINARY_DIR}/fsclipping.pc @ONLY)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/fsclipping.pc
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)
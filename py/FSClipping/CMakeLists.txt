set(package_name FSClipping)

find_package(SWIG 4.0 REQUIRED COMPONENTS python)

cmake_policy(SET CMP0078 NEW)
cmake_policy(SET CMP0086 NEW)

include(UseSWIG)

set(swig_sources_list ${CMAKE_SOURCE_DIR}/include/FSClipping/_FSClippedMesh.i)

foreach(source ${swig_sources_list})
  set_property(SOURCE ${source} PROPERTY CPLUSPLUS on)
  set_property(SOURCE ${source} PROPERTY SWIG_MODULE_NAME ${package_name})
  set_property(SOURCE ${source} PROPERTY USE_TARGET_INCLUDE_DIRECTORIES true)
endforeach()

set(swig_module_target ${package_name}_swig)

set(CMAKE_SWIG_FLAGS -fcompact -fvirtual)

if(${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.20.0")
  set(SWIG_USE_SWIG_DEPENDENCIES true)
endif()

if(${SWIG_VERSION} VERSION_LESS "4.1.0")
  list(APPEND CMAKE_SWIG_FLAGS -py3)
endif()

if(${SWIG_VERSION} VERSION_GREATER_EQUAL "4.0.1")
  # generate Python docstrings for SWIG-wrapped functions while disabling some
  # warnings
  list(APPEND CMAKE_SWIG_FLAGS -doxygen -w560)
endif()

swig_add_library(
  ${swig_module_target}
  TYPE MODULE
  LANGUAGE python
  SOURCES ${swig_sources_list})

set_target_properties(
  ${swig_module_target}
  PROPERTIES
    OUTPUT_NAME
      ${package_name}
)

target_link_libraries(
  ${swig_module_target}
  PRIVATE
    fsclipping::fsclipping
    Python::Module
)

add_library(fsclipping::python ALIAS ${swig_module_target})

# additional Python files
foreach(python_file IN ITEMS __init__.py)
  configure_file(${python_file} ${CMAKE_CURRENT_BINARY_DIR}/${python_file}
                 COPYONLY)
endforeach()

set_target_properties(
  ${swig_module_target} PROPERTIES INSTALL_RPATH
                                   "$ORIGIN/../../../;${CMAKE_INSTALL_RPATH}")

cmake_path(
  APPEND
  python_install_dir
  "lib"
  "python${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}"
  "site-packages"
  "${package_name}"
)

install(
  TARGETS
    ${swig_module_target}
  LIBRARY
    DESTINATION
      ${python_install_dir}
    COMPONENT
      FSClipping_Python
)

install(
  FILES
    ${CMAKE_CURRENT_BINARY_DIR}/__init__.py
    ${CMAKE_CURRENT_BINARY_DIR}/FSClipping.py
  DESTINATION
    ${python_install_dir}
  COMPONENT
    FSClipping_Python
)

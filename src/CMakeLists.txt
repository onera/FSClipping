# --------------------------------------------------------------
# src/CMakeLists.txt – builds the shared library
# --------------------------------------------------------------

# --------------------------------------------------------------
# Shared library definition
# --------------------------------------------------------------
add_library(fsclipping SHARED   FSClippingUtil.cpp
                                FSClippingFace.cpp
                                FSFaceSeparator.cpp
                                FSFaceMatcher.cpp
                                FSCell2NodeBuilder.cpp
                                FSPolyFaceBuilder.cpp
                                FSBoundaryFaceProvider.cpp
                                FSTopologyAssembler.cpp
                                FSMeshReconstruction.cpp
                                FSClippingInterface.cpp
                                FSFaceExchange.cpp
                                FSClippingInterfacePar.cpp
                                FSClippedMeshParams.cpp
                                FSClippedMesh.cpp)

# Create an alias so that other targets can refer to it as fsclipping::fsclipping
add_library(fsclipping::fsclipping ALIAS fsclipping)

# --------------------------------------------------------------
# Public include directories (project headers + FSDM + Python)
# --------------------------------------------------------------
target_include_directories(fsclipping
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>   # for the installed package
        ${FSDM_INCLUDE_DIRS}
        ${Python_INCLUDE_DIRS}
)

# --------------------------------------------------------------
# Link against external libraries
# --------------------------------------------------------------
target_link_libraries(fsclipping
    PUBLIC
        PkgConfig::FSDM          # libfsdm.so / .dll
        MPI::MPI_CXX
)

# --------------------------------------------------------------
# Public headers (to be installed)
# --------------------------------------------------------------
set(HEADERS
    ${CMAKE_CURRENT_SOURCE_DIR}/../include/FSClipping/FSBoundaryFace.h
    ${CMAKE_CURRENT_SOURCE_DIR}/../include/FSClipping/FSFaceMatcher.h
    ${CMAKE_CURRENT_SOURCE_DIR}/../include/FSClipping/FSCellAdress.h
    ${CMAKE_CURRENT_SOURCE_DIR}/../include/FSClipping/FSClippingUtil.h
    ${CMAKE_CURRENT_SOURCE_DIR}/../include/FSClipping/FSFace.h
    ${CMAKE_CURRENT_SOURCE_DIR}/../include/FSClipping/FSFaceSeparator.h
    ${CMAKE_CURRENT_SOURCE_DIR}/../include/FSClipping/FSCell2NodeBuilder.h
    ${CMAKE_CURRENT_SOURCE_DIR}/../include/FSClipping/FSPolyFaceBuilder.h
    ${CMAKE_CURRENT_SOURCE_DIR}/../include/FSClipping/FSBoundaryFaceProvider.h
    ${CMAKE_CURRENT_SOURCE_DIR}/../include/FSClipping/FSTopologyAssembler.h
    ${CMAKE_CURRENT_SOURCE_DIR}/../include/FSClipping/FSMeshReconstruction.h
    ${CMAKE_CURRENT_SOURCE_DIR}/../include/FSClipping/FSClippingInterface.h
    ${CMAKE_CURRENT_SOURCE_DIR}/../include/FSClipping/FSFaceExchange.h
    ${CMAKE_CURRENT_SOURCE_DIR}/../include/FSClipping/FSClippingInterfacePar.h
    ${CMAKE_CURRENT_SOURCE_DIR}/../include/FSClipping/FSClippedMeshParams.h
    ${CMAKE_CURRENT_SOURCE_DIR}/../include/FSClipping/FSClippedMeshParams.i
    ${CMAKE_CURRENT_SOURCE_DIR}/../include/FSClipping/FSClippedMesh.h
    ${CMAKE_CURRENT_SOURCE_DIR}/../include/FSClipping/_FSClippedMesh.i
)

# Expose the list of public headers to the install command
set_target_properties(fsclipping PROPERTIES
    PUBLIC_HEADER "${HEADERS}"
)

# --------------------------------------------------------------
# Installation of the library
# --------------------------------------------------------------
install(TARGETS fsclipping
        EXPORT  fsclippingTargets               # (optional – for a CMake package config)
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}   # .so / .dll
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}   # .a (if you also create a static lib)
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}   # executables linked to the lib (unlikely here)
        PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/FSClipping
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

# --------------------------------------------------------------
# (Optional) Export a CMake package configuration so other projects can
# find this library with `find_package(fsclipping CONFIG)`
# --------------------------------------------------------------
install(EXPORT fsclippingTargets
        FILE fsclippingTargets.cmake
        NAMESPACE fsclipping::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/fsclipping)
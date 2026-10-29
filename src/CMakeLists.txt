add_library(fsclipping SHARED   FSClippingUtil.cpp
                                FSClippingFace.cpp
                                FSFaceSeparator.cpp
                                FSFaceMatcher.cpp
                                FSBoundaryFaceProvider.cpp
                                FSMeshReconstruction.cpp
                                FSClippingEngine.cpp
                                FSFaceExchange.cpp
                                FSCellSource.cpp
                                FSCellClipped.cpp
                                FSNodeRegistry.cpp
                                FSPolyCell.cpp
                                FSTopologySink.cpp
                                FSClippedMeshParams.cpp
                                FSClippedMesh.cpp)
add_library(fsclipping::fsclipping ALIAS fsclipping)
fs_enable_warnings(fsclipping)

target_include_directories(fsclipping
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

target_link_libraries(fsclipping
    PUBLIC
        PkgConfig::FSDM
        MPI::MPI_CXX
)

# Installed headers, including the SWIG interfaces
set(HEADERS
    ${CMAKE_CURRENT_SOURCE_DIR}/../include/FSClipping/FSBoundaryFace.h
    ${CMAKE_CURRENT_SOURCE_DIR}/../include/FSClipping/FSFaceMatcher.h
    ${CMAKE_CURRENT_SOURCE_DIR}/../include/FSClipping/FSClippingUtil.h
    ${CMAKE_CURRENT_SOURCE_DIR}/../include/FSClipping/FSClippingFace.h
    ${CMAKE_CURRENT_SOURCE_DIR}/../include/FSClipping/FSFace.h
    ${CMAKE_CURRENT_SOURCE_DIR}/../include/FSClipping/FSFaceSeparator.h
    ${CMAKE_CURRENT_SOURCE_DIR}/../include/FSClipping/FSBoundaryFaceProvider.h
    ${CMAKE_CURRENT_SOURCE_DIR}/../include/FSClipping/FSMeshReconstruction.h
    ${CMAKE_CURRENT_SOURCE_DIR}/../include/FSClipping/FSNodeKey.h
    ${CMAKE_CURRENT_SOURCE_DIR}/../include/FSClipping/FSClippingEngine.h
    ${CMAKE_CURRENT_SOURCE_DIR}/../include/FSClipping/FSFaceExchange.h
    ${CMAKE_CURRENT_SOURCE_DIR}/../include/FSClipping/FSCellSource.h
    ${CMAKE_CURRENT_SOURCE_DIR}/../include/FSClipping/FSCellClipped.h
    ${CMAKE_CURRENT_SOURCE_DIR}/../include/FSClipping/FSNodeRegistry.h
    ${CMAKE_CURRENT_SOURCE_DIR}/../include/FSClipping/FSPolyCell.h
    ${CMAKE_CURRENT_SOURCE_DIR}/../include/FSClipping/FSTopologySink.h
    ${CMAKE_CURRENT_SOURCE_DIR}/../include/FSClipping/FSTopologyData.h
    ${CMAKE_CURRENT_SOURCE_DIR}/../include/FSClipping/FSClippedMeshParams.h
    ${CMAKE_CURRENT_SOURCE_DIR}/../include/FSClipping/FSClippedMeshParams.i
    ${CMAKE_CURRENT_SOURCE_DIR}/../include/FSClipping/FSClippedMesh.h
    ${CMAKE_CURRENT_SOURCE_DIR}/../include/FSClipping/_FSClippedMesh.i
)
set_target_properties(fsclipping PROPERTIES PUBLIC_HEADER "${HEADERS}")

install(TARGETS fsclipping
        EXPORT  fsclippingTargets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/FSClipping
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

# CMake package for find_package(fsclipping CONFIG)
set(config_install_dir ${CMAKE_INSTALL_LIBDIR}/cmake/fsclipping)

install(EXPORT fsclippingTargets
        FILE fsclippingTargets.cmake
        NAMESPACE fsclipping::
        DESTINATION ${config_install_dir})

include(CMakePackageConfigHelpers)
configure_package_config_file(${PROJECT_SOURCE_DIR}/tools/fsclippingConfig.cmake.in
                              ${CMAKE_CURRENT_BINARY_DIR}/fsclippingConfig.cmake
                              INSTALL_DESTINATION ${config_install_dir})
# 0.x: any minor bump may break the API
write_basic_package_version_file(${CMAKE_CURRENT_BINARY_DIR}/fsclippingConfigVersion.cmake
                                 VERSION ${PROJECT_VERSION}
                                 COMPATIBILITY SameMinorVersion)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/fsclippingConfig.cmake
              ${CMAKE_CURRENT_BINARY_DIR}/fsclippingConfigVersion.cmake
        DESTINATION ${config_install_dir})

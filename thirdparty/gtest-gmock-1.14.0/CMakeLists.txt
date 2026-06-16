add_library(gtest_mpi SHARED)
target_sources(gtest_mpi PRIVATE gmock-gtest-all.cc)
target_include_directories(gtest_mpi SYSTEM PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(gtest_mpi MPI::MPI_CXX)

add_library(GTest::gtest_mpi ALIAS gtest_mpi)

// test/TestUtils.hpp

#include <cstdlib>
#include <string>

#ifndef FS_TEST_DATA_DIR
#error "FS_TEST_DATA_DIR not defined"
#endif

const std::string testDir = FS_TEST_DATA_DIR;


inline std::string GetTestDataDir()
{
  const char* env = std::getenv("FS_TEST_DATA_DIR");
  if(env)
    return std::string(env);

  return std::string(FS_TEST_DATA_DIR);
}

inline std::string MeshPath(const std::string& file)
{
  return GetTestDataDir() + "/" + file;
}
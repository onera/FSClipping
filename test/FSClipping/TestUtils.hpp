// test/TestUtils.hpp

#include <cstdlib>
#include <string>
#include <atomic>
#include <ctime>
#include <unistd.h>
#include <filesystem>

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

// -----------------------------------------------------------------------------
// TempDir: RAII temporary directory
// Creates a unique directory on construction and removes it (recursively) on
// destruction.  Safe across MPI ranks: each process gets its own directory
// thanks to a combination of PID, timestamp and an atomic counter.
// -----------------------------------------------------------------------------
class TempDir
{
public:
  TempDir()
  {
    static std::atomic<int> counter{0};
    std::filesystem::path temp = std::filesystem::temp_directory_path();
    std::string name = "fsc_test_" + std::to_string(static_cast<long long>(::getpid())) + "_" +
                       std::to_string(std::time(nullptr)) + "_" +
                       std::to_string(counter.fetch_add(1, std::memory_order_relaxed));
    path_ = temp / name;
    std::filesystem::create_directories(path_);
  }

  ~TempDir()
  {
    std::error_code ec;
    std::filesystem::remove_all(path_, ec);
  }

  // Non-copyable, non-movable to preserve uniqueness
  TempDir(const TempDir&) = delete;
  TempDir& operator=(const TempDir&) = delete;
  TempDir(TempDir&&) = delete;
  TempDir& operator=(TempDir&&) = delete;

  const std::string& path() const { return path_; }
  operator const std::string&() const { return path_; }

private:
  std::string path_;
};

inline std::string TempMeshPath(const std::string& file, const TempDir& tmpDir)
{
  return tmpDir.path() + "/" + file;
}
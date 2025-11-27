#include <cstdlib>

#include <mpi.h>

#include "gtest/gtest.h"

class RootOnlyPrinter : public ::testing::EmptyTestEventListener {
 public:
  explicit RootOnlyPrinter(bool is_root_process) :
      _is_root_process{is_root_process},
      _printer{::testing::UnitTest::GetInstance()->listeners().Release(
          ::testing::UnitTest::GetInstance()->listeners().default_result_printer())},
      _xml_printer{::testing::UnitTest::GetInstance()->listeners().Release(
          ::testing::UnitTest::GetInstance()->listeners().default_xml_generator())} {}

  ~RootOnlyPrinter() override {
    delete _printer;
    delete _xml_printer;
  }

  RootOnlyPrinter(RootOnlyPrinter&) = delete;
  RootOnlyPrinter(RootOnlyPrinter&&) = delete;
  RootOnlyPrinter& operator=(RootOnlyPrinter&&) = delete;
  RootOnlyPrinter& operator=(const RootOnlyPrinter&) = delete;

  void OnTestCaseStart(const ::testing::TestCase& test_case) override {
    if (_is_root_process) {
      _printer->OnTestCaseStart(test_case);
    }
  }

  void OnTestStart(const ::testing::TestInfo& test_info) override {
    if (_is_root_process) {
      _printer->OnTestStart(test_info);
    }
  }

  void OnTestPartResult(const ::testing::TestPartResult& test_part_result) override {
    if (_is_root_process) {
      _printer->OnTestPartResult(test_part_result);
    }
  }
  void OnTestEnd(const ::testing::TestInfo& test_info) override {
    if (_is_root_process) {
      _printer->OnTestEnd(test_info);
    }
  }

  void OnTestIterationEnd(const ::testing::UnitTest& unit_test, int iteration) override {
    if (_is_root_process) {
      _printer->OnTestIterationEnd(unit_test, iteration);
      if (_xml_printer != nullptr) {
        _xml_printer->OnTestIterationEnd(unit_test, iteration);
      }
    }
  }

 private:
  ::testing::TestEventListener* _printer;
  ::testing::TestEventListener* _xml_printer;
  bool _is_root_process;
};

int main(int argc, char** argv) {
  auto required = MPI_THREAD_SERIALIZED;
  auto provided = -1;
  if (MPI_Init_thread(&argc, &argv, required, &provided) != MPI_SUCCESS || provided < required) {
    return EXIT_FAILURE;
  }

  int my_rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  const bool is_root_process = my_rank == 0;

  int status{0};

  // gtest wants to print on the master process no matter what.
  // in case of us just listing the tests: only pass the parameters on rank 0 (only print there) and exit
  if (std::count(argv, argv + static_cast<std::size_t>(argc), std::string{"--gtest_list_tests"}) > 0) {
    if (is_root_process) {
      ::testing::InitGoogleTest(&argc, argv);
      status = RUN_ALL_TESTS();
    }
  } else {
    // normal
    ::testing::InitGoogleTest(&argc, argv);

    // RootOnlyPrinter takes control of the default and xml printer
    ::testing::TestEventListeners& listeners = ::testing::UnitTest::GetInstance()->listeners();
    listeners.Append(new RootOnlyPrinter(is_root_process));  // NOLINT(cppcoreguidelines-owning-memory)
    status = RUN_ALL_TESTS();
  }

  std::ignore = MPI_Finalize();
  return status;
}

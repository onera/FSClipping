#include "FSClipping/FSBoundaryFaceProvider.h"
#include "FSClipping/FSFaceExchange.h"
#include "FSClipping/FSFaceMatcher.h"
#include "FSClipping/FSTopologyAssembler.h"
#include "TestUtilsParallel.hpp"
#include "gtest/gtest.h"

_FS_BEGIN_NAMESPACE

// Level 1 — Three-process test with a dedicated clipper proc.
//
// Setup (mirrors the real user pattern):
//   globalClac = FSClac(MPI_COMM_WORLD)       3 procs
//   meshID     = globalClac.GetProcID()        0, 1, or 2
//   clac       = FSClac(); globalClac.DivideIntoGroups(meshID, clac)
//
//   Proc 0 : mesh1 owner → extracts boundary faces → FSFaceExchange::Send to proc 2
//   Proc 1 : mesh2 owner → extracts boundary faces → FSFaceExchange::Send to proc 2
//   Proc 2 : clipper     → FSFaceExchange::Receive both → FSFaceMatcher → matches
//
// No cross-mesh data ever lives on proc 0 or proc 1.
TEST(FSClippingTestInterfaceParallel3, ComputeMatchesOnClipperProc)
{
  FSClac globalClac(MPI_COMM_WORLD);
  ASSERT_EQ(globalClac.GetNProcs(), 3) << "This test requires exactly 3 MPI processes";

  FS_intT meshID = globalClac.GetProcID();
  FSClac clac;
  globalClac.DivideIntoGroups(meshID, clac);

  const FS_intT clipperRank = 2;
  const FS_floatT tol = 1e-8;

  if(meshID == 0) {
    auto mesh = LoadMeshWithClac(clac, MeshPath("input/cube_hexa_coarse.grid"));
    FSMeshFaceExtractor extractor;
    auto extraction = FSBoundaryFaceProvider::Extract(mesh, extractor, 2, tol);
    FSFaceExchange::Send(globalClac, clipperRank, extraction.faces);
    std::cout << "[Proc 0] sent " << extraction.faces.size() << " faces (mesh1) to proc 2\n";

  } else if(meshID == 1) {
    auto mesh = LoadMeshWithClac(clac, MeshPath("input/cube_hexa_fine.grid"));
    FSMeshFaceExtractor extractor;
    auto extraction = FSBoundaryFaceProvider::Extract(mesh, extractor, 1, tol);
    FSFaceExchange::Send(globalClac, clipperRank, extraction.faces);
    std::cout << "[Proc 1] sent " << extraction.faces.size() << " faces (mesh2) to proc 2\n";

  } else {
    auto subjectFaces = FSFaceExchange::Receive(globalClac, 0);
    auto clippedFaces = FSFaceExchange::Receive(globalClac, 1);
    std::cout << "[Proc 2] received " << subjectFaces.size() << " subject faces (mesh1) and "
              << clippedFaces.size() << " clipped faces (mesh2)\n";

    // clac is a 1-proc communicator → BVH built and queried locally, no MPI
    FSFaceMatcher matcher(clac, subjectFaces, clippedFaces, tol);
    std::vector<FSFaceMatch> matches;
    matcher.ComputeMatches(matches);
    std::cout << "[Proc 2] computed " << matches.size() << " face matches\n";
    EXPECT_GT(matches.size(), 0u);
  }
}

_FS_END_NAMESPACE

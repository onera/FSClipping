#include "FSClipping/FSBoundaryFaceProvider.h"
#include "FSClipping/FSClippingInterfacePar.h"
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
//   Proc 2 : clipper     → Receive both → ComputeMatches → InjectTopology
//                       → BuildSurfaceTopo
//
// FSFaceExchange transmits both vertex geometry and FSFaceConnectivity
// (owner/neighbor cell IDs and types) so that proc 2 can fully assemble
// the surface topology without owning either mesh.
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
    // --- Mesh 1 proc: extract boundary faces and send geometry + topology ---
    auto mesh = LoadMeshWithClac(clac, MeshPath("input/cube_hexa_coarse.grid"));
    FSMeshFaceExtractor extractor;
    auto extraction = FSBoundaryFaceProvider::Extract(mesh, extractor, 2, tol);
    FSFaceExchange::Send(globalClac, clipperRank, extraction.faces);
    std::cout << "[Proc 0] sent " << extraction.faces.size() << " faces (mesh1) to proc 2\n";

  } else if(meshID == 1) {
    // --- Mesh 2 proc: extract boundary faces and send geometry + topology ---
    auto mesh = LoadMeshWithClac(clac, MeshPath("input/cube_hexa_fine.grid"));
    FSMeshFaceExtractor extractor;
    auto extraction = FSBoundaryFaceProvider::Extract(mesh, extractor, 1, tol);
    FSFaceExchange::Send(globalClac, clipperRank, extraction.faces);
    std::cout << "[Proc 1] sent " << extraction.faces.size() << " faces (mesh2) to proc 2\n";

  } else {
    // --- Clipper proc: receive both face sets, compute matches, build surface topo ---
    auto subjectReceived = FSFaceExchange::Receive(globalClac, 0);
    auto clippedReceived = FSFaceExchange::Receive(globalClac, 1);

    std::cout << "[Proc 2] received " << subjectReceived.faces.size()
              << " subject faces (mesh1) and " << clippedReceived.faces.size()
              << " clipped faces (mesh2)\n";

    // --- ComputeMatches (clac is 1-proc → BVH purely local) ---
    FSFaceMatcher matcher(clac, subjectReceived.faces, clippedReceived.faces, tol);
    std::vector<FSFaceMatch> matches;
    matcher.ComputeMatches(matches);
    std::cout << "[Proc 2] computed " << matches.size() << " face matches\n";
    ASSERT_GT(matches.size(), 0u);

    // --- Inject FSDM topology into matches so BuildSurfaceTopo has cell IDs ---
    FSFaceExchange::InjectTopology(matches, subjectReceived);

    // --- Build surface topology (faceKeys not needed for surface-only) ---
    FSTopologyAssembler assembler(tol);
    const std::unordered_set<GeomFaceKey, GeomFaceKeyHash> emptyFaceKeys;
    FSTopologyData surfaceTopo = assembler.BuildSurfaceTopo(matches, emptyFaceKeys);

    std::cout << "[Proc 2] surface topo: " << surfaceTopo.globalCoords.size()
              << " nodes\n";

    EXPECT_GT(surfaceTopo.globalCoords.size(), 0u);
  }
}

TEST(FSCLippingTestInterfaceParllel3, SurfaceInterface)
{
  FSClac globalClac(MPI_COMM_WORLD);
  FS_intT meshID = globalClac.GetProcID();
  FSClac clac;
  globalClac.DivideIntoGroups(meshID, clac);

  FSMesh mesh;
  const FS_floatT tol = 1e-8;
  FS_intT marker;

  if(meshID == 0) {
    mesh = LoadMeshWithClac(clac, MeshPath("input/cube_hexa_coarse.grid"));
    marker = 2;

  } else if(meshID == 1) {
    mesh = LoadMeshWithClac(clac, MeshPath("input/cube_hexa_fine.grid"));
    marker = 1;
  }

  FSClippingInterfacePar surfaceInterface(globalClac, clac, tol, marker);
  surfaceInterface.BuildSurfaceInterface(mesh);
}

_FS_END_NAMESPACE

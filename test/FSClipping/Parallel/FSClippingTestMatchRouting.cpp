// Match-routing tests runs at exactly 5 MPI procs.
// Layout: rank 0 = clipper, ranks 1..2 = mesh 1, ranks 3..4 = mesh 2.
// These tests target the gather/scatter routing path

#include "FSClipping/FSClippingEngine.h"
#include "FSClipping/FSClippingFace.h"
#include "FSClipping/FSFaceExchange.h"
#include "FSClipping/FSFaceMatcher.h"
#include "TestUtilsParallel.hpp"
#include "gtest/gtest.h"

#include <algorithm>
#include <deque>
#include <vector>

_FS_BEGIN_NAMESPACE

namespace {

// Proc layout: 1 clipper + 2 procs per mesh.
constexpr FS_intT kClipperProc = 0;
constexpr FS_intT kNProcs = 5;

// meshID in {0 = clipper, 1 = mesh A, 2 = mesh B}.
FS_intT RoutingMeshID(FS_intT procId) { return procId == kClipperProc ? 0 : (procId <= 2 ? 1 : 2); }

// Index of the proc inside its own mesh group (0 or 1).
FS_intT RankInMesh(FS_intT procId) { return RoutingMeshID(procId) == 1 ? procId - 1 : procId - 3; }

// One axis-aligned unit quad in the z = 0 plane, spanning [x0, x0+1] x [y0, y0+1].
FSFloatArrayT MakeQuad(FS_floatT x0, FS_floatT y0)
{
  FSFloatArrayT coords(FSMESH_NNODES_QUAD4, FS_3D);
  const FS_floatT pts[FSMESH_NNODES_QUAD4][FS_3D] = {
    {x0, y0, 0}, {x0 + 1, y0, 0}, {x0 + 1, y0 + 1, 0}, {x0, y0 + 1, 0}};
  for(FS_intT i = 0; i < FSMESH_NNODES_QUAD4; ++i)
    for(FS_intT d = 0; d < FS_3D; ++d)
      coords(i, d) = pts[i][d];
  return coords;
}


struct FaceSet {
  std::deque<FSFaceConnectivity> connectivity;
  std::vector<FSClippingFace> faces;
};

// Build one face whose owner cell has the given local ID, on the given proc.
void AddFace(FaceSet& set, FS_intT ownerCell, FS_intT ownerProc, FS_floatT x0, FS_floatT y0)
{
  FSCellFace owner(ownerProc, ownerCell, FSMeshEnums::CellType::CT_Hexa8, FSMeshEnums::CellType::CT_Quad4, 0);
  FSCellFace neighbor;
  set.connectivity.emplace_back(owner, neighbor);
  set.faces.emplace_back(FSFace(set.connectivity.back()), ownerCell, MakeQuad(x0, y0), -1);
}

// Geometry layout shared by the tests below.
//   mesh A, rank 0 (proc 1) : cells 7, 8   quads (0,0) (1,0)
//   mesh A, rank 1 (proc 2) : cells 7, 8   quads (0,1) (1,1)   <-- same IDs
//   mesh B, rank 0 (proc 3) : cells 7, 8   quads (0,0) (1,0)
//   mesh B, rank 1 (proc 4) : cells 7, 8   quads (0,1) (1,1)   <-- same IDs
constexpr FS_intT kCellLow = 7;
constexpr FS_intT kCellHigh = 8;

FaceSet BuildLocalFaces(FS_intT procId)
{
  const FS_intT rank = RankInMesh(procId);
  const FS_floatT y0 = (rank == 0) ? 0.0 : 1.0;

  FaceSet set;
  AddFace(set, kCellLow, procId, 0.0, y0);
  AddFace(set, kCellHigh, procId, 1.0, y0);
  return set;
}

} // namespace

TEST(FSClippingTestMatchRouting, OverlappingCellIdsBetweenProcsOfSameMesh)
{
  FSClac globalClac(MPI_COMM_WORLD);
  ASSERT_EQ(globalClac.GetNProcs(), kNProcs);

  const FS_intT procId = globalClac.GetProcID();
  const FS_intT meshID = RoutingMeshID(procId);

  FSClac localClac;
  globalClac.DivideIntoGroups(meshID, localClac);

  if(procId == kClipperProc) {
    FSClippingEngine::RunMatcherProc(globalClac, localClac, 1e-12);
    return;
  }

  FaceSet set = BuildLocalFaces(procId);
  FSFaceExchange::GatherSend(globalClac, kClipperProc, meshID, set.faces);
  std::vector<FSFaceMatch> matches = FSMatchExchange::ScatterReceive(globalClac, kClipperProc);

  EXPECT_EQ(matches.size(), 2u);
}

TEST(FSClippingTestMatchRouting, NoCellReceivesDuplicateMatches)
{
  FSClac globalClac(MPI_COMM_WORLD);
  ASSERT_EQ(globalClac.GetNProcs(), kNProcs);

  const FS_intT procId = globalClac.GetProcID();
  const FS_intT meshID = RoutingMeshID(procId);

  FSClac localClac;
  globalClac.DivideIntoGroups(meshID, localClac);

  if(procId == kClipperProc) {
    FSClippingEngine::RunMatcherProc(globalClac, localClac, 1e-12);
    return;
  }

  FaceSet set = BuildLocalFaces(procId);
  FSFaceExchange::GatherSend(globalClac, kClipperProc, meshID, set.faces);
  std::vector<FSFaceMatch> matches = FSMatchExchange::ScatterReceive(globalClac, kClipperProc);

  std::vector<FS_intT> owners;
  owners.reserve(matches.size());
  for(const auto& m : matches)
    owners.push_back(m.elemOwner1);
  std::sort(owners.begin(), owners.end());

  EXPECT_EQ(std::adjacent_find(owners.begin(), owners.end()), owners.end());
}

TEST(FSClippingTestMatchRouting, MatchesAreSplitEvenlyAcrossProcs)
{
  FSClac globalClac(MPI_COMM_WORLD);
  ASSERT_EQ(globalClac.GetNProcs(), kNProcs);

  const FS_intT procId = globalClac.GetProcID();
  const FS_intT meshID = RoutingMeshID(procId);

  FSClac localClac;
  globalClac.DivideIntoGroups(meshID, localClac);

  if(procId == kClipperProc) {
    FSClippingEngine::RunMatcherProc(globalClac, localClac, 1e-12);
    return;
  }

  FaceSet set = BuildLocalFaces(procId);
  FSFaceExchange::GatherSend(globalClac, kClipperProc, meshID, set.faces);
  std::vector<FSFaceMatch> matches = FSMatchExchange::ScatterReceive(globalClac, kClipperProc);

  FS_intT localCount = static_cast<FS_intT>(matches.size());
  FS_intT totalCount = 0;
  localClac.Sum(&localCount, &totalCount, 1);

  EXPECT_EQ(totalCount, 4);

  EXPECT_EQ(localCount, 2);
}

_FS_END_NAMESPACE

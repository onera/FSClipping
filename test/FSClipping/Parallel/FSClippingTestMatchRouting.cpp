// Match-routing tests — run at exactly 5 MPI procs.
//
// Layout: rank 0 = clipper, ranks 1..2 = mesh 1, ranks 3..4 = mesh 2.
//
// These tests target the gather/scatter routing path
// (FSFaceExchange::GatherReceiveAll + FSMatchExchange::ScatterSend), which the
// point-to-point tests in FSClippingTestMatchesSendReceive.cpp never exercise:
// there, each proc is its own mesh group and matches are sent to hard-coded
// ranks, so no routing table is built.
//
// Cell IDs are LOCAL to a proc: two procs OF THE SAME MESH may legitimately use
// the same ID for two different cells. GatherReceiveAll keys its routing table
// (GatheredFaces::cellToGlobalProc) on the bare cell ID, and the two procs of one
// mesh share a single GatheredFaces, so the second proc overwrites the first and
// every match for that ID is scattered to one (wrong) proc. Downstream,
// FSTopologyAssembler indexes cell2node on the same bare ID, merging two distinct
// cells into one — the "super Poly3D" spanning two cells.
//
// Note that IDs colliding ACROSS the two meshes are harmless: mesh A and mesh B
// get separate GatheredFaces (and therefore separate tables), so that case is
// deliberately NOT what these tests exercise.

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

// Faces carry FSDM connectivity, which FSFaceExchange::GatherSend requires
// (it dereferences topo()._faceFSDM). The connectivity objects must outlive the
// faces, since FSFace stores a raw pointer into them — hence the deque, whose
// references are stable across push_back.
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
  set.faces.emplace_back(FSFace(set.connectivity.back()), ownerCell, MakeQuad(x0, y0));
}

// Geometry layout shared by the tests below.
//
// Both meshes tile the same 2x2 unit-square region, so every face of mesh A is
// identical to exactly one face of mesh B, and each proc must get back exactly
// the matches of the 2 cells it owns.
//
// The two procs of EACH mesh deliberately reuse the SAME local cell IDs (7 and
// 8) for different cells — exactly what the asymmetric partition produces on the
// fine mesh, where procs 1 and 2 both own cells numbered 201 and 202:
//
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

// ── Overlapping cell IDs between the two procs of one mesh ───────────────────
//
// Each proc owns 2 cells, each of which has exactly one counterpart in the other
// mesh, so each proc must receive exactly 2 matches. With the routing table keyed
// on the bare cell ID, the entries for IDs 7 and 8 written by rank 0 are
// overwritten by rank 1, so all 4 matches of the mesh are scattered to rank 1:
// rank 0 receives 0 matches and rank 1 receives 4.
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

  EXPECT_EQ(matches.size(), 2u) << "proc " << procId << " owns 2 cells but received " << matches.size()
                                << " matches — matches of the other proc of the same mesh were routed here";
}

// ── No cell may collect more matches than it has interface faces ─────────────
//
// Each cell here has exactly one interface face, hence exactly one match. Two
// cells collapsing onto one ID show up as the same elemOwner1 appearing twice
// with different intersection geometry — the shape that later merges into a
// single over-sized Poly3D in FSTopologyAssembler.
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

  EXPECT_EQ(std::adjacent_find(owners.begin(), owners.end()), owners.end())
    << "proc " << procId << " received two matches for the same cell ID: two distinct cells "
    << "sharing a local ID were merged";
}

// ── Matches are distributed, not concentrated on one proc ────────────────────
//
// Summed over the 2 procs of one mesh the count is right (no match is lost — they
// are misrouted, not dropped), so the total alone cannot catch the bug. Asserting
// the total AND the per-proc split together distinguishes a routing bug from a
// matching bug: here the total stays 4 while the split degenerates to 0/4.
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

  // The 4 faces of this mesh all match, and none is lost in transit.
  EXPECT_EQ(totalCount, 4) << "mesh " << meshID << " lost or duplicated matches during scatter";

  // ...and they must not all pile up on a single proc.
  EXPECT_EQ(localCount, 2) << "proc " << procId << " holds " << localCount << " of the mesh's " << totalCount
                           << " matches instead of its own 2";
}

_FS_END_NAMESPACE

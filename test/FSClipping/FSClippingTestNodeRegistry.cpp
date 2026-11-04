#include "gtest/gtest.h"
#include <memory>
#include <vector>

#include "FSClipping/FSNodeRegistry.h"
#include "MockMixedMesh.hpp"

_FS_BEGIN_NAMESPACE

namespace {

constexpr FS_floatT kTol = 1e-6;
const FS_intT kNOld = FS_intT(mock::kCoords.size());

class FSNodeRegistryTest : public ::testing::Test
{
protected:
  FSClac clac;
  FSMesh mesh = mock::BuildMixedMesh(clac);
  std::vector<FSCellSource> sources = BuildCellSources(mesh);
  std::vector<FSFaceMatch> matches;
  BoundaryExtraction be;
  std::unique_ptr<FSCellClipped> clipped;
  std::unique_ptr<FSNodeRegistry> nodes;

  const FSNodeRegistry& Build(FS_floatT tol = kTol)
  {
    for(const auto& m : matches)
      be.cell2BoundaryFaces[m.elemOwner1].push_back(0);
    clipped = std::make_unique<FSCellClipped>(matches, be);
    nodes = std::make_unique<FSNodeRegistry>(mesh, sources, *clipped, tol);
    return *nodes;
  }

  void Add(FSMeshEnums::CellType t, std::vector<FSVec3> poly, std::vector<FS_intT> ids)
  {
    matches.push_back(mock::MakeMatch(mock::Offset(mesh, t), t, std::move(poly), std::move(ids)));
  }

  // Bottom face of the hexa split in two by the segment P(.5,0,0) - R(.5,1,0)
  void AddSplitBottomFace()
  {
    Add(FSMeshEnums::CT_Hexa8, {{0, 0, 0}, {0.5, 0, 0}, {0.5, 1, 0}, {0, 1, 0}}, {10, 100, 101, 13});
    Add(FSMeshEnums::CT_Hexa8, {{0.5, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0.5, 1, 0}}, {100, 11, 12, 101});
  }
};

std::vector<NodeId> ToVector(std::span<const NodeId> s) { return {s.begin(), s.end()}; }

} // namespace

TEST_F(FSNodeRegistryTest, NoMatchKeepsOldNodes)
{
  const auto& reg = Build();
  EXPECT_EQ(reg.Size(), kNOld);
  for(FS_intT n = 0; n < kNOld; n++) {
    EXPECT_EQ(reg.Old(n), n);
    EXPECT_TRUE(reg.Coord(n).Equal(mock::kCoords[n])) << "node " << n;
    EXPECT_EQ(reg.ClipperId(n), -1);
  }
}

TEST_F(FSNodeRegistryTest, CornerVerticesResolveToOldNodes)
{
  AddSplitBottomFace();
  const auto& reg = Build();

  EXPECT_EQ(reg.Polygon(0)[0], 0);
  EXPECT_EQ(reg.Polygon(0)[3], 3);
  EXPECT_EQ(reg.Polygon(1)[1], 1);
  EXPECT_EQ(reg.Polygon(1)[2], 2);
  EXPECT_EQ(reg.ClipperId(0), 10);
  EXPECT_EQ(reg.ClipperId(1), 11);
  EXPECT_EQ(reg.ClipperId(2), 12);
  EXPECT_EQ(reg.ClipperId(3), 13);
}

TEST_F(FSNodeRegistryTest, EdgeVerticesBecomeNewNodes)
{
  AddSplitBottomFace();
  const auto& reg = Build();

  ASSERT_EQ(reg.Size(), kNOld + 2);
  EXPECT_EQ(ToVector(reg.Polygon(0)), (std::vector<NodeId>{0, kNOld, kNOld + 1, 3}));
  EXPECT_TRUE(reg.Coord(kNOld).Equal(FSVec3(0.5, 0, 0)));
  EXPECT_TRUE(reg.Coord(kNOld + 1).Equal(FSVec3(0.5, 1, 0)));
  EXPECT_EQ(reg.ClipperId(kNOld), 100);
  EXPECT_EQ(reg.ClipperId(kNOld + 1), 101);
}

// Old nodes keep their global number, new vertices follow all old nodes, keyed by their clipper id
TEST_F(FSNodeRegistryTest, GlobalIdKeepsOldNumbersThenClipperIds)
{
  AddSplitBottomFace();
  const auto& reg = Build();
  for(NodeId n = 0; n < kNOld; n++)
    EXPECT_EQ(reg.GlobalId(n), n) << "node " << n;
  EXPECT_EQ(reg.GlobalId(kNOld), kNOld + 100);
  EXPECT_EQ(reg.GlobalId(kNOld + 1), kNOld + 101);
}

TEST_F(FSNodeRegistryTest, SharedClipperIdResolvedOnce)
{
  AddSplitBottomFace();
  const auto& reg = Build();
  EXPECT_EQ(ToVector(reg.Polygon(1)), (std::vector<NodeId>{kNOld, 1, 2, kNOld + 1}));
}

TEST_F(FSNodeRegistryTest, UntouchedNodesHaveNoClipperId)
{
  AddSplitBottomFace();
  const auto& reg = Build();
  for(FS_intT n = 4; n < kNOld; n++)
    EXPECT_EQ(reg.ClipperId(n), -1) << "node " << n;
}

TEST_F(FSNodeRegistryTest, ToleranceOnPerturbedVertex)
{
  Add(FSMeshEnums::CT_Hexa8, {{kTol / 10, 0, 0}, {1, 0, 0}, {1, 1, 0}}, {10, 11, 12});
  Add(FSMeshEnums::CT_Hexa8, {{10 * kTol, 1, 0}, {1, 1, 0}, {1, 0, 0}}, {13, 12, 11});
  const auto& reg = Build();

  EXPECT_EQ(ToVector(reg.Polygon(0)), (std::vector<NodeId>{0, 1, 2}));
  EXPECT_EQ(ToVector(reg.Polygon(1)), (std::vector<NodeId>{kNOld, 2, 1}));
  EXPECT_EQ(reg.Size(), kNOld + 1);
}

TEST_F(FSNodeRegistryTest, NearestOldNodeWins)
{
  // Nodes 0 and 1 are both within tol of (.6,0,0); node 1 is closer
  Add(FSMeshEnums::CT_Hexa8, {{0.6, 0, 0}, {1, 1, 0}, {0, 1, 0}}, {10, 12, 13});
  const auto& reg = Build(0.8);
  EXPECT_EQ(reg.Polygon(0)[0], 1);
}

TEST_F(FSNodeRegistryTest, OnlyOwnerNodesAreCandidates)
{
  // Node 8 (pyramid apex) is not a hexa node, so it is not found from a hexa match
  Add(FSMeshEnums::CT_Hexa8, {{0.5, 0.5, 1.5}, {0, 0, 1}, {1, 0, 1}}, {108, 104, 105});
  const auto& reg = Build();
  EXPECT_EQ(ToVector(reg.Polygon(0)), (std::vector<NodeId>{kNOld, 4, 5}));
}

TEST_F(FSNodeRegistryTest, OwnerPoly3DAndTetra)
{
  Add(FSMeshEnums::CT_Poly3D, {{1, 0, 1}, {1, 1, 1}, {0.5, 0.5, 1.5}}, {205, 206, 208});
  Add(FSMeshEnums::CT_Tetra4, {{1, 0, 1}, {1.5, 0.5, 1.5}, {1, 1, 1}}, {205, 209, 206});
  const auto& reg = Build();

  EXPECT_EQ(ToVector(reg.Polygon(0)), (std::vector<NodeId>{5, 6, 8}));
  EXPECT_EQ(ToVector(reg.Polygon(1)), (std::vector<NodeId>{5, 9, 6}));
  EXPECT_EQ(reg.Size(), kNOld);
  EXPECT_EQ(reg.ClipperId(9), 209);
}

TEST_F(FSNodeRegistryTest, PolygonSizes)
{
  AddSplitBottomFace();
  Add(FSMeshEnums::CT_Tetra4, {{1, 0, 1}, {1.5, 0.5, 1.5}, {1, 1, 1}}, {205, 209, 206});
  const auto& reg = Build();
  for(FS_intT m = 0; m < FS_intT(matches.size()); m++)
    EXPECT_EQ(reg.Polygon(m).size(), matches[m].clippedPoly3D.size());
}

_FS_END_NAMESPACE

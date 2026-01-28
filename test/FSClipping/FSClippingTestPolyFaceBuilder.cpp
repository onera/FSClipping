#include "FSClipping/FSPolyFaceBuilder.h"

#include <gtest/gtest.h>



TEST(FSPolyFaceBuilder, SingleCellSingleFace)
{
  // ------------------------------------------------------------
  // 1. Build a Cell2NodeBuilder
  // ------------------------------------------------------------
  FSCell2NodeBuilder cellBuilder;

  // Global cell ID = 10
  // Simple tetrahedron
  FS_intT cellId = 10;

  FSIntArrayT cell2Node(1, 4);
  cell2Node(0, 0) = 0;
  cell2Node(0, 1) = 1;
  cell2Node(0, 2) = 2;
  cell2Node(0, 3) = 3;

  FSFloatArrayT coords(4, 3);
  coords(0, 0) = 0;
  coords(0, 1) = 0;
  coords(0, 2) = 0;
  coords(1, 0) = 1;
  coords(1, 1) = 0;
  coords(1, 2) = 0;
  coords(2, 0) = 0;
  coords(2, 1) = 1;
  coords(2, 2) = 0;
  coords(3, 0) = 0;
  coords(3, 1) = 0;
  coords(3, 2) = 1;

  cellBuilder.AddOldCellNodes(
    cellId,
    FSMeshEnums::CellType::CT_Tetra4
    // cell2Node, coords
  );

  // ------------------------------------------------------------
  // 2. Create a face match with a clipped polygon (triangle)
  // ------------------------------------------------------------
  FSFaceMatch match;
  match.type = FSFaceMatch::INTERSECTING;
  match.elemOwner1 = cellId;

  match.clippedPoly3D = {
    FSVec3(0, 0, 0),
    FSVec3(1, 0, 0),
    FSVec3(0, 1, 0)};

  std::vector<FSFaceMatch> matches = {match};

  cellBuilder.AddClippedPolygon(match.elemOwner1, match.clippedPoly3D);
  cellBuilder.BuildGlobalNumbering();

  // ------------------------------------------------------------
  // 3. Build polyhedral faces
  // ------------------------------------------------------------
  FSPolyFaceBuilder faceBuilder(matches, cellBuilder);

  FSMeshPolyFaceStorage polyFaces;
  faceBuilder.Build(reinterpret_cast<FSMeshPolyFaceStorage&>(polyFaces));

  // ------------------------------------------------------------
  // 4. Checks
  // ------------------------------------------------------------

  // Single cell
  EXPECT_EQ(polyFaces.GetNumPolys(), 1);

  // Single face
  EXPECT_EQ(polyFaces.GetNFaces(0), 1);

  // Triangular face
  EXPECT_EQ(polyFaces.GetNFaceNodes(0, 0), 3);

  // Expected local node indices
  const auto& nodes = polyFaces.GetFaceNodes(0, 0);

  ASSERT_EQ(nodes.Size(), 3);

  // Must match the local indices of the tetrahedron
  EXPECT_EQ(nodes[0], 0);
  EXPECT_EQ(nodes[1], 1);
  EXPECT_EQ(nodes[2], 2);
}

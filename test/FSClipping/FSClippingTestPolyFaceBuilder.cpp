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

  cellBuilder.AddCellNodes(
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
  std::unordered_set<GeomFaceKey, GeomFaceKeyHash> faceKeys;

  cellBuilder.AddClippedPolygon(match.elemOwner1, match.clippedPoly3D);
  cellBuilder.BuildGlobalNumbering();

  // ------------------------------------------------------------
  // 3. Build polyhedral faces
  // ------------------------------------------------------------
  FSPolyFaceBuilder faceBuilder(matches, cellBuilder, faceKeys);
  faceBuilder.CollectMatchesFaces();

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

TEST(FSPolyFaceBuilder, TwoCellsSharedFace)
{
  FSCell2NodeBuilder cellBuilder;

  // Cell 1
  FS_intT c1 = 1;
  cellBuilder.AddCellNodes(c1, FSMeshEnums::CellType::CT_Tetra4);

  // Cell 2
  FS_intT c2 = 2;
  cellBuilder.AddCellNodes(c2, FSMeshEnums::CellType::CT_Tetra4);

  std::vector<FSVec3> face = {
    FSVec3(0, 0, 0),
    FSVec3(1, 0, 0),
    FSVec3(0, 1, 0)};

  FSFaceMatch m1, m2;
  m1.type = FSFaceMatch::INTERSECTING;
  m1.elemOwner1 = c1;
  m1.clippedPoly3D = face;

  m2.type = FSFaceMatch::INTERSECTING;
  m2.elemOwner1 = c2;
  m2.clippedPoly3D = face;

  std::vector<FSFaceMatch> matches = {m1, m2};
  std::unordered_set<GeomFaceKey, GeomFaceKeyHash> keys;

  cellBuilder.AddClippedPolygon(c1, face);
  cellBuilder.AddClippedPolygon(c2, face);
  cellBuilder.BuildGlobalNumbering();

  FSPolyFaceBuilder builder(matches, cellBuilder, keys);

  builder.CollectMatchesFaces();

  const auto& cellFaces = builder.CellFaces();

  ASSERT_EQ(cellFaces.size(), 2);
  EXPECT_EQ(cellFaces[0].size(), 1);
  EXPECT_EQ(cellFaces[1].size(), 1);

  const auto& f0 = cellFaces[0][0].nodeIds;
  const auto& f1 = cellFaces[1][0].nodeIds;

  EXPECT_EQ(f0.size(), f1.size());
  for(size_t i = 0; i < f0.size(); ++i)
    EXPECT_EQ(f0[i], f1[i]);
}

TEST(FSPolyFaceBuilder, FaceReorienting)
{
  FSCell2NodeBuilder cellBuilder;

  FS_intT cellId = 0;

  // Tetra explicite
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

  cellBuilder.AddVolumeCellNodes(cellId,
                                 FSMeshEnums::CellType::CT_Tetra4,
                                 cell2Node,
                                 coords);

  std::vector<FSVec3> faceCorrect = {
    FSVec3(0, 0, 0),
    FSVec3(1, 0, 0),
    FSVec3(0, 1, 0)};

  std::vector<FSVec3> faceWrong = {
    FSVec3(0, 1, 0),
    FSVec3(0, 0, 0),
    FSVec3(1, 0, 0)};

  FSFaceMatch match;
  match.type = FSFaceMatch::INTERSECTING;
  match.elemOwner1 = cellId;
  match.clippedPoly3D = faceWrong;

  std::vector<FSFaceMatch> matches = {match};
  std::unordered_set<GeomFaceKey, GeomFaceKeyHash> keys;

  cellBuilder.AddClippedPolygon(cellId, faceWrong);
  cellBuilder.BuildGlobalNumbering();

  FSPolyFaceBuilder builder(matches, cellBuilder, keys);

  builder.CollectMatchesFaces();

  auto before = builder.CellFaces()[0][0].nodeIds;

  builder.Reorienting();

  auto after = builder.CellFaces()[0][0].nodeIds;

  ASSERT_EQ(before.size(), after.size());

  bool reversed = true;
  for(size_t i = 0; i < before.size(); ++i) {
    if(before[i] != after[before.size() - 1 - i]) {
      reversed = false;
      break;
    }
  }

  EXPECT_TRUE(reversed);
}
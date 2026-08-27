#include "FSClipping/FSCell2NodeBuilder.h"

#include <gtest/gtest.h>

TEST(FSCell2NodeBuilder, SingleCellWithClippedPolygon)
{
  FSCell2NodeBuilder builder(1e-12);

  // Identity square cell
  FSIntArrayT cell2Node(1, 8);
  for(int i = 0; i < 8; ++i)
    cell2Node(0, i) = i;

  FSFloatArrayT coords(8, 3);
  coords(0, 0) = 0;
  coords(0, 1) = 0;
  coords(0, 2) = 0;

  coords(1, 0) = 1;
  coords(1, 1) = 0;
  coords(1, 2) = 0;

  coords(2, 0) = 1;
  coords(2, 1) = 1;
  coords(2, 2) = 0;

  coords(3, 0) = 0;
  coords(3, 1) = 1;
  coords(3, 2) = 0;

  coords(4, 0) = 0;
  coords(4, 1) = 0;
  coords(4, 2) = 1;

  coords(5, 0) = 1;
  coords(5, 1) = 0;
  coords(5, 2) = 1;

  coords(6, 0) = 1;
  coords(6, 1) = 1;
  coords(6, 2) = 1;

  coords(7, 0) = 0;
  coords(7, 1) = 1;
  coords(7, 2) = 1;

  builder.AddCellNodes(0, FSMeshEnums::CellType::CT_Hexa8);

  // Clipped polygon on the plane z = 1
  std::vector<FSVec3> poly = {FSVec3(0.5, 0.0, 1.0), FSVec3(1.0, 0.5, 1.0), FSVec3(0.5, 1.0, 1.0),
                              FSVec3(0.0, 0.5, 1.0)};

  builder.AddClippedPolygon(0, poly);
  builder.BuildGlobalNumbering();

  const auto& cellData = builder.CellData();
  ASSERT_EQ(cellData.size(), 1);

  const auto& data = cellData.at(0);
  EXPECT_EQ(data.coords.size(), 4); // 4 nodes after clipping on the border
  EXPECT_EQ(data.nodeIds.size(), 4);

  const auto& globalCoords = builder.GlobalCoords();
  EXPECT_EQ(globalCoords.size(), 4);

  // Check unicity
  for(size_t i = 0; i < globalCoords.size(); ++i)
    for(size_t j = i + 1; j < globalCoords.size(); ++j)
      EXPECT_GT((globalCoords[i] - globalCoords[j]).L2Norm(), 1e-12);
}

TEST(FSCell2NodeBuilder, SharedNodesBetweenCells)
{
  FSCell2NodeBuilder builder(1e-12);

  std::vector<FSVec3> poly = {FSVec3(0.0, 0.0, 0.0), FSVec3(1.0, 0.0, 0.0), FSVec3(1.0, 1.0, 0.0),
                              FSVec3(0.0, 1.0, 0.0)};

  builder.AddCellNodes(0, FSMeshEnums::CellType::CT_Hexa8);
  builder.AddCellNodes(1, FSMeshEnums::CellType::CT_Hexa8);

  builder.AddClippedPolygon(0, poly);
  builder.AddClippedPolygon(1, poly);

  builder.BuildGlobalNumbering();

  const auto& globalCoords = builder.GlobalCoords();

  EXPECT_EQ(globalCoords.size(), 4);
}

TEST(FSCell2NodeBuilder, NodeMergingWithTolerance)
{
  FSCell2NodeBuilder builder(1e-6);

  std::vector<FSVec3> poly = {FSVec3(0.0, 0.0, 0.0), FSVec3(1.0, 0.0, 0.0), FSVec3(1.0, 1.0, 0.0),
                              FSVec3(0.0, 1.0, 0.0)};

  std::vector<FSVec3> polyPerturbed = {FSVec3(0.0 + 1e-8, 0.0, 0.0), FSVec3(1.0, 0.0 + 1e-8, 0.0),
                                       FSVec3(1.0, 1.0, 0.0), FSVec3(0.0, 1.0, 0.0)};

  builder.AddCellNodes(0, FSMeshEnums::CellType::CT_Hexa8);
  builder.AddCellNodes(1, FSMeshEnums::CellType::CT_Hexa8);


  builder.AddClippedPolygon(0, poly);
  builder.AddClippedPolygon(1, polyPerturbed);

  builder.BuildGlobalNumbering();

  EXPECT_EQ(builder.GlobalCoords().size(), 4);
}

TEST(FSCell2NodeBuilder, CellIdMapping)
{
  FSCell2NodeBuilder builder;

  std::vector<FSVec3> poly1 = {FSVec3(0.0, 0.0, 0.0), FSVec3(1.0, 0.0, 0.0), FSVec3(1.0, 1.0, 0.0)};

  std::vector<FSVec3> poly2 = {FSVec3(0.0, 0.0, 0.0), FSVec3(-1.0, 0.0, 0.0), FSVec3(-1.0, -1.0, 0.0)};


  builder.AddCellNodes(42, FSMeshEnums::CellType::CT_Tri3);
  builder.AddCellNodes(7, FSMeshEnums::CellType::CT_Tri3);


  builder.AddClippedPolygon(42, poly1);
  builder.AddClippedPolygon(7, poly2);

  builder.BuildGlobalNumbering();

  auto l0 = builder.LocalCellIndex(42);
  auto l1 = builder.LocalCellIndex(7);

  EXPECT_NE(l0, l1);
  EXPECT_EQ(builder.GlobalCellId(l0), 42);
  EXPECT_EQ(builder.GlobalCellId(l1), 7);
}

TEST(FSCell2NodeBuilder, LocalIndexConsistency)
{
  FSCell2NodeBuilder builder;

  std::vector<FSVec3> poly = {FSVec3(0, 0, 0), FSVec3(1, 0, 0), FSVec3(1, 1, 0)};

  builder.AddCellNodes(0, FSMeshEnums::CellType::CT_Tri3);
  builder.AddClippedPolygon(0, poly);
  builder.BuildGlobalNumbering();

  const auto& data = builder.CellData().at(0);

  for(const auto& p : data.coords) {
    NodeKey key(p, 1e-12);
    EXPECT_TRUE(data.localIndex.find(key) != data.localIndex.end());
  }
}

TEST(NodeGlobalNumbers, MapsLocalIndicesToClipperIds)
{
  FSCell2NodeBuilder builder(1e-12);

  const std::vector<FSVec3> poly = {FSVec3(0, 0, 0), FSVec3(1, 0, 0), FSVec3(0, 1, 0)};
  const std::vector<FS_intT> clipperIds = {42, 7, 13};

  builder.AddCellNodes(0, FSMeshEnums::CellType::CT_Hexa8);
  builder.AddClippedPolygon(0, poly, clipperIds);
  builder.BuildGlobalNumbering();

  const FSIntArrayT numbers = builder.NodeGlobalNumbers();
  ASSERT_EQ(numbers.Size(), 3);

  // For every local node, the returned ID must be the clipper ID of the
  // geometrically matching point.
  const auto& coords = builder.GlobalCoords();
  for(FS_intT l = 0; l < numbers.Size(); ++l) {
    bool found = false;
    for(std::size_t i = 0; i < poly.size(); ++i) {
      if((coords[l] - poly[i]).L2Norm() < 1e-12) {
        EXPECT_EQ(numbers[l], clipperIds[i]);
        found = true;
      }
    }
    EXPECT_TRUE(found);
  }
}

TEST(NodeGlobalNumbers, SharedNodeKeepsSingleId)
{
  FSCell2NodeBuilder builder(1e-12);

  // Two cells sharing the point (1,0,0) with ID 99 in both.
  builder.AddCellNodes(0, FSMeshEnums::CellType::CT_Hexa8);
  builder.AddCellNodes(1, FSMeshEnums::CellType::CT_Hexa8);
  builder.AddClippedPolygon(0, {FSVec3(0, 0, 0), FSVec3(1, 0, 0)}, {1, 99});
  builder.AddClippedPolygon(1, {FSVec3(1, 0, 0), FSVec3(2, 0, 0)}, {99, 2});
  builder.BuildGlobalNumbering();

  const FSIntArrayT numbers = builder.NodeGlobalNumbers();
  ASSERT_EQ(numbers.Size(), 3); // 3 unique points

  FS_intT count99 = 0;
  for(FS_intT l = 0; l < numbers.Size(); ++l)
    if(numbers[l] == 99)
      ++count99;
  EXPECT_EQ(count99, 1);
}
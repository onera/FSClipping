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

  builder.AddOldCellNodes(
    0, FSMeshEnums::CellType::CT_Hexa8);

  // Clipped polygon on the plane z = 1
  std::vector<FSVec3> poly = {
    FSVec3(0.5, 0.0, 1.0),
    FSVec3(1.0, 0.5, 1.0),
    FSVec3(0.5, 1.0, 1.0),
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

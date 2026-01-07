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

  builder.addOldCellNodes(
    0, FSMeshEnums::CellType::CT_Hexa8, cell2Node, coords);

  // Clipped polygon on the plane z = 1
  std::vector<FSVec3> poly = {
    FSVec3(0.5, 0.0, 1.0),
    FSVec3(1.0, 0.5, 1.0),
    FSVec3(0.5, 1.0, 1.0),
    FSVec3(0.0, 0.5, 1.0)};

  builder.addClippedPolygon(0, poly);
  builder.buildGlobalNumbering();

  const auto& cellData = builder.cellData();
  ASSERT_EQ(cellData.size(), 1);

  const auto& data = cellData.at(0);
  EXPECT_EQ(data.coords.size(), 12); // 8 old nodes + 4 new nodes
  EXPECT_EQ(data.nodeIds.size(), 12);

  const auto& globalCoords = builder.globalCoords();
  EXPECT_EQ(globalCoords.size(), 12);

  // Vérifier unicité
  for(size_t i = 0; i < globalCoords.size(); ++i)
    for(size_t j = i + 1; j < globalCoords.size(); ++j)
      EXPECT_GT((globalCoords[i] - globalCoords[j]).L2Norm(), 1e-12);
}

#ifndef FSCELL2NODEBUILDER_H
#define FSCELL2NODEBUILDER_H

#include <FSArray.h>
#include <FSTypes.h>
#include <FSVec3.h>
#include "FSMeshEnums.h"

struct Cell2NodeData {
  std::vector<FSVec3> coords;
  std::vector<FS_intT> nodeIds;
};

class FSCell2NodeBuilder
{
public:
  explicit FSCell2NodeBuilder(FS_floatT tol = 1e-12) : tol_(tol) {};

  void addOldCellNodes(FS_intT cellId,
                       FSMeshEnums::CellType cellType,
                       const FSIntArrayT& cell2Node,
                       const FSFloatArrayT& coords);

  void addClippedPolygon(FS_intT cellId,
                         const std::vector<FSVec3>& poly); // we need to check if the nodes is not already in cellData_

  void buildGlobalNumbering();

  const std::unordered_map<FS_intT, Cell2NodeData>& cellData() const noexcept { return cellData_; };

  const std::vector<FSVec3>& globalCoords() const noexcept { return globalCoords_; };

  FSIntArrayT cell2Node();

private:
  FS_floatT tol_;
  std::unordered_map<FS_intT, Cell2NodeData> cellData_;
  std::vector<FSVec3> globalCoords_; // the list of all the coords

  bool exists(const std::vector<FSVec3>& nodes,
              const FSVec3& p) const;
};


#endif
#include "FSClipping/FSCell2NodeBuilder.h"
#include "FSCellInfo.h"

bool FSCell2NodeBuilder::exists(const std::vector<FSVec3>& nodes,
                                const FSVec3& p) const
{
  for(const auto& q : nodes)
    if((p - q).L2Norm() <= tol_)
      return true;
  return false;
}

void FSCell2NodeBuilder::addOldCellNodes(FS_intT cellId,
                                         FSMeshEnums::CellType cellType,
                                         const FSIntArrayT& cell2Node,
                                         const FSFloatArrayT& coords)
{
  auto [it, inserted] = cellData_.try_emplace(cellId);

  if(!inserted)
    return;

  auto& nodes = it->second.coords;

  const FS_intT nCellNodes = FSCellInfo::cNNodes[cellType];
  nodes.reserve(nCellNodes);

  for(FS_intT node = 0; node < nCellNodes; ++node) {
    FS_intT idx = cell2Node(cellId, node);
    nodes.emplace_back(coords(idx, 0),
                       coords(idx, 1),
                       coords(idx, 2));
  }
}


void FSCell2NodeBuilder::addClippedPolygon(
  FS_intT cellId,
  const std::vector<FSVec3>& poly)
{
  auto& nodes = cellData_.at(cellId).coords;

  for(const auto& p : poly) {
    if(!exists(nodes, p))
      nodes.push_back(p);
  }
}

void FSCell2NodeBuilder::buildGlobalNumbering()
{
  std::vector<FSVec3> uniqueCoords;
  std::vector<FS_intT> ids;

  for(auto& [cellId, data] : cellData_) {
    data.nodeIds.resize(data.coords.size());

    for(size_t i = 0; i < data.coords.size(); ++i) {
      const auto& p = data.coords[i];

      auto it = std::find_if(
        uniqueCoords.begin(), uniqueCoords.end(),
        [&](const FSVec3& q) {
          return (p - q).L2Norm() < tol_;
        }); // O(n^2) complexity

      if(it == uniqueCoords.end()) {
        FS_intT newId = uniqueCoords.size();
        uniqueCoords.push_back(p);
        data.nodeIds[i] = newId;
      } else {
        data.nodeIds[i] = std::distance(uniqueCoords.begin(), it);
      }
    }
  }

  globalCoords_ = std::move(uniqueCoords);
}

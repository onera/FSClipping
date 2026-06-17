#include "FSClipping/FSCell2NodeBuilder.h"
#include "FSCellInfo.h"

NodeKey::NodeKey(const FSVec3& p, FS_floatT tol)
{
  ix = static_cast<FS_intT>(std::llround(p[0] / tol));
  iy = static_cast<FS_intT>(std::llround(p[1] / tol));
  iz = static_cast<FS_intT>(std::llround(p[2] / tol));
}

void FSCell2NodeBuilder::SetAllCellOnTheBorder()
{
  for(auto& cellData : cellData_)
    cellData.second.onTheBorder = true;
  numCellOnTheBorder_ = cellData_.size();
}

bool FSCell2NodeBuilder::exists(const std::vector<FSVec3>& nodes, const FSVec3& p) const
{
  for(const auto& q : nodes)
    if((p - q).L2Norm() <= tol_)
      return true;
  return false;
}

void FSCell2NodeBuilder::BuildCellIdMapping()
{
  const FS_intT nCells = static_cast<FS_intT>(cellData_.size());
  cellIds_.clear();
  cellIds_.reserve(nCells);

  for(const auto& [cellId, _] : cellData_)
    cellIds_.push_back(cellId);

  std::sort(cellIds_.begin(), cellIds_.end());

  cellId2L_.clear();
  cellId2L_.reserve(nCells);

  for(FS_intT i = 0; i < nCells; ++i)
    cellId2L_.emplace(cellIds_[i], i);
}

FS_intT FSCell2NodeBuilder::LocalCellIndex(FS_intT globalId) const
{
  auto it = cellId2L_.find(globalId);
  if(it == cellId2L_.end())
    FSError.SetAndPrintAndExit("FSCell2NodeBuilder : Unknown cellId");
  return it->second;
}

FS_intT FSCell2NodeBuilder::GlobalCellId(FS_intT localId) const { return cellIds_[localId]; }

void FSCell2NodeBuilder::AddCellNodes(FS_intT cellId, FSMeshEnums::CellType cellType)
{
  auto [it, inserted] = cellData_.try_emplace(cellId);

  if(!inserted)
    return;

  auto& nodes = it->second.coords;

  const FS_intT nCellNodes = FSCellInfo::cNNodes[cellType];
  nodes.reserve(nCellNodes);
}

void FSCell2NodeBuilder::AddClippedPolygon(FS_intT cellId, const std::vector<FSVec3>& poly)
{
  auto it = cellData_.find(cellId);
  if(it == cellData_.end())
    FSError.SetAndPrintAndExit(
      "FSCell2NodeBuilder::AddClippedPolygon You are trying to add a clipped polygon to an inexistent cell.");

  auto& nodes = it->second.coords;
  for(const auto& p : poly) {
    if(!exists(nodes, p))
      nodes.push_back(p);
  }
}

void FSCell2NodeBuilder::BuildGlobalNumbering()
{
  // Pass 1: assign a unique global index to each coordinate.
  // Builds globalCoords_ (deduplicated list) and coordToNode_ (O(1) lookup),
  // then stores the resulting index in cell2Node.nodeIds.
  for(auto& [_, data] : cellData_) {
    data.nodeIds.resize(data.coords.size());

    for(std::size_t i = 0; i < data.coords.size(); ++i) {
      const auto& p = data.coords[i];
      NodeKey key(p, tol_);

      auto it = coordToNode_.find(key);
      if(it == coordToNode_.end()) {
        FS_intT id = globalCoords_.size();
        globalCoords_.push_back(p);
        coordToNode_[key] = id;
        data.nodeIds[i] = id;
      } else {
        data.nodeIds[i] = it->second;
      }
    }
  }

  // Pass 2: build the per-cell coord→localIndex map used by FSPolyFaceBuilder
  // to resolve face node indices within the cell's local coordinate array.
  for(auto& [cellId, data] : cellData_) {
    data.localIndex.reserve(data.coords.size());

    for(std::size_t i = 0; i < data.coords.size(); ++i) {
      NodeKey key(data.coords[i], tol_);
      data.localIndex.emplace(key, i);
    }
  }
  BuildCellIdMapping();
}

void FSCell2NodeBuilder::AddVolumeCellNodes(FS_intT cellId, FSMeshEnums::CellType cellType,
                                            const FSIntArrayT& cell2Node, const FSFloatArrayT& coords)
{
  auto [it, inserted] = cellData_.try_emplace(cellId);

  auto& nodes = it->second.coords;

  if(nodes.empty() && !inserted)
    FSError.SetAndPrintAndExit("FSCell2NodeBuilder::AddVolumeCellNodes You have an empty array of nodes. The match "
                               "builder should have added some nodes.");

  const FS_intT nCellNodes = FSCellInfo::cNNodes[cellType];

  for(FS_intT node = 0; node < nCellNodes; ++node) {
    FS_intT idx = cell2Node(cellId, node);
    auto x = coords(idx, 0);
    auto y = coords(idx, 1);
    auto z = coords(idx, 2);
    FSVec3 vecNode{x, y, z};
    NodeKey key(vecNode, tol_);
    if(!it->second.localIndex.contains(key))
      nodes.emplace_back(coords(idx, 0), coords(idx, 1), coords(idx, 2));
  }
}

FSIntRegisterT FSCell2NodeBuilder::Cell2NodePoly3D()
{
  FSIntRegisterT cell2Node;
  const FS_intT nCells = static_cast<FS_intT>(cellData_.size());
  cell2Node.Init(numCellOnTheBorder_);
  FS_intT iterCellBorder = 0;

  // Pass 1: count nodes per border cell
  for(FS_intT localIdx = 0; localIdx < nCells; ++localIdx) {
    const auto& data = cellData_.at(GlobalCellId(localIdx));
    if(data.onTheBorder) {
      cell2Node.Count(iterCellBorder, static_cast<FS_intT>(data.nodeIds.size()));
      iterCellBorder++;
    }
  }
  cell2Node.Prepare();
  iterCellBorder = 0;

  // Pass 2: fill node indices
  for(FS_intT localIdx = 0; localIdx < nCells; ++localIdx) {
    const auto& data = cellData_.at(GlobalCellId(localIdx));
    if(data.onTheBorder) {
      for(FS_intT nodeId : data.nodeIds)
        cell2Node.Add(iterCellBorder, nodeId);
      iterCellBorder++;
    }
  }
  return cell2Node;
}

FSIntArrayT FSCell2NodeBuilder::Cell2NodeInner(FSMeshEnums::CellType type)
{
  FSIntArrayT cell2Node;
  const FS_intT nCells = static_cast<FS_intT>(cellData_.size());
  FS_intT i = 0;
  cell2Node.Resize(nCells - numCellOnTheBorder_, FSCellInfo::GetNNodes(type));

  for(FS_intT localIdx = 0; localIdx < nCells; ++localIdx) {
    FS_intT globalId = GlobalCellId(localIdx);
    if(!cellData_[globalId].onTheBorder) {
      const auto& data = cellData_.at(GlobalCellId(localIdx));
      for(std::size_t j = 0; j < data.nodeIds.size(); j++)
        cell2Node(i, j) = data.nodeIds[j];
      i++;
    }
  }

  return cell2Node;
}

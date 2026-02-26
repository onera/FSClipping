#include "FSClipping/FSCell2NodeBuilder.h"
#include "FSCellInfo.h"

NodeKey::NodeKey(const FSVec3& p,
                 FS_floatT tol)
{
  ix = static_cast<FS_intT>(std::llround(p[0] / tol));
  iy = static_cast<FS_intT>(std::llround(p[1] / tol));
  iz = static_cast<FS_intT>(std::llround(p[2] / tol));
}


bool FSCell2NodeBuilder::exists(const std::vector<FSVec3>& nodes,
                                const FSVec3& p) const
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

FS_intT FSCell2NodeBuilder::GlobalCellId(FS_intT localId) const
{
  return cellIds_[localId];
}

void FSCell2NodeBuilder::AddCellNodes(FS_intT cellId,
                                      FSMeshEnums::CellType cellType)
{
  auto [it, inserted] = cellData_.try_emplace(cellId);

  if(!inserted)
    return;

  auto& nodes = it->second.coords;

  const FS_intT nCellNodes = FSCellInfo::cNNodes[cellType];
  nodes.reserve(nCellNodes);
}

void FSCell2NodeBuilder::AddClippedPolygon(FS_intT cellId,
                                           const std::vector<FSVec3>& poly)
{

  auto& nodes = cellData_.at(cellId).coords;

  for(const auto& p : poly) {
    if(!exists(nodes, p)) // we need to check if the nodes is not already in cellData_
      nodes.push_back(p);
  }
}

void FSCell2NodeBuilder::BuildGlobalNumbering()
{
  // First part :
  // We restart the global numbering of the nodes. So the first node
  // that we met will be at the 0 indice in the globalCoords_ vector.
  // We also construct the coordToNode_ mapping which is the mapping
  // coord -> localIndex in GlobalCoords container.
  // And finally, we construct the cell2Node connectiviy with the new Id.
  for(auto& [_, data] : cellData_) {
    data.nodeIds.resize(data.coords.size());

    for(std::size_t i = 0; i < data.coords.size(); ++i) {
      const auto& p = data.coords[i];
      NodeKey key(p, tol_);

      auto it = coordToNode_.find(key); // O(1) the reason of the NodeKey struct
      if(it == coordToNode_.end()) {
        FS_intT id = globalCoords_.size(); // new node id corresponding to the actual size of globalCoords
        globalCoords_.push_back(p);
        coordToNode_[key] = id; // add the local index in the globalCoords container
        data.nodeIds[i] = id;   // same
      } else {
        data.nodeIds[i] = it->second; // we add to cell 2 node a node already existing in globalCoords/coordToNode
      }
    }
  }

  // Second part :
  // We keep the local/inter index of the nodes within the container
  // cellData_[elem].nodeIds. Very important later for
  // the FSPolyFaceBuilder and the face2node container
  for(auto& [cellId, data] : cellData_) {
    data.localIndex.reserve(data.coords.size());

    for(std::size_t i = 0; i < data.coords.size(); ++i) {
      NodeKey key(data.coords[i], tol_);
      data.localIndex.emplace(key, i);
    }
  }
  BuildCellIdMapping();
}

void FSCell2NodeBuilder::AddVolumeCellNodes(FS_intT cellId,
                                            FSMeshEnums::CellType cellType,
                                            const FSIntArrayT& cell2Node,
                                            const FSFloatArrayT& coords)
{
  auto [it, inserted] = cellData_.try_emplace(cellId);

  if(inserted)
    FSError.SetAndPrintAndExit("FSCell2NodeBuilder::AddOldCellNodes You're trying to insert a cell which is not on the border.");

  auto& nodes = it->second.coords;

  if(nodes.empty())
    FSError.SetAndPrintAndExit("FSCell2NodeBuilder::AddOldCellNodes You have an empty array of nodes. The match builder should have added some nodes.");

  // std::cout << nodes.size() << " surface node already add for the cell " << cellId << " : " << "\n";
  // for(const auto& n : nodes)
  //   std::cout << n[0] << " " << n[1] << " " << n[2] << "\t";
  // std::cout << "\n";
  const FS_intT nCellNodes = FSCellInfo::cNNodes[cellType];
  // nodes.reserve(nCellNodes);

  // This part add the coordinates if the nodes in the element border of the previous mesh.
  // We don't need it if we want just to reconstruct only the new faces border.
  for(FS_intT node = 0; node < nCellNodes; ++node) {
    FS_intT idx = cell2Node(cellId, node);
    auto x = coords(idx, 0);
    auto y = coords(idx, 1);
    auto z = coords(idx, 2);
    FSVec3 vecNode{x, y, z};
    NodeKey key(vecNode, tol_);
    if(!(it->second.localIndex.contains(key))) { // we check that the coord does not already exist
      nodes.emplace_back(coords(idx, 0),
                         coords(idx, 1),
                         coords(idx, 2));
    }
  }

  // std::cout << nodes.size() << " volume node already add for the cell " << cellId << " : " << "\n";
  // for(const auto& n : nodes)
  //   std::cout << n[0] << " " << n[1] << " " << n[2] << "\t";
  // std::cout << "\n";
}

void FSCell2NodeBuilder::AddVolumeCellNodesInner(FS_intT cellId,
                                                 FSMeshEnums::CellType cellType,
                                                 const FSIntArrayT& cell2Node,
                                                 const FSFloatArrayT& coords)
{

  auto [it, inserted] = cellData_.try_emplace(cellId);

  if(!inserted)
    FSError.SetAndPrintAndExit("FSCell2NodeBuilder::AddOldCellNodesInner You're trying to insert a cell which is in the border.");

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

FSIntRegisterT FSCell2NodeBuilder::Cell2Node()
{
  // Return the final cell2Node connectivity for FSDM, in particular
  // to fit with the parameter of FSMeshData::InitUnstructCells.
  FSIntRegisterT cell2Node;
  // const FS_intT nCells = static_cast<FS_intT>(numCellsBorder_);
  const FS_intT nCells = static_cast<FS_intT>(cellData_.size());
  cell2Node.Init(nCells);

  // --- PASS 1 : Count ---
  for(FS_intT localIdx = 0; localIdx < nCells; ++localIdx) {
    const auto& data = cellData_.at(GlobalCellId(localIdx));
    cell2Node.Count(localIdx,
                    static_cast<FS_intT>(data.nodeIds.size()));
  }

  cell2Node.Prepare();

  // --- PASS 2 : Add ---
  for(FS_intT localIdx = 0; localIdx < nCells; ++localIdx) {
    const auto& data = cellData_.at(GlobalCellId(localIdx));
    for(FS_intT nodeId : data.nodeIds)
      cell2Node.Add(localIdx, nodeId);
  }

  return cell2Node;
}

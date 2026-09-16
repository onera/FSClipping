#include "FSClipping/FSMeshReconstruction.h"
#include <FSMeshData.h>

#include <unordered_map>

FSMesh FSMeshReconstruction::Build(const FSUnstructMeshData& meshDataOriginal, FSTopologyData& topologyData)
{
  FSMesh meshClipped = FSMesh(&clac_);

  const auto& globalCoords = topologyData.globalCoords;
  auto& cell2NodePoly2D = topologyData.cell2NodePoly2D;
  auto& cell2NodePoly3D = topologyData.cell2NodePoly3D;
  auto& cell2NodeInner = topologyData.cell2NodeInner;
  auto& polyFaces = topologyData.polyFaces;

  // Remap the global numbering of the node if distributed
  const bool distributed = clac_.GetNProcs() > 1;
  std::vector<FS_intT> ownedLocalRows;
  FSIntArrayT ownedGlobalNumbers;

  if(distributed)
    RemapToDistributedNumbering(topologyData, ownedLocalRows, ownedGlobalNumbers);
  else {
    ownedLocalRows.resize(globalCoords.size());
    for(std::size_t i = 0; i < globalCoords.size(); ++i)
      ownedLocalRows[i] = static_cast<FS_intT>(i);
  }
  const FS_intT nOwned = static_cast<FS_intT>(ownedLocalRows.size());

  // Begin the initialization of the new clippedMesh
  meshClipped.BeginInitialization();
  FSMeshData* meshDataClippedPtr = meshClipped.GetMeshData();
  FSUnstructMeshData& unstructMeshData = meshDataClippedPtr->GetUnstructCells();

  meshClipped.InitUnstructNodes(nOwned);
  if(!cell2NodeInner.IsEmpty()) {
    for(auto& [type, cell2Node] : topologyData.cell2NodeInner)
      meshClipped.InitUnstructCells(type, cell2Node);
  }

  meshDataClippedPtr->InitUnstructCells(FSMeshEnums::CT_Poly2D, cell2NodePoly2D);
  meshDataClippedPtr->InitUnstructCells(FSMeshEnums::CT_Poly3D, cell2NodePoly3D);
  meshDataClippedPtr->InitUnstructCellFaces(FSMeshEnums::CT_Poly3D, polyFaces);
  meshClipped.EndInitialization();

  FSQuantityDescArrayT coordsDesc(3);
  coordsDesc[0] = FSQuantityDesc(FSDataName::Coordinates(), FSDataName::Coordinate().X());
  coordsDesc[1] = FSQuantityDesc(FSDataName::Coordinates(), FSDataName::Coordinate().Y());
  coordsDesc[2] = FSQuantityDesc(FSDataName::Coordinates(), FSDataName::Coordinate().Z());

  FSFloatArrayT coordInterface(nOwned, 3);
  for(FS_intT i = 0; i < nOwned; i++) {
    const auto& p = globalCoords[ownedLocalRows[i]];
    coordInterface(i, 0) = p[0];
    coordInterface(i, 1) = p[1];
    coordInterface(i, 2) = p[2];
  }

  bool success = unstructMeshData.SetCoordinates3D(coordsDesc, coordInterface);
  if(!success)
    FSError.SetAndPrintAndExit("FSMeshReconstruction : Error while set coordinates");

  const FSString globalNumberName = FSMeshEnums::AttributeTypeToString(FSMeshEnums::AT_GlobalNumber);

  if(distributed)
    unstructMeshData.InitCellAttribute(globalNumberName, FSMeshEnums::CT_Node, ownedGlobalNumbers);
  else {
    FS_intT currentOffset = 0;
    FSIntArrayT cellTypeArray_2 = meshDataClippedPtr->GetUnstructCells().GetCellTypesArray();
    for(FSIntArrayT::ConstIterator cellType = cellTypeArray_2.BeginConst(); cellType.IsValid(); cellType.Next()) {
      if(*cellType == FSMeshEnums::CellType::CT_Node) {
        meshDataClippedPtr->GetUnstructCells().InitGlobalCellNumber((FSMeshEnums::CellType)*cellType, currentOffset);
        currentOffset += meshDataClippedPtr->GetUnstructCells().GetNCells((FSMeshEnums::CellType)*cellType);
      }
    }
  }
  // Poly2D global numbers come from the clipper's match numbering
  if(!topologyData.poly2DGlobalNumbers.IsEmpty())
    unstructMeshData.InitCellAttribute(globalNumberName, FSMeshEnums::CT_Poly2D, topologyData.poly2DGlobalNumbers);

  CopyCellAttributes(meshDataOriginal, topologyData, unstructMeshData);

  success = meshClipped.IsInitialized();
  if(!success)
    FSError.SetAndPrintAndExit("FSMeshReconstruction : Error while initialization the new mesh");

  return meshClipped;
}

void FSMeshReconstruction::RemapToDistributedNumbering(FSTopologyData& topo, std::vector<FS_intT>& ownedLocalRows,
                                                       FSIntArrayT& ownedGlobalNumbers)
{
  const FS_intT nProcs = clac_.GetNProcs();
  const FS_intT nLocal = static_cast<FS_intT>(topo.globalCoords.size());

  // 1 - Exchange every proc's node list as coordinates quantized by tol_
  std::vector<FS_int32T> counts(nProcs);
  FS_int32T nLocal32 = static_cast<FS_int32T>(nLocal);
  clac_.AllGather(&nLocal32, 1, counts.data(), 1);

  FS_int32T maxCount = 1;
  for(const auto c : counts)
    maxCount = std::max(maxCount, c);

  std::vector<FS_int32T> sendKeys(static_cast<std::size_t>(maxCount) * 3, 0);
  for(FS_intT i = 0; i < nLocal; ++i) {
    const NodeKey k(topo.globalCoords[i], tol_);
    sendKeys[3 * i] = k.ix;
    sendKeys[3 * i + 1] = k.iy;
    sendKeys[3 * i + 2] = k.iz;
  }
  std::vector<FS_int32T> allKeys(static_cast<std::size_t>(maxCount) * 3 * nProcs);
  clac_.AllGather(sendKeys.data(), 3 * maxCount, allKeys.data(), 3 * maxCount);

  // 2 - Numbering in scan order yields contiguous per-proc ranges, as FSDM requires
  std::unordered_map<NodeKey, FS_intT> keyToDistributed;
  keyToDistributed.reserve(static_cast<std::size_t>(maxCount) * nProcs);
  const FS_intT myRank = clac_.GetProcID();

  ownedLocalRows.clear();
  for(FS_intT p = 0; p < nProcs; ++p) {
    for(FS_int32T k = 0; k < counts[p]; ++k) {
      const std::size_t base = (static_cast<std::size_t>(p) * maxCount + k) * 3;
      NodeKey key(FSVec3(0, 0, 0), 1.0);
      key.ix = allKeys[base];
      key.iy = allKeys[base + 1];
      key.iz = allKeys[base + 2];
      const auto [it, inserted] = keyToDistributed.try_emplace(key, static_cast<FS_intT>(keyToDistributed.size()));
      if(inserted && p == myRank)
        ownedLocalRows.push_back(k); // local row k is owned by this proc
    }
  }

  // 3 - Distributing numbers consecutive within each proc's range by construction
  ownedGlobalNumbers.Resize(static_cast<FS_intT>(ownedLocalRows.size()));
  for(std::size_t i = 0; i < ownedLocalRows.size(); ++i)
    ownedGlobalNumbers[static_cast<FS_intT>(i)] =
      keyToDistributed.at(NodeKey(topo.globalCoords[ownedLocalRows[i]], tol_));

  // 4 - Remap the connectivity in place: local row index → distributed number
  auto localToDistributed = [&](FS_intT localRow) {
    return keyToDistributed.at(NodeKey(topo.globalCoords[localRow], tol_));
  };
  auto remapRegister = [&](FSIntRegisterT& reg) {
    for(FS_intT i = 0; i < reg.Size(); ++i)
      for(FS_intT j = 0; j < reg.GetNItems(i); ++j)
        reg(i, j) = localToDistributed(reg(i, j));
  };
  remapRegister(topo.cell2NodePoly2D);
  remapRegister(topo.cell2NodePoly3D);

  for(auto& [type, cell2Node] : topo.cell2NodeInner) {
    for(FS_intT c = 0; c < cell2Node.Size(0); ++c)
      for(FS_intT n = 0; n < cell2Node.Size(1); ++n)
        cell2Node(c, n) = localToDistributed(cell2Node(c, n));
  }
}

void FSMeshReconstruction::CopyCellAttributes(const FSUnstructMeshData& meshDataOriginal, const FSTopologyData& topo,
                                              FSUnstructMeshData& meshDataNew)
{
  for(const auto& [cellType, parentArray] : topo.cellParent) {
    const FS_intT numCells = parentArray.Size();

    if(numCells == 0)
      continue;

    // Resolve parent cell type
    FSMeshEnums::CellType parentCellType = cellType;

    if(FSMeshEnums::IsPolyCellType(cellType)) {
      const auto& parentTypes = topo.cellParentType.at(cellType);

      parentCellType = FSMeshEnums::Int2CellType(parentTypes[0]);

#ifndef NDEBUG
      for(FS_intT i = 1; i < parentTypes.Size(); ++i) {
        if(parentTypes[i] != parentTypes[0]) {
          FSError.SetAndPrintAndExit("CopyCellAttributes: heterogeneous parent types detected");
        }
      }
#endif
    }

    const auto& attribNames = meshDataOriginal.GetCellAttributes(parentCellType);

    const FS_intT offset = meshDataOriginal.GetCellPool(parentCellType)->GetOffset();

    // A Poly2D born from a match on a volume/surface pair split by the partitioner has a parent owned by another proc
    static const FSString attrCADGroupIDName = FSEnums::AttributeTypeToString(FSEnums::AT_CADGroupID);
    const bool isPoly2D = (FSMeshEnums::Int2CellType(cellType) == FSMeshEnums::CellType::CT_Poly2D);
    const bool hasCarriedMarkers = isPoly2D && topo.poly2DMarkers.Size() == numCells;

    for(FSStringArrayT::ConstIterator AI = attribNames.BeginConst(); AI.IsValid(); ++AI) {
      const auto& valuesOrig = meshDataOriginal.GetCellAttribute(*AI, parentCellType);

      FSIntArrayT valuesNew(numCells);

      if(hasCarriedMarkers && *AI == attrCADGroupIDName) {
        for(FS_intT i = 0; i < numCells; ++i)
          valuesNew[i] = topo.poly2DMarkers[i];
      } else {
        for(FS_intT i = 0; i < numCells; ++i) {
          const FS_intT localId = parentArray[i] - offset;
          if(localId < 0 || localId >= valuesOrig.Size()) {
            // Parent owned by another proc and no carried value for this attribute.
            valuesNew[i] = -1;
            continue;
          }
          valuesNew[i] = valuesOrig[localId];
        }
      }

      meshDataNew.InitCellAttribute(*AI, cellType, valuesNew);

      if(meshDataOriginal.HasCellAttributeValueNames(*AI)) {
        meshDataNew.SetCellAttributeValueNames(*AI, meshDataOriginal.GetCellAttributeValueNames(*AI));
      }
    }
  }
}

void FSMeshReconstruction::CopyAttributes(const FSMesh& meshOriginal, FSMesh& meshClipped)
{
  const auto& attributes = meshOriginal.GetAttributeNames();
  for(FSStringArrayT::ConstIterator AI = attributes.BeginConst(); AI.IsValid(); ++AI) {
    const auto attributeName = meshOriginal.GetAttribute(*AI);
    meshClipped.SetAttribute(*AI, attributeName);
  }
}

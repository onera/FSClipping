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

  // Level 2 parallelism: several procs share this mesh. FSDM's global numbering
  // requires each node to be owned by exactly one proc and cell2node to carry
  // the distributed global numbers — assemble that from the clipper IDs.
  const bool distributed = !topologyData.nodeGlobalNumbers.IsEmpty() && clac_.GetNProcs() > 1;
  std::vector<FS_intT> ownedLocalRows;
  FSIntArrayT ownedGlobalNumbers;
  if(distributed) {
    RemapToDistributedNumbering(topologyData, ownedLocalRows, ownedGlobalNumbers);
  } else {
    ownedLocalRows.resize(globalCoords.size());
    for(std::size_t i = 0; i < globalCoords.size(); ++i)
      ownedLocalRows[i] = static_cast<FS_intT>(i);
  }
  const FS_intT nOwned = static_cast<FS_intT>(ownedLocalRows.size());

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

  if(!topologyData.nodeGlobalNumbers.IsEmpty()) {
    // Level 2 parallelism: the clipper-assigned IDs are the stable GlobalNumber
    // of each node — identical on every proc that sees the same interface
    // point, so FSDM can identify shared entities across the sub-communicator.
    const FSString globalNumberName = FSMeshEnums::AttributeTypeToString(FSMeshEnums::AT_GlobalNumber);
    if(distributed) {
      unstructMeshData.InitCellAttribute(globalNumberName, FSMeshEnums::CT_Node, ownedGlobalNumbers);
    } else {
      unstructMeshData.InitCellAttribute(globalNumberName, FSMeshEnums::CT_Node,
                                         topologyData.nodeGlobalNumbers);
    }
    if(!topologyData.poly2DGlobalNumbers.IsEmpty())
      unstructMeshData.InitCellAttribute(globalNumberName, FSMeshEnums::CT_Poly2D,
                                         topologyData.poly2DGlobalNumbers);
  } else {
    // Sequential mode: local numbering is globally valid.
    FS_intT currentOffset = 0;
    FSIntArrayT cellTypeArray_2 = meshDataClippedPtr->GetUnstructCells().GetCellTypesArray();
    for(FSIntArrayT::ConstIterator cellType = cellTypeArray_2.BeginConst(); cellType.IsValid(); cellType.Next()) {
      if(*cellType == FSMeshEnums::CellType::CT_Node) {
        meshDataClippedPtr->GetUnstructCells().InitGlobalCellNumber((FSMeshEnums::CellType)*cellType, currentOffset);
        currentOffset += meshDataClippedPtr->GetUnstructCells().GetNCells((FSMeshEnums::CellType)*cellType);
      }
    }
  }

  CopyCellAttributes(meshDataOriginal, topologyData, unstructMeshData);

  success = meshClipped.IsInitialized();
  if(!success)
    FSError.SetAndPrintAndExit("FSMeshReconstruction : Error while initialization the new mesh");

  return meshClipped;
}

void FSMeshReconstruction::RemapToDistributedNumbering(FSTopologyData& topo,
                                                       std::vector<FS_intT>& ownedLocalRows,
                                                       FSIntArrayT& ownedGlobalNumbers)
{
  const FS_intT nProcs = clac_.GetNProcs();
  const FS_intT nLocal = topo.nodeGlobalNumbers.Size();

  // 1. Exchange the clipper-ID list of every proc (in local row order).
  //    Fixed-size AllGather with padding since counts differ across procs.
  std::vector<FS_int32T> counts(nProcs);
  FS_int32T nLocal32 = static_cast<FS_int32T>(nLocal);
  clac_.AllGather(&nLocal32, 1, counts.data(), 1);

  FS_int32T maxCount = 1;
  for(const auto c : counts)
    maxCount = std::max(maxCount, c);

  std::vector<FS_int32T> sendIds(maxCount, -1);
  for(FS_intT i = 0; i < nLocal; ++i)
    sendIds[i] = static_cast<FS_int32T>(topo.nodeGlobalNumbers[i]);
  std::vector<FS_int32T> allIds(static_cast<std::size_t>(maxCount) * nProcs);
  clac_.AllGather(sendIds.data(), maxCount, allIds.data(), maxCount);

  // 2. Deterministic ownership and contiguous numbering: scanning procs in
  //    ascending rank, the first proc listing a clipper ID owns it. Numbering
  //    in scan order yields contiguous per-proc ranges, as FSDM requires.
  std::unordered_map<FS_intT, FS_intT> clipperToDistributed;
  clipperToDistributed.reserve(static_cast<std::size_t>(maxCount) * nProcs);
  const FS_intT myRank = clac_.GetProcID();

  ownedLocalRows.clear();
  for(FS_intT p = 0; p < nProcs; ++p) {
    for(FS_int32T k = 0; k < counts[p]; ++k) {
      const FS_intT id = allIds[static_cast<std::size_t>(p) * maxCount + k];
      const auto [it, inserted] = clipperToDistributed.try_emplace(id, static_cast<FS_intT>(clipperToDistributed.size()));
      if(inserted && p == myRank)
        ownedLocalRows.push_back(k); // local row k is owned by this proc
    }
  }

  // 3. Stable GlobalNumber attribute of the owned nodes (= their clipper IDs).
  ownedGlobalNumbers.Resize(static_cast<FS_intT>(ownedLocalRows.size()));
  for(std::size_t i = 0; i < ownedLocalRows.size(); ++i)
    ownedGlobalNumbers[static_cast<FS_intT>(i)] = topo.nodeGlobalNumbers[ownedLocalRows[i]];

  // 4. Remap the connectivity in place: local row index → distributed number.
  auto remapRegister = [&](FSIntRegisterT& reg) {
    for(FS_intT i = 0; i < reg.Size(); ++i)
      for(FS_intT j = 0; j < reg.GetNItems(i); ++j)
        reg(i, j) = clipperToDistributed.at(topo.nodeGlobalNumbers[reg(i, j)]);
  };
  remapRegister(topo.cell2NodePoly2D);
  remapRegister(topo.cell2NodePoly3D);
}

void FSMeshReconstruction::CopyCellAttributes(const FSUnstructMeshData& meshDataOriginal, const FSTopologyData& topo,
                                              FSUnstructMeshData& meshDataNew)
{
  for(const auto& [cellType, parentArray] : topo.cellParent) {
    const FS_intT numCells = parentArray.Size();

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

    for(FSStringArrayT::ConstIterator AI = attribNames.BeginConst(); AI.IsValid(); ++AI) {
      const auto& valuesOrig = meshDataOriginal.GetCellAttribute(*AI, parentCellType);

      FSIntArrayT valuesNew(numCells);

      for(FS_intT i = 0; i < numCells; ++i) {
        const FS_intT parentId = parentArray[i];
        valuesNew[i] = valuesOrig[parentId - offset];
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

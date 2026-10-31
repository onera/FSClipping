#include "FSClipping/FSMeshReconstruction.h"
#include "FSClipping/FSNodeKey.h"
#include <FSMeshData.h>

#include <cassert>
#include <map>
#include <set>
#include <sstream>
#include <unordered_map>

namespace {

// A Poly cell may have parents of different types: resolve the parent type cell by cell
FSMeshEnums::CellType ParentTypeOf(const FSTopologyData& topo, FSMeshEnums::CellType outType, FS_intT i)
{
  if(!FSMeshEnums::IsPolyCellType(outType))
    return outType;
  return FSMeshEnums::Int2CellType(topo.cellParentType.at(outType)[i]);
}

// FSDM rejects an empty values array: an empty pool gets its attribute from the fill overload
void InitIntCellAttribute(FSUnstructMeshData& meshData, const FSString& name, FSMeshEnums::CellType t,
                          FSIntArrayT& values)
{
  if(values.Size() == 0)
    meshData.InitCellAttribute(name, t, 1, -1);
  else
    meshData.InitCellAttribute(name, t, values);
}

} // namespace

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

void FSMeshReconstruction::UniteAttributeNames(std::map<FSMeshEnums::CellType, std::set<FSString>>& namesOf)
{
  const FS_intT nProcs = clac_.GetNProcs();
  if(nProcs == 1)
    return;

  // One "type name" line per attribute, names hold no whitespace
  std::string local;
  for(const auto& [type, names] : namesOf)
    for(const FSString& name : names)
      local += std::to_string(static_cast<FS_intT>(type)) + ' ' + name + '\n';

  FS_int32T localSize = static_cast<FS_int32T>(local.size());
  std::vector<FS_int32T> sizes(nProcs);
  clac_.AllGather(&localSize, 1, sizes.data(), 1);

  for(FS_intT p = 0; p < nProcs; ++p) {
    if(sizes[p] == 0)
      continue;
    std::string buf = (p == clac_.GetProcID()) ? local : std::string(sizes[p], '\0');
    clac_.Broadcast(buf.data(), sizes[p], p);

    std::istringstream in(buf);
    FS_intT type;
    std::string name;
    while(in >> type >> name)
      namesOf[FSMeshEnums::Int2CellType(type)].insert(FSString(name));
  }
}

void FSMeshReconstruction::CopyCellAttributes(const FSUnstructMeshData& meshDataOriginal, const FSTopologyData& topo,
                                              FSUnstructMeshData& meshDataNew)
{
  // Pool offset of each parent type present, per output type
  std::map<FSMeshEnums::CellType, std::map<FSMeshEnums::CellType, FS_intT>> offsetsOf;
  // Union of the attribute names of those parent types, per output type
  std::map<FSMeshEnums::CellType, std::set<FSString>> attribNamesOf;

  for(const auto& [cellType, parentArray] : topo.cellParent) {
    const FSMeshEnums::CellType outType = FSMeshEnums::CellType(cellType);
    auto& offsetOf = offsetsOf[outType];
    auto& attribNames = attribNamesOf[outType];

    for(FS_intT i = 0; i < parentArray.Size(); ++i) {
      const FSMeshEnums::CellType t = ParentTypeOf(topo, outType, i);
      if(offsetOf.contains(t))
        continue;
      // Undefined parent type, or no pool of that type on this proc: nothing to read
      const FSCellPool* pool = meshDataOriginal.GetCellPool(t);
      if(!pool)
        continue;
      offsetOf.emplace(t, pool->GetOffset());
      const FSStringArrayT names = meshDataOriginal.GetCellAttributes(t);
      for(FSStringArrayT::ConstIterator AI = names.BeginConst(); AI.IsValid(); ++AI)
        attribNames.insert(*AI);
    }
  }

  // Same attributes on every proc, including those with no cell of a type
  UniteAttributeNames(attribNamesOf);

  static const FSIntArrayT noParent;
  for(auto& [outType, attribNames] : attribNamesOf) {
    // The type has no pool in the new mesh on this proc: nothing to carry the attribute
    if(!meshDataNew.GetCellPool(outType))
      continue;

    const auto parentIt = topo.cellParent.find(outType);
    const FSIntArrayT& parentArray = (parentIt != topo.cellParent.end()) ? parentIt->second : noParent;
    const FS_intT numCells = parentArray.Size();
    const auto& offsetOf = offsetsOf[outType];

    // Never inherit the global number: several clipped cells may share one parent
    static const FSString attrGlobalNumberName = FSMeshEnums::AttributeTypeToString(FSMeshEnums::AT_GlobalNumber);
    attribNames.erase(attrGlobalNumberName);

    // A Poly2D born from a match on a volume/surface pair split by the partitioner has a parent owned by another proc
    static const FSString attrCADGroupIDName = FSEnums::AttributeTypeToString(FSEnums::AT_CADGroupID);
    const bool isPoly2D = (outType == FSMeshEnums::CellType::CT_Poly2D);
    const bool hasCarriedMarkers = isPoly2D && topo.poly2DMarkers.Size() == numCells;

    for(const FSString& name : attribNames) {
      FSIntArrayT valuesNew(numCells);
      const bool useCarriedMarkers = hasCarriedMarkers && name == attrCADGroupIDName;

      // nullptr when this parent type does not carry the attribute
      std::map<FSMeshEnums::CellType, const FSIntArrayT*> valuesOf;
      for(const auto& [t, offset] : offsetOf)
        valuesOf[t] =
          meshDataOriginal.HasCellAttribute(name, t) ? &meshDataOriginal.GetCellAttribute(name, t) : nullptr;

      for(FS_intT i = 0; i < numCells; ++i) {
        // A carried marker wins; -1 means none was carried, so read the parent
        if(useCarriedMarkers && topo.poly2DMarkers[i] >= 0) {
          valuesNew[i] = topo.poly2DMarkers[i];
          continue;
        }
        const auto it = valuesOf.find(ParentTypeOf(topo, outType, i));
        const FSIntArrayT* valuesOrig = (it != valuesOf.end()) ? it->second : nullptr;
        const FS_intT localId = valuesOrig ? parentArray[i] - offsetOf.at(it->first) : -1;
        // Parent owned by another proc, without local pool, or attribute missing on its type: no value to copy
        if(!valuesOrig || localId < 0 || localId >= valuesOrig->Size()) {
          valuesNew[i] = -1;
          continue;
        }
        valuesNew[i] = (*valuesOrig)[localId];
      }

      InitIntCellAttribute(meshDataNew, name, outType, valuesNew);

      if(meshDataOriginal.HasCellAttributeValueNames(name)) {
        meshDataNew.SetCellAttributeValueNames(name, meshDataOriginal.GetCellAttributeValueNames(name));
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

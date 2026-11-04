#include "FSClipping/FSMeshReconstruction.h"
#include <FSCellInfo.h>
#include <FSMeshData.h>

#include <cassert>
#include <map>
#include <set>
#include <unordered_map>

namespace {

// A Poly cell may have parents of different types: resolve the parent type cell by cell
FSMeshEnums::CellType ParentTypeOf(const FSTopologyData& topo, FSMeshEnums::CellType outType, FS_intT i)
{
  if(!FSMeshEnums::IsPolyCellType(outType))
    return outType;
  return FSMeshEnums::Int2CellType(topo.cellParentType.at(outType)[i]);
}

// Pool offset of each parent type of `outType` cells that has a pool on this proc
std::map<FSMeshEnums::CellType, FS_intT> ParentOffsets(const FSUnstructMeshData& meshData, const FSTopologyData& topo,
                                                       FSMeshEnums::CellType outType)
{
  std::map<FSMeshEnums::CellType, FS_intT> offsetOf;
  const auto parentIt = topo.cellParent.find(outType);
  if(parentIt == topo.cellParent.end())
    return offsetOf;
  for(FS_intT i = 0; i < parentIt->second.Size(); ++i) {
    const FSMeshEnums::CellType t = ParentTypeOf(topo, outType, i);
    if(offsetOf.contains(t))
      continue;
    // Undefined parent type, or no pool of that type on this proc: nothing to read
    if(const FSCellPool* pool = meshData.GetCellPool(t))
      offsetOf.emplace(t, pool->GetOffset());
  }
  return offsetOf;
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

  // Same pools and attributes on every proc, fixed before any collective init
  const FSPoolLayout layout = PoolLayout(meshDataOriginal);

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
  // InitUnstructCells is collective: loop over the shared layout, never over the local types
  for(const auto& [type, names] : layout) {
    if(FSMeshEnums::IsPolyCellType(type))
      continue;
    FSIntArrayT& cell2Node = cell2NodeInner[type];
    // FSDM wants a 0 x nNodes connectivity for a pool empty on this proc
    if(!cell2Node.IsInitialized())
      cell2Node.Resize(0, FSNNodes(type));
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

  // Node global numbers only if the original mesh has them
  const FSString globalNumberName = FSMeshEnums::AttributeTypeToString(FSMeshEnums::AT_GlobalNumber);
  if(meshDataOriginal.HasCellAttribute(globalNumberName, FSMeshEnums::CT_Node)) {
    if(distributed)
      unstructMeshData.InitCellAttribute(globalNumberName, FSMeshEnums::CT_Node, ownedGlobalNumbers);
    else
      unstructMeshData.InitGlobalCellNumber(FSMeshEnums::CT_Node, 0);
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

  // 1 - Exchange every proc's node identities: exact, so no tolerance and one integer per node
  std::vector<FS_int32T> counts(nProcs);
  FS_int32T nLocal32 = static_cast<FS_int32T>(nLocal);
  clac_.AllGather(&nLocal32, 1, counts.data(), 1);

  FS_int32T maxCount = 1;
  for(const auto c : counts)
    maxCount = std::max(maxCount, c);

  const FSIntArrayT& ids = topo.nodeGlobalNumbers;
  assert(ids.Size() == nLocal);
  std::vector<FS_intT> sendIds(maxCount, -1);
  for(FS_intT i = 0; i < nLocal; ++i)
    sendIds[i] = ids[i];
  std::vector<FS_intT> allIds(static_cast<std::size_t>(maxCount) * nProcs);
  clac_.AllGather(sendIds.data(), maxCount, allIds.data(), maxCount);

  // 2 - The first proc that holds a node owns it; scan order yields contiguous per-proc ranges, as FSDM requires
  std::unordered_map<FS_intT, FS_intT> idToDistributed;
  idToDistributed.reserve(static_cast<std::size_t>(maxCount) * nProcs);
  const FS_intT myRank = clac_.GetProcID();

  ownedLocalRows.clear();
  for(FS_intT p = 0; p < nProcs; ++p) {
    for(FS_int32T k = 0; k < counts[p]; ++k) {
      const FS_intT id = allIds[static_cast<std::size_t>(p) * maxCount + k];
      const auto [it, inserted] = idToDistributed.try_emplace(id, static_cast<FS_intT>(idToDistributed.size()));
      if(inserted && p == myRank)
        ownedLocalRows.push_back(k); // local row k is owned by this proc
    }
  }

  // 3 - Distributing numbers consecutive within each proc's range by construction
  ownedGlobalNumbers.Resize(static_cast<FS_intT>(ownedLocalRows.size()));
  for(std::size_t i = 0; i < ownedLocalRows.size(); ++i)
    ownedGlobalNumbers[static_cast<FS_intT>(i)] = idToDistributed.at(ids[ownedLocalRows[i]]);

  // 4 - Remap the connectivity in place: local row index → distributed number
  auto localToDistributed = [&](FS_intT localRow) { return idToDistributed.at(ids[localRow]); };
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

FSPoolLayout FSMeshReconstruction::PoolLayout(const FSUnstructMeshData& meshData)
{
  static const FSString attrGlobalNumberName = FSMeshEnums::AttributeTypeToString(FSMeshEnums::AT_GlobalNumber);
  FSPoolLayout layout;
  // Poly pools always exist: a clipped cell may come from any proc
  auto& poly2DNames = layout[FSMeshEnums::CT_Poly2D];
  auto& poly3DNames = layout[FSMeshEnums::CT_Poly3D];
  for(const FSMeshEnums::CellType t : meshData.GetCellTypes()) {
    if(t == FSMeshEnums::CT_Node)
      continue;
    auto& names = layout[t];
    const FSStringArrayT parentNames = meshData.GetCellAttributes(t);
    for(FSStringArrayT::ConstIterator AI = parentNames.BeginConst(); AI.IsValid(); ++AI) {
      // Never inherit the global number: several clipped cells may share one parent
      if(*AI == attrGlobalNumberName)
        continue;
      names.insert(*AI);
      // A clipped cell inherits from a parent of the same dimension
      if(FSMeshEnums::IsUnstructSurfaceCellType(t))
        poly2DNames.insert(*AI);
      else if(FSMeshEnums::IsUnstructVolumeCellType(t))
        poly3DNames.insert(*AI);
    }
  }
  return layout;
}

void FSMeshReconstruction::CopyCellAttributes(const FSUnstructMeshData& meshDataOriginal, const FSTopologyData& topo,
                                              FSUnstructMeshData& meshDataNew)
{
  static const FSIntArrayT noParent;
  for(const auto& [outType, attribNames] : PoolLayout(meshDataOriginal)) {
    // The type has no pool in the new mesh: nothing to carry the attribute
    if(!meshDataNew.GetCellPool(outType))
      continue;

    const auto parentIt = topo.cellParent.find(outType);
    const FSIntArrayT& parentArray = (parentIt != topo.cellParent.end()) ? parentIt->second : noParent;
    const FS_intT numCells = parentArray.Size();
    const auto offsetOf = ParentOffsets(meshDataOriginal, topo, outType);

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

#include "FSClipping/FSMeshReconstruction.h"
#include <FSMeshData.h>

FSMesh FSMeshReconstruction::Build(const FSUnstructMeshData& meshDataOriginal, FSTopologyData& topologyData)
{
  FSMesh meshClipped = FSMesh(&clac_);

  const auto& globalCoords = topologyData.globalCoords;
  auto& cell2NodePoly2D = topologyData.cell2NodePoly2D;
  auto& cell2NodePoly3D = topologyData.cell2NodePoly3D;
  auto& cell2NodeInner = topologyData.cell2NodeInner;
  auto& polyFaces = topologyData.polyFaces;

  meshClipped.BeginInitialization();
  FSMeshData* meshDataClippedPtr = meshClipped.GetMeshData();
  FSUnstructMeshData& unstructMeshData = meshDataClippedPtr->GetUnstructCells();

  meshClipped.InitUnstructNodes(globalCoords.size());
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

  FSFloatArrayT coordInterface(globalCoords.size(), 3);
  for(std::size_t i = 0; i < globalCoords.size(); i++) {
    coordInterface(i, 0) = globalCoords[i][0];
    coordInterface(i, 1) = globalCoords[i][1];
    coordInterface(i, 2) = globalCoords[i][2];
  }

  bool success = unstructMeshData.SetCoordinates3D(coordsDesc, coordInterface);
  if(!success)
    FSError.SetAndPrintAndExit("FSMeshReconstruction : Error while set coordinates");

  FS_intT currentOffset = 0;
  FSIntArrayT cellTypeArray_2 = meshDataClippedPtr->GetUnstructCells().GetCellTypesArray();
  for(FSIntArrayT::ConstIterator cellType = cellTypeArray_2.BeginConst(); cellType.IsValid(); cellType.Next()) {
    if(*cellType == FSMeshEnums::CellType::CT_Node) {
      meshDataClippedPtr->GetUnstructCells().InitGlobalCellNumber((FSMeshEnums::CellType)*cellType, currentOffset);
      currentOffset += meshDataClippedPtr->GetUnstructCells().GetNCells((FSMeshEnums::CellType)*cellType);
    }
  }

  CopyCellAttributes(meshDataOriginal, topologyData, unstructMeshData);

  success = meshClipped.IsInitialized();
  if(!success)
    FSError.SetAndPrintAndExit("FSMeshReconstruction : Error while initialization the new mesh");

  return meshClipped;
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

#include "FSClipping/FSMeshReconstruction.h"
#include <FSMeshData.h>

FSMeshData FSMeshReconstruction::Build(const FSUnstructMeshData& meshDataOriginal, FSTopologyData& topologyData)
{

  FSMeshData meshDataInterface = FSMeshData(&clac_);
  // init pointer to meshdata
  // FSMeshData* meshDataInterfacePtr = &meshDataInterface;
  // init reference to unstructured mesh data
  FSUnstructMeshData& unstructMeshData = meshDataInterface.GetUnstructCells();

  const auto& globalCoords = topologyData.globalCoords;
  auto& cell2NodePoly2D = topologyData.cell2NodePoly2D;
  auto& cell2NodePoly3D = topologyData.cell2NodePoly3D;
  auto& cell2NodeInner = topologyData.cell2NodeInner;
  auto& polyFaces = topologyData.polyFaces;

  // begin initialization
  meshDataInterface.BeginInitialization();

  meshDataInterface.InitUnstructNodes(globalCoords.size());
  meshDataInterface.InitUnstructCells(FSMeshEnums::CT_Poly3D, cell2NodePoly3D);
  meshDataInterface.InitUnstructCells(FSMeshEnums::CT_Poly2D, cell2NodePoly2D);

  if(!cell2NodeInner.IsEmpty()) {
    for(auto& [type, cell2Node] : topologyData.cell2NodeInner)
      meshDataInterface.InitUnstructCells(type, cell2Node);
  }
  meshDataInterface.InitUnstructCellFaces(FSMeshEnums::CT_Poly3D, polyFaces);
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

  CopyCellAttributes(meshDataOriginal, topologyData, unstructMeshData);
  //  end initialization
  meshDataInterface.EndInitialization();
  // check if initialization is complete
  success = meshDataInterface.IsInitialized();

  if(!success)
    FSError.SetAndPrintAndExit("FSMeshReconstruction : Error while initialization the new mesh");

  return meshDataInterface;
}


// void FSMeshReconstruction::CopyCellAttributes(const FSUnstructMeshData& meshDataOriginal,
//                                               const FSTopologyData& topo,
//                                               FSUnstructMeshData& meshDataNew)
//{
//   for(const auto& [cellType, parentArray] : topo.cellParent) {
//     const FS_intT numCells = parentArray.Size();
//     FS_intT parentCellType = cellType;
//     if(FSMeshEnums::IsPolyCellType(cellType)) {
//       auto intArrayT = topo.cellParentType.at(cellType);
//       parentCellType = intArrayT[0];
//     }
//     const auto& attribNames = meshDataOriginal.GetCellAttributes(FSMeshEnums::Int2CellType(parentCellType));
//     FS_intT offSet = meshDataOriginal.GetCellPool(FSMeshEnums::Int2CellType(parentCellType))->GetOffset();
//
//     for(FSStringArrayT::ConstIterator AI = attribNames.BeginConst(); AI.IsValid(); ++AI) {
//       const auto& valuesOrig = meshDataOriginal.GetCellAttribute(*AI, FSMeshEnums::Int2CellType(parentCellType));
//       FSIntArrayT valuesNew(numCells);
//
//       for(FS_intT i = 0; i < numCells; ++i) {
//         FS_intT parent = parentArray[i] - offSet;
//         valuesNew[i] = valuesOrig[parent];
//       }
//
//       meshDataNew.InitCellAttribute(*AI, cellType, valuesNew);
//
//       if(meshDataOriginal.HasCellAttributeValueNames(*AI)) {
//         meshDataNew.SetCellAttributeValueNames(
//           *AI,
//           meshDataOriginal.GetCellAttributeValueNames(*AI));
//       }
//     }
//   }
// }

void FSMeshReconstruction::CopyCellAttributes(const FSUnstructMeshData& meshDataOriginal,
                                              const FSTopologyData& topo,
                                              FSUnstructMeshData& meshDataNew)
{
  for(const auto& [cellType, parentArray] : topo.cellParent) {
    const FS_intT numCells = parentArray.Size();

    // -----------------------------
    // Resolve parent cell type
    // -----------------------------
    FSMeshEnums::CellType parentCellType = cellType;

    if(FSMeshEnums::IsPolyCellType(cellType)) {
      const auto& parentTypes = topo.cellParentType.at(cellType);

      parentCellType = FSMeshEnums::Int2CellType(parentTypes[0]);

#ifndef NDEBUG
      for(FS_intT i = 1; i < parentTypes.Size(); ++i) {
        if(parentTypes[i] != parentTypes[0]) {
          FSError.SetAndPrintAndExit(
            "CopyCellAttributes: heterogeneous parent types detected");
        }
      }
#endif
    }

    // -----------------------------
    // Access original data
    // -----------------------------
    const auto& attribNames =
      meshDataOriginal.GetCellAttributes(parentCellType);

    const FS_intT offset =
      meshDataOriginal.GetCellPool(parentCellType)->GetOffset();

    // -----------------------------
    // Copy attributes
    // -----------------------------
    for(FSStringArrayT::ConstIterator AI = attribNames.BeginConst(); AI.IsValid(); ++AI) {
      const auto& valuesOrig =
        meshDataOriginal.GetCellAttribute(*AI, parentCellType);

      FSIntArrayT valuesNew(numCells);

      for(FS_intT i = 0; i < numCells; ++i) {
        const FS_intT parentId = parentArray[i];
        valuesNew[i] = valuesOrig[parentId - offset];
      }

      meshDataNew.InitCellAttribute(*AI, cellType, valuesNew);

      if(meshDataOriginal.HasCellAttributeValueNames(*AI)) {
        meshDataNew.SetCellAttributeValueNames(
          *AI,
          meshDataOriginal.GetCellAttributeValueNames(*AI));
      }
    }
  }
}
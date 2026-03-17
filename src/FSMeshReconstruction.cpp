#include "FSClipping/FSMeshReconstruction.h"
#include <FSMeshData.h>

FSMeshData FSMeshReconstruction::Build(FSTopologyData& topologyData)
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

  FS_intT currentOffset = 0;
  FSIntArrayT cellTypeArray_2 = meshDataInterface.GetUnstructCells().GetCellTypesArray();
  for(FSIntArrayT::ConstIterator cellType = cellTypeArray_2.BeginConst(); cellType.IsValid(); cellType.Next()) {
    // meshDataInterface.GetUnstructCells().InitGlobalCellNumber((FSMeshEnums::CellType)*cellType, currentOffset);
    currentOffset += meshDataInterface.GetUnstructCells().GetNCells((FSMeshEnums::CellType)*cellType);
  }

  // end initialization
  meshDataInterface.EndInitialization();
  // check if initialization is complete
  success = meshDataInterface.IsInitialized();

  if(!success)
    FSError.SetAndPrintAndExit("FSMeshReconstruction : Error while initialization the new mesh");

  return meshDataInterface;
}
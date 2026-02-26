#include "FSClipping/FSClippingInterface.h"
#include "FSClipping/FSBoundaryFaceProvider.h"
#include "FSClipping/FSFaceMatcher.h"
#include "FSClipping/FSMeshReconstruction.h"
#include "FSClipping/FSTopologyBuilder.h"
#include <FSMeshData.h>

_FS_BEGIN_NAMESPACE

FSMeshData FSClippingInterface::BuildInterface(FSMesh& mesh1, FSMesh& mesh2, FSClac& clacInterface)
{

  FSMeshFaceExtractor ex1, ex2;

  std::vector<FSClippingFace> clipFaces1 = FSBoundaryFaceProvider::Extract(mesh1, ex1, boundaryMarkerMesh1_);
  std::vector<FSClippingFace> clipFaces2 = FSBoundaryFaceProvider::Extract(mesh2, ex2, boundaryMarkerMesh2_);

  /// 2- Compute the matches and the corresponding intersections
  FSFaceMatcher matcher(clac2_, clipFaces1, clipFaces2, tol_);
  std::vector<FSFaceMatch> matches;
  matcher.ComputeMatches(matches);

  /// 3- Reconstruct the connectivity cell2Node and face2Node of the new border
  FSTopologyBuilder topologyBuilder(tol_);
  auto topologySurfacique = topologyBuilder.BuildSurfaceTopo(matches);
  // auto topologyVolumique = topologyBuilder(matches, oldCell2Node, oldCell2Face2Node, oldCoord, CellPool)
  //  auto globalCoords = topology.globalCoords;
  //  auto cell2Node = topology.cell2Node;
  //  auto polyFaces = topology.polyFaces;

  /// 4- Reconstruct the interface between the two mesh border
  FSMeshReconstruction FSMeshReconstruction(clacInterface);
  auto meshDataInterface = FSMeshReconstruction.Build(topologySurfacique);

  return meshDataInterface;
}


// void collectCellType2Cell2Node(FSMesh& mesh,
//                                const std::unordered_map<FS_intT, std::set<FS_intT> >& bdryElemType,
//                                FSCellType2IntArrayT& cellType2Cell2NodeInner,
//                                FSCellType2IntArrayT& cellType2Cell2NodeBdry)
//{
//   const auto cellTypes = mesh.GetCellTypes();
//
//   for(const auto& t : cellTypes) {
//     const FSMeshEnums::CellType type = t.Value();
//
//     if(!FSMeshEnums::IsUnstructVolumeCellType(type))
//       continue;
//
//     const auto& cell2NodeOld = mesh.GetCell2Node(t);
//
//     const FS_intT numNodes = FSCellInfo::NNodes(t);
//     const FS_intT size = cell2NodeOld.Size();
//     const FS_intT offset = cell2NodeOld.Offset();
//     const FS_intT numCells = size / numNodes;
//
//     std::vector<FS_intT> innerData;
//     std::vector<FS_intT> bdryData;
//
//     innerData.reserve(size);
//     bdryData.reserve(size);
//
//     const auto it = bdryElemType.find(type);
//     const std::set<FS_intT>* bdrySet = nullptr;
//
//     if(it != bdryElemType.end())
//       bdrySet = &it->second;
//
//     for(FS_intT i = 0; i < numCells; ++i) {
//       const FS_intT cellIndex = offset + i;
//       const bool isBdry = (bdrySet && bdrySet->contains(cellIndex));
//
//       auto& target = isBdry ? bdryData : innerData;
//
//       for(FS_intT n = 0; n < numNodes; ++n)
//         target.emplace_back(cell2NodeOld[i * numNodes + n]);
//     }
//
//     // --- Build FSIntArrayT INNER ---
//     if(!innerData.empty()) {
//       FSIntArrayT innerArray(
//         innerData.size() / numNodes,
//         numNodes);
//
//       innerArray.CopyData(innerData.data(), innerData.size());
//       cellType2Cell2NodeInner.Insert(type, innerArray);
//     }
//
//     // --- Build FSIntArrayT BDRY ---
//     if(!bdryData.empty()) {
//       FSIntArrayT bdryArray(
//         bdryData.size() / numNodes,
//         numNodes);
//       // bdryArray.SetOffset()
//       bdryArray.CopyData(bdryData.data(), bdryData.size());
//       cellType2Cell2NodeBdry.Insert(type, bdryArray);
//     }
//   }
// }

FSMeshData FSClippingInterface::BuildInterfaceTmp(FSMesh& mesh1, FSMesh& mesh2, FSClac& clacInterface)
{

  /// 1- Boundary extraction
  FSMeshFaceExtractor ex1, ex2;

  auto clipFaces1 = FSBoundaryFaceProvider::Extract(mesh1, ex1, boundaryMarkerMesh1_);
  auto clipFaces2 = FSBoundaryFaceProvider::Extract(mesh2, ex2, boundaryMarkerMesh2_);

  auto bdryCells1 = FSBoundaryFaceProvider::ExtractOwnerCells(clipFaces1);
  auto bdryCells2 = FSBoundaryFaceProvider::ExtractOwnerCells(clipFaces2);

  /// 2- Compute the matches and the corresponding intersections
  FSFaceMatcher matcher(clac2_, clipFaces1, clipFaces2, tol_);
  std::vector<FSFaceMatch> matches;
  matcher.ComputeMatches(matches);

  /// 3- Split volume connectivity and get coordinates
  // FSCellType2IntArrayT cellType2Cell2NodeInner;
  // FSCellType2IntArrayT cellType2Cell2NodeBdry;
  // collectCellType2Cell2Node(mesh1, bdryCells1, cellType2Cell2NodeInner, cellType2Cell2NodeBdry);

  FSFloatArrayT oldCoords;
  FS_intT nodeOffset;
  FSQuantityDescArrayT coordDesc;
  mesh1.GetMeshData()->GetUnstructCells().GetCoordinates3D(coordDesc, oldCoords, nodeOffset);

  /// 4- Compute the new topology of the border, first the surface one with the matches and then the volume
  FSTopologyBuilder topologyBuilder(tol_);
  FSTopologyData surfaceTopology = topologyBuilder.BuildSurfaceTopo(matches);

  FSTopologyData volumeTopology;
  const auto& cellType = mesh1.GetCellTypes();
  for(const auto& t : cellType) {
    if(FSMeshEnums::IsUnstructVolumeCellType(t)) {
      const auto& cell2Node1 = mesh1.GetCell2Node(t);
      const auto& cellPool1 = mesh1.GetMeshData()->GetUnstructCells().GetCellPool(t);

      volumeTopology = topologyBuilder.BuildVolumeTopo(cell2Node1, bdryCells1.at(t), *cellPool1, oldCoords);
    }
  }

  auto globalCoords = volumeTopology.globalCoords;
  auto cell2Node = volumeTopology.cell2Node;
  auto polyFaces = volumeTopology.polyFaces;

  // 4- Reconstruct the interface between the two mesh border
  FSMeshReconstruction FSMeshReconstruction(clacInterface);
  auto meshDataInterface = FSMeshReconstruction.Build(volumeTopology);

  return meshDataInterface;
}


_FS_END_NAMESPACE

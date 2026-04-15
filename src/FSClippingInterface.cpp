#include "FSClipping/FSClippingInterface.h"
#include "FSClipping/FSBoundaryFaceProvider.h"
#include "FSClipping/FSFaceMatcher.h"
#include "FSClipping/FSMeshReconstruction.h"
#include "FSClipping/FSTopologyAssembler.h"
#include <FSMeshData.h>

_FS_BEGIN_NAMESPACE

FSMeshData FSClippingInterface::BuildSurfaceInterface(FSMesh& mesh1, FSMesh& mesh2)
{
  // 1. Extract boundary faces
  FSMeshFaceExtractor extractor1, extractor2;

  auto boundaryExtraction1 = FSBoundaryFaceProvider::Extract(mesh1, extractor1, boundaryMarkerMesh1_, tol_);
  auto boundaryExtraction2 = FSBoundaryFaceProvider::Extract(mesh2, extractor2, boundaryMarkerMesh2_, tol_);

  auto clipFaces1 = boundaryExtraction1.faces;
  auto clipFaces2 = boundaryExtraction2.faces;

  // 2. Compute geometric matches between boundary faces
  FSFaceMatcher matcher(clac2_, clipFaces1, clipFaces2, tol_);
  std::vector<FSFaceMatch> matches;
  matcher.ComputeMatches(matches);

  // 3. Build surface topology from matches
  FSTopologyAssembler topologyBuilder(tol_);
  FSTopologyData surfaceTopology =
    topologyBuilder.BuildSurfaceTopo(matches, boundaryExtraction1.faceKeys);

  // 4. Reconstruct surface mesh
  FSMeshReconstruction meshReconstruction(clac1_);
  return meshReconstruction.Build(mesh1.GetMeshData()->GetUnstructCells(), surfaceTopology);
}

FSMeshData FSClippingInterface::BuildVolumeInterface(FSMesh& mesh1, FSMesh& mesh2)
{
  /// 1- Boundary extraction
  FSMeshFaceExtractor ex1, ex2;
  FSTopologyData volumeTopology;

  auto boundaryExtraction1 = FSBoundaryFaceProvider::Extract(mesh1, ex1, boundaryMarkerMesh1_, tol_);
  auto boundaryExtraction2 = FSBoundaryFaceProvider::Extract(mesh2, ex2, boundaryMarkerMesh2_, tol_);

  auto clipFaces1 = boundaryExtraction1.faces;
  auto clipFaces2 = boundaryExtraction2.faces;

  std::unordered_map<FS_intT, std::set<FS_intT> > bdry3DCellsType1 = FSBoundaryFaceProvider::ExtractOwner3DCells(clipFaces1);
  std::unordered_map<FS_intT, std::set<FS_intT> > bdry3DCellsType2 = FSBoundaryFaceProvider::ExtractOwner3DCells(clipFaces2);

  std::unordered_map<FS_intT, std::set<FS_intT> > bdry2DCellsType1 = FSBoundaryFaceProvider::Extract2DCells(clipFaces1);
  std::unordered_map<FS_intT, std::set<FS_intT> > bdry2DCellsType2 = FSBoundaryFaceProvider::Extract2DCells(clipFaces2);

  /// 2- Compute the matches and the corresponding intersections
  FSFaceMatcher matcher(clac2_, clipFaces1, clipFaces2, tol_);
  std::vector<FSFaceMatch> matches;
  matcher.ComputeMatches(matches);

  FSFloatArrayT oldCoords;
  FS_intT nodeOffset;
  FSQuantityDescArrayT coordDesc;
  mesh1.GetMeshData()->GetUnstructCells().GetCoordinates3D(coordDesc, oldCoords, nodeOffset);

  /// 3- Compute the new surface topology of the border
  FSTopologyAssembler topologyBuilder(tol_);
  FSTopologyData surfaceTopology = topologyBuilder.BuildSurfaceTopo(matches, boundaryExtraction1.faceKeys);

  /// 4- Then compute the topology of all the meshes
  const auto& cellType = mesh1.GetCellTypes();
  for(const auto& t : cellType) {
    if(FSMeshEnums::IsUnstructVolumeCellType(t)) {
      const auto& cell2Node1 = mesh1.GetCell2Node(t);
      const auto& cellPool1 = mesh1.GetMeshData()->GetUnstructCells().GetCellPool(t);
      const auto& bdryCellPool1 = bdry3DCellsType1.at(t);

      topologyBuilder.BuildVolumeTopo(cell2Node1, bdryCellPool1, *cellPool1, oldCoords, volumeTopology);
      // volumeTopology = topolyBuilder.BuildVolumeTopo(cell2Node, surfaceTopology); <- idealement
    }
  }

  topologyBuilder.AppendUnclippedSurfaces(mesh1, bdry2DCellsType1, oldCoords, volumeTopology);

  volumeTopology.cellParent[FSMeshEnums::CellType::CT_Poly2D] = surfaceTopology.cellParent[FSMeshEnums::CellType::CT_Poly2D];
  volumeTopology.cellParentType[FSMeshEnums::CellType::CT_Poly2D] = surfaceTopology.cellParentType[FSMeshEnums::CellType::CT_Poly2D];

  // 5- Reconstruct the entire new mesh
  FSMeshReconstruction FSMeshReconstruction(clac1_);
  auto meshDataInterface = FSMeshReconstruction.Build(mesh1.GetMeshData()->GetUnstructCells(), volumeTopology);

  return meshDataInterface;
}


_FS_END_NAMESPACE

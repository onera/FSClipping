#include "FSClipping/FSClippingInterface.h"
#include "FSClipping/FSBoundaryFaceProvider.h"
#include "FSClipping/FSFaceMatcher.h"
#include "FSClipping/FSMeshReconstruction.h"
#include "FSClipping/FSTopologyAssembler.h"
#include <FSMeshData.h>

_FS_BEGIN_NAMESPACE

FSMeshData FSClippingInterface::BuildSurfaceInterface(FSMesh& mesh1, FSMesh& mesh2, FSClac& clacInterface)
{
  // 1. Extract boundary faces
  FSMeshFaceExtractor extractor1, extractor2;

  auto clipFaces1 = FSBoundaryFaceProvider::Extract(
    mesh1, extractor1, boundaryMarkerMesh1_);

  auto clipFaces2 = FSBoundaryFaceProvider::Extract(
    mesh2, extractor2, boundaryMarkerMesh2_);

  // 2. Compute geometric matches between boundary faces
  FSFaceMatcher matcher(clac2_, clipFaces1, clipFaces2, tol_);
  std::vector<FSFaceMatch> matches;
  matcher.ComputeMatches(matches);

  // 3. Build surface topology from matches
  FSTopologyAssembler topologyBuilder(tol_);
  FSTopologyData surfaceTopology =
    topologyBuilder.BuildSurfaceTopo(matches);

  // 4. Reconstruct surface mesh
  FSMeshReconstruction meshReconstruction(clacInterface);
  return meshReconstruction.Build(surfaceTopology);
}

FSMeshData FSClippingInterface::BuildVolumeInterface(FSMesh& mesh1, FSMesh& mesh2, FSClac& clacInterface)
{

  /// 1- Boundary extraction
  FSMeshFaceExtractor ex1, ex2;

  auto clipFaces1 = FSBoundaryFaceProvider::Extract(mesh1, ex1, boundaryMarkerMesh1_);
  auto clipFaces2 = FSBoundaryFaceProvider::Extract(mesh2, ex2, boundaryMarkerMesh2_);

  std::unordered_map<FS_intT, std::set<FS_intT> > bdryCells1 = FSBoundaryFaceProvider::ExtractOwnerCells(clipFaces1);
  std::unordered_map<FS_intT, std::set<FS_intT> > bdryCells2 = FSBoundaryFaceProvider::ExtractOwnerCells(clipFaces2);

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
  FSTopologyData surfaceTopology = topologyBuilder.BuildSurfaceTopo(matches);

  /// 4- Then compute the topology of all the meshes
  FSTopologyData volumeTopology;
  const auto& cellType = mesh1.GetCellTypes();
  for(const auto& t : cellType) {
    if(FSMeshEnums::IsUnstructVolumeCellType(t)) {
      const auto& cell2Node1 = mesh1.GetCell2Node(t);
      const auto& cellPool1 = mesh1.GetMeshData()->GetUnstructCells().GetCellPool(t);
      const auto& bdryCellPool1 = bdryCells1.at(t);

      volumeTopology = topologyBuilder.BuildVolumeTopo(cell2Node1, bdryCellPool1, *cellPool1, oldCoords);
      // volumeTopology = topolyBuilder.BuildVolumeTopo(surfaceTopology);
    }
  }

  // 5- Reconstruct the entire new mesh
  FSMeshReconstruction FSMeshReconstruction(clacInterface);
  auto meshDataInterface = FSMeshReconstruction.Build(volumeTopology);

  return meshDataInterface;
}


_FS_END_NAMESPACE

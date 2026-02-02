#include "FSClipping/FSClippingInterface.h"
#include "FSClipping/FSBoundaryFaceProvider.h"
#include "FSClipping/FSFaceMatcher.h"
#include "FSClipping/FSMeshReconstruction.h"
#include "FSClipping/FSTopologyBuilder.h"
#include <FSMeshData.h>

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
  auto topology = topologyBuilder.Build(matches);

  auto globalCoords = topology.globalCoords;
  auto cell2Node = topology.cell2Node;
  auto polyFaces = topology.polyFaces;

  /// 4- Reconstruct the interface between the two mesh border
  FSMeshReconstruction FSMeshReconstruction(clacInterface);
  auto meshDataInterface = FSMeshReconstruction.Build(topology);

  return meshDataInterface;
}
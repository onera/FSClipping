#include "FSClipping/FSClippingInterface.h"
#include "FSClipping/FSBoundaryFaceProvider.h"
#include "FSClipping/FSClippingEngine.h"
#include "FSClipping/FSFaceMatcher.h"

_FS_BEGIN_NAMESPACE

FSMesh FSClippingInterface::BuildInterface(FSMesh& mesh1, FSMesh& mesh2)
{
  // 1. Extract boundary faces
  FSMeshFaceExtractor extractor1, extractor2;

  auto boundaryExtraction1 = FSBoundaryFaceProvider::Extract(mesh1, extractor1, boundaryMarkerMesh1_, tol_);
  auto boundaryExtraction2 = FSBoundaryFaceProvider::Extract(mesh2, extractor2, boundaryMarkerMesh2_, tol_);

  // 2. Compute geometric matches between boundary faces
  FSFaceMatcher matcher(clac2_, boundaryExtraction1.faces, boundaryExtraction2.faces, tol_);
  std::vector<FSFaceMatch> matches;
  matcher.ComputeMatches(matches);
  if(matches.empty()) {
    FSError.SetAndExit("FSClippingInterface: No matches were found");
  }
  FSFaceMatcher::AssignGlobalNodeIds(matches, tol_);

  // 3. Build topology and reconstruct the mesh
  FSClippingEngine engine(tol_);
  FSTopologyData topology = engine.BuildTopology(mesh1, boundaryExtraction1, matches);
  return engine.Reconstruct(clac1_, mesh1, topology);
}

_FS_END_NAMESPACE

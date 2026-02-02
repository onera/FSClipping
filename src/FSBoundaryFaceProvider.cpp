#include "FSClipping/FSBoundaryFaceProvider.h"
#include "FSClipping/FSFaceSeparator.h"

_FS_BEGIN_NAMESPACE


std::vector<FSClippingFace>
FSBoundaryFaceProvider::Extract(FSMesh& mesh,
                                FSMeshFaceExtractor& ex,
                                FS_intT boundaryMarker)
{

  // mesh.GetMeshData()->GetUnstructCells().CreateLocalNumbering();
  ex.PrepareFaceConnectivity(mesh.GetMeshData()->GetUnstructCells(), true,
                             false);
  ex.PrepareFaceNodeCoordinates(FSQuantityDescArrayT());


  auto bdryFaces =
    FSFaceSeparator::SeparateBoundariesFaceWithMarker(mesh, ex, boundaryMarker);

  std::vector<FSClippingFace> clipFaces;
  // clipFaces.reserve(bdryFaces.size());

  FSFloatArrayT faceCoordinates;
  for(const auto& f : bdryFaces) {
    FS_intT faceIndex = ex.GetFaceIndex(*(f._faceFSDM));
    ex.GetFaceNodeCoordinates(faceIndex, faceCoordinates);
    clipFaces.emplace_back(f, faceIndex, faceCoordinates);
  }

  return clipFaces;
}
_FS_BEGIN_NAMESPACE
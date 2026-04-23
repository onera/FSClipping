#include "FSClipping/FSBoundaryFaceProvider.h"
#include "FSClipping/FSFaceSeparator.h"

_FS_BEGIN_NAMESPACE

BoundaryExtraction
FSBoundaryFaceProvider::Extract(FSMesh& mesh,
                                FSMeshFaceExtractor& ex,
                                FS_intT boundaryMarker,
                                FS_floatT tol)
{

  // mesh.GetMeshData()->GetUnstructCells().CreateLocalNumbering();
  assert(mesh.GetMeshData()->GetUnstructCells().HasLocalNumbering());
  BoundaryExtraction result;

  ex.PrepareFaceConnectivity(mesh.GetMeshData()->GetUnstructCells(), true,
                             false);
  ex.PrepareFaceNodeCoordinates(FSQuantityDescArrayT());

  // The FSDM domain local bdry faces (owning element (FSDM volume cell) also
  // local)
  std::vector<FSBoundaryFace> bdryFaces;
  // Locally used, but external (i.e. located somewhere remote) FSDM bdry faces.
  // std::vector<FSFace> extBdryFacesFSDM;
  // Separated local FSDM bdry faces, i.e. owner element somewhere remote.
  // std::vector<FSBoundaryFace> sepBdryFacesFSDM;
  // the FSDM domain local faces (between two local elements, respectively)
  // std::vector<FSFace> localFacesFSDM;
  // the FSDM halo-touching faces, i.e. owning element is local, but neighboring
  // element somewhere remote
  // std::vector<FSFace> haloFacesFSDM;

  FSFaceSeparator::SeparateFaces(mesh, ex, bdryFaces, result, boundaryMarker);

  std::vector<FSClippingFace> clipFaces;
  // clipFaces.reserve(bdryFaces.size());


  FSFloatArrayT faceCoordinates;
  for(const auto& f : bdryFaces) {
    FS_intT faceIndex = ex.GetFaceIndex(*(f._faceFSDM));
    ex.GetFaceNodeCoordinates(faceIndex, faceCoordinates);

    GeomFaceKey key(faceCoordinates, tol);

    result.faceKeys.insert(key);
    result.faces.emplace_back(f, faceIndex, faceCoordinates);

    // volume (owner)
    result.volumeCells[f._faceFSDM->mOwner.mCellType]
      .insert(f._faceFSDM->mOwner.mCell);

    // surface (neighbor)
    result.surfaceCells[f._faceFSDM->mNeighbor.mCellType]
      .insert(f._faceFSDM->mNeighbor.mCell);
  }

  return result;
}
_FS_BEGIN_NAMESPACE
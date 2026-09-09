#include "FSClipping/FSBoundaryFaceProvider.h"
#include "FSClipping/FSFaceSeparator.h"

_FS_BEGIN_NAMESPACE

BoundaryExtraction FSBoundaryFaceProvider::Extract(FSMesh& mesh, FSMeshFaceExtractor& ex, FS_intT boundaryMarker,
                                                   FS_floatT tol)
{
  const bool hasLocalNumbering = mesh.GetMeshData()->GetUnstructCells().HasLocalNumbering();
  if(hasLocalNumbering)
    mesh.GetMeshData()->GetUnstructCells().CreateLocalNumbering();

  BoundaryExtraction result;

  const FS_intT procID = mesh.GetClac()->GetProcID();

  ex.PrepareFaceConnectivity(mesh.GetMeshData()->GetUnstructCells());

  ex.PrepareFaceNodeCoordinates(FSQuantityDescArrayT());

  std::vector<FSBoundaryFace> bdryFaces;
  FSFaceSeparator::SeparateFaces(mesh, ex, bdryFaces, boundaryMarker);

  FSFloatArrayT faceCoordinates;
  for(const auto& f : bdryFaces) {
    FS_intT faceIndex = ex.GetFaceIndex(*(f._faceFSDM));
    ex.GetFaceNodeCoordinates(faceIndex, faceCoordinates);

    GeomFaceKey key(faceCoordinates, tol);
    // On est vraiment sur qu'on ne peut pas utiliser les tables de hachages de FSDM ?

    result.faceKeys.insert(key);
    result.faces.emplace_back(f, faceIndex, faceCoordinates);

    result.volumeCells[f._faceFSDM->mOwner.mCellType].insert(f._faceFSDM->mOwner.mCell);
    // Only a local surface cell may be listed here: the set is used to drop the clipped cells from this proc's
    // surface pools, and a remote cell index would collide with an unrelated local one.
    if(f._faceFSDM->mNeighbor.mCellProcID == procID)
      result.surfaceCells[f._faceFSDM->mNeighbor.mCellType].insert(f._faceFSDM->mNeighbor.mCell);
  }

  return result;
}
_FS_END_NAMESPACE

#include "FSClipping/FSBoundaryFaceProvider.h"
#include "FSClipping/FSFaceSeparator.h"

_FS_BEGIN_NAMESPACE

BoundaryExtraction FSBoundaryFaceProvider::Extract(FSMesh& mesh, FSMeshFaceExtractor& ex, FS_intT boundaryMarker)
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

    result.faces.emplace_back(f, faceIndex, faceCoordinates);

    auto [it, inserted] = result.cell2BoundaryFaces.try_emplace(f._faceFSDM->mOwner.mCell);
    it->second.push_back(f._faceFSDM->mOwner.mFaceIndex);

    // Only a local surface cell may be listed here
    if(f._faceFSDM->mNeighbor.mCellProcID == procID)
      result.surfaceCells.insert(f._faceFSDM->mNeighbor.mCell);
  }

  return result;
}
_FS_END_NAMESPACE

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


std::unordered_map<FS_intT, std::set<FS_intT> >
FSBoundaryFaceProvider::ExtractOwnerCells(const std::vector<FSClippingFace>& faces)
{
  std::unordered_map<FS_intT, std::set<FS_intT> > bdryElemTypes;
  for(const auto& f : faces) {
    auto elemIndex = f.topo()._faceFSDM->mOwner.mCell;
    auto elemType = f.topo()._faceFSDM->mOwner.mCellType;
    if(bdryElemTypes.contains(elemType)) {
      bdryElemTypes.at(elemType).insert(elemIndex);
    } else {
      std::set<FS_intT> index = {elemIndex};
      bdryElemTypes.insert({elemType, index});
    }
  }
  return bdryElemTypes;
}
_FS_BEGIN_NAMESPACE
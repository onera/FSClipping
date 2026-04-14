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
  ex.PrepareFaceConnectivity(mesh.GetMeshData()->GetUnstructCells(), true,
                             false);
  ex.PrepareFaceNodeCoordinates(FSQuantityDescArrayT());

  auto bdryFaces =
    FSFaceSeparator::SeparateBoundariesFaceWithMarker(mesh, ex, boundaryMarker);

  std::vector<FSClippingFace> clipFaces;
  // clipFaces.reserve(bdryFaces.size());

  BoundaryExtraction result;
  FSFloatArrayT faceCoordinates;
  for(const auto& f : bdryFaces) {
    FS_intT faceIndex = ex.GetFaceIndex(*(f._faceFSDM));
    ex.GetFaceNodeCoordinates(faceIndex, faceCoordinates);

    GeomFaceKey key(faceCoordinates, tol);

    result.faceKeys.insert(key);
    result.faces.emplace_back(f, faceIndex, faceCoordinates);
  }

  return result;
}


std::unordered_map<FS_intT, std::set<FS_intT> >
FSBoundaryFaceProvider::ExtractOwner3DCells(const std::vector<FSClippingFace>& faces)
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

std::unordered_map<FS_intT, std::set<FS_intT> >
FSBoundaryFaceProvider::Extract2DCells(const std::vector<FSClippingFace>& faces)
{
  std::unordered_map<FS_intT, std::set<FS_intT> > bdryCellTypes;
  for(const auto& f : faces) {
    auto cellIndex = f.topo()._faceFSDM->mNeighbor.mCell;
    auto cellType = f.topo()._faceFSDM->mNeighbor.mCellType;
    if(bdryCellTypes.contains(cellType)) {
      bdryCellTypes.at(cellType).insert(cellIndex);
    } else {
      std::set<FS_intT> index = {cellIndex};
      bdryCellTypes.insert({cellType, index});
    }
  }
  return bdryCellTypes;
}
_FS_BEGIN_NAMESPACE
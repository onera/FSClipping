#include "FSClipping/FSFaceSeparator.h"


_FS_BEGIN_NAMESPACE

namespace FSFaceSeparator {

void SeparateFaces(const FSMesh& fsmesh, const FSMeshFaceExtractor& faceExtractor,
                   std::vector<FSBoundaryFace>& bdryFaces, const FS_intT markerBoundaryClipped)
{
  assert(bdryFaces.empty());

  FS_intT procID = fsmesh.GetClac()->GetProcID();
  const FS_intT numAllFaces = faceExtractor.GetNLocalFaces();
  assert(fsmesh.GetMeshData());
  const FSUnstructMeshData& meshData = fsmesh.GetMeshData()->GetUnstructCells();

  const FSString attrCADGroupIDName = FSEnums::AttributeTypeToString(FSEnums::AT_CADGroupID);
  const FSString attrGlobalNumberName = FSMeshEnums::AttributeTypeToString(FSMeshEnums::AT_GlobalNumber);

  for(FS_intT i = 0; i < numAllFaces; ++i) {
    const auto& face = faceExtractor.GetFaceConnectivity(i);
    assert(face.mOwner.mCellProcID == procID);
    assert(face.mOwner.mCell >= 0);
    const FSCellPool* ownerCellPool = meshData.GetCellPool(face.mOwner.mCellType);
    assert(ownerCellPool);
    assert(ownerCellPool->IsOwned(face.mOwner.mCell));

    if(face.mOwner.mFaceIndex >= 0) // owner cell is a volume cell -- and
                                    // faceIdx is the face's index in this cell
    {
      if(face.mNeighbor.mCellProcID == procID) // neighbor cell is local, too
      {
        if(FSMeshEnums::IsUnstructSurfaceCellType(face.mNeighbor.mCellType)) {
          const FSCellPool* boundaryCellPool = meshData.GetCellPool(face.mNeighbor.mCellType);

          FS_intT marker = -1;
          if(boundaryCellPool->HasCellAttribute(attrCADGroupIDName)) {
            const FSIntArrayT& boundaryMarkers = boundaryCellPool->GetCellAttribute(attrCADGroupIDName);
            marker = boundaryMarkers(face.mNeighbor.mCell);
          }

          FS_intT fsdmUniqueBdryFaceID = -1;
          if(boundaryCellPool->HasCellAttribute(attrGlobalNumberName)) {
            const FSIntArrayT& globalBoundaryFacesIDs = boundaryCellPool->GetCellAttribute(attrGlobalNumberName);
            fsdmUniqueBdryFaceID = globalBoundaryFacesIDs(face.mNeighbor.mCell);
          }
          if(marker == markerBoundaryClipped)
            bdryFaces.emplace_back(face, marker, fsdmUniqueBdryFaceID);
        }
      }
    }
  }

  // order the faces according to FSDM cell number of the owner and neighbor
  // cell, since their original order in the face extractor is implementation
  // defined (but currently depends on the hash of the node number(s)
  // referenced)
  std::sort(bdryFaces.begin(), bdryFaces.end());
}

} // namespace FSFaceSeparator
_FS_END_NAMESPACE

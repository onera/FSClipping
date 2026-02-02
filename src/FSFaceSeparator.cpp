#include "FSClipping/FSFaceSeparator.h"

_FS_BEGIN_NAMESPACE

namespace FSFaceSeparator {

void SeparateFaces(const FSMesh& fsmesh,
                   const FSMeshFaceExtractor& faceExtractor,
                   std::vector<FSFace>& localFacesFSDM,
                   std::vector<FSFace>& haloFacesFSDM,
                   std::vector<FSBoundaryFace>& bdryFacesFSDM,
                   std::vector<FSFace>& extBdryFacesFSDM,
                   std::vector<FSBoundaryFace>& sepBdryFacesFSDM)
{
  assert(localFacesFSDM.empty() && haloFacesFSDM.empty() &&
         bdryFacesFSDM.empty() && extBdryFacesFSDM.empty() &&
         sepBdryFacesFSDM.empty());

  FS_intT procID = fsmesh.GetClac()->GetProcID();
  const FS_intT numAllFaces = faceExtractor.GetNLocalFaces();
  assert(fsmesh.GetMeshData());
  const FSUnstructMeshData& meshData = fsmesh.GetMeshData()->GetUnstructCells();

  std::vector<FS_intT> localCellIDsFSDM;
  std::vector<FSCellAddress> remoteCellIDsFSDM;

  localCellIDsFSDM.reserve(2 * numAllFaces);

  const FSString attrCADGroupIDName =
    FSEnums::AttributeTypeToString(FSEnums::AT_CADGroupID);
  const FSString attrGlobalNumberName =
    FSMeshEnums::AttributeTypeToString(FSMeshEnums::AT_GlobalNumber);

  FS_intT faceCounter = 0;
  for(FS_intT i = 0; i < numAllFaces; ++i) {
    const auto& face = faceExtractor.GetFaceConnectivity(i);
    assert(face.mOwner.mCellProcID == procID);
    assert(face.mOwner.mCell >= 0);
    const FSCellPool* ownerCellPool =
      meshData.GetCellPool(face.mOwner.mCellType);
    // IgnoreUnused(ownerCellPool);
    assert(ownerCellPool);
    assert(ownerCellPool->IsOwned(face.mOwner.mCell));

    if(face.mNeighbor.mCell < 0)
      throw std::runtime_error(
        "CODA only supports FSMesh instances as input for which the "
        "FSMeshFaceExtractor delivers a valid neighbor "
        "cell for each face. In particular, the complete boundary of the "
        "computational domain must be represented by "
        "surface cells.");

    if(FSMeshEnums::IsUnstructSurfaceCellType(face.mOwner.mCellType) &&
       face.mNeighbor.mCellProcID == procID)
      throw std::runtime_error(
        "CODA only supports FSMesh instances as input for which the "
        "FSMeshFaceExtractor does not deliver an "
        "(unstructured) surface cell as the owner of a face, unless the "
        "neighbor cell is remote. In particular, "
        "CODA does not support pseudo cells mapping a 'large' surface cell "
        "to multiple 'small' boundary faces.");

    ++faceCounter;
    if(face.mOwner.mFaceIndex >= 0) // owner cell is a volume cell -- and
                                    // faceIdx is the face's index in this cell
    {
      localCellIDsFSDM.push_back(face.mOwner.mCell);
      if(face.mNeighbor.mCellProcID == procID) // neighbor cell is local, too
      {
        if(FSMeshEnums::IsUnstructSurfaceCellType(face.mNeighbor.mCellType)) {
          const FSCellPool* boundaryCellPool =
            meshData.GetCellPool(face.mNeighbor.mCellType);

          FS_intT marker = -1;
          if(boundaryCellPool->HasCellAttribute(attrCADGroupIDName)) {
            const FSIntArrayT& boundaryMarkers =
              boundaryCellPool->GetCellAttribute(attrCADGroupIDName);
            marker = boundaryMarkers(face.mNeighbor.mCell);
          }

          FS_intT fsdmUniqueBdryFaceID = -1;
          if(boundaryCellPool->HasCellAttribute(attrGlobalNumberName)) {
            const FSIntArrayT& globalBoundaryFacesIDs =
              boundaryCellPool->GetCellAttribute(attrGlobalNumberName);
            fsdmUniqueBdryFaceID = globalBoundaryFacesIDs(face.mNeighbor.mCell);
          }

          bdryFacesFSDM.emplace_back(face, marker, fsdmUniqueBdryFaceID);
        } else // i.e. neighbor no unstructured surface cell -> should be a
               // volume cell...
        {
          assert(
            FSMeshEnums::IsUnstructVolumeCellType(face.mNeighbor.mCellType));
          localCellIDsFSDM.push_back(face.mNeighbor.mCell);
          localFacesFSDM.emplace_back(face);
        }
      } else // i.e. neighborCell located on a different domain
        if(FSMeshEnums::IsUnstructSurfaceCellType(face.mNeighbor.mCellType)) {
          extBdryFacesFSDM.emplace_back(face);
        } else // remote neighbor is a volume cell
        {
          assert(
            FSMeshEnums::IsUnstructVolumeCellType(face.mNeighbor.mCellType));
          remoteCellIDsFSDM.emplace_back(face.mNeighbor.mCellProcID,
                                         face.mNeighbor.mCell);
          haloFacesFSDM.emplace_back(face);
        }
    } else // i.e. faceIDx < 0, indicating that ownerCell is not a volume cell
    {
      assert(FSMeshEnums::IsUnstructSurfaceCellType(face.mOwner.mCellType));
      const FSCellPool* boundaryCellPool =
        meshData.GetCellPool(face.mOwner.mCellType);
      assert(boundaryCellPool);

      FS_intT marker = -1;
      if(boundaryCellPool->HasCellAttribute(attrCADGroupIDName)) {
        const FSIntArrayT& boundaryMarkers =
          boundaryCellPool->GetCellAttribute(attrCADGroupIDName);
        marker = boundaryMarkers(face.mOwner.mCell);
      }
      assert(marker != -1 &&
             "assuming that there is no marker '-1' in the mesh, is there?");

      FS_intT fsdmUniqueBdryFaceID = -1;
      if(boundaryCellPool->HasCellAttribute(attrGlobalNumberName)) {
        const FSIntArrayT& globalBoundaryFacesIDs =
          boundaryCellPool->GetCellAttribute(attrGlobalNumberName);
        fsdmUniqueBdryFaceID = globalBoundaryFacesIDs(face.mOwner.mCell);
      }
      sepBdryFacesFSDM.emplace_back(face, marker, fsdmUniqueBdryFaceID);
    }
  }

  // order the faces according to FSDM cell number of the owner and neighbor
  // cell, since their original order in the face extractor is implementation
  // defined (but currently depends on the hash of the node number(s)
  // referenced)
  std::sort(localFacesFSDM.begin(), localFacesFSDM.end());
  std::sort(haloFacesFSDM.begin(), haloFacesFSDM.end());
  std::sort(bdryFacesFSDM.begin(), bdryFacesFSDM.end());
  std::sort(extBdryFacesFSDM.begin(), extBdryFacesFSDM.end());
  std::sort(sepBdryFacesFSDM.begin(), sepBdryFacesFSDM.end());

  // As an FSDM cell ID may have been added multiple times, make each cell ID
  // appear once and only once.
  std::sort(localCellIDsFSDM.begin(), localCellIDsFSDM.end());
  localCellIDsFSDM.erase(
    std::unique(localCellIDsFSDM.begin(), localCellIDsFSDM.end()),
    localCellIDsFSDM.end());
  //_localCellIDsFSDM.assign(localCellIDsFSDM.begin(), localCellIDsFSDM.end());

  std::sort(remoteCellIDsFSDM.begin(), remoteCellIDsFSDM.end());
  remoteCellIDsFSDM.erase(
    std::unique(remoteCellIDsFSDM.begin(), remoteCellIDsFSDM.end()),
    remoteCellIDsFSDM.end());
  //_remoteCellIDsFSDM.assign(remoteCellIDsFSDM.begin(),
  // remoteCellIDsFSDM.end());

  // check there are as many local(ly) separated boundary faces as there are
  // external boundary faces somewhere remote
  auto missingBdryFaces = static_cast<int>(
    extBdryFacesFSDM.size() - sepBdryFacesFSDM.size()); // possibly != 0
  if(fsmesh.GetClac()->NProcs() > 1)
    fsmesh.GetClac()->SumInPlace(&missingBdryFaces, 1);
  assert(missingBdryFaces == 0);
}

bool PrepareNodalData(const FSMesh& fsmesh,
                      FSMeshFaceExtractor& faceExtractor)
{
  FSQuantityDescArrayT nodalDataDesc;
  // we want to extract coordinates and maybe also velocities and *maybe* even
  // derivative information thereof
  if(!fsmesh.DetermineCoordinates(nodalDataDesc)) {
    auto msg = FSError.Get();
    msg += " in FSMesh::DetermineCoordinates() called during "
           "FaceBasedMeshAdapter::PrepareNodalData()";
    throw std::runtime_error(msg.c_str());
  }

  bool okFlag = faceExtractor.PrepareFaceNodalData(nodalDataDesc);

  return okFlag;
}

std::vector<FSBoundaryFace> SeparateBoundariesFaceWithMarker(
  FSMesh& fsmesh, FSMeshFaceExtractor& faceExtractor, const FS_intT marker)
{

  assert(fsmesh.GetMeshData()->GetUnstructCells().HasLocalNumbering());

  // The FSDM domain local bdry faces (owning element (FSDM volume cell) also
  // local)
  std::vector<FSBoundaryFace> bdryFacesFSDM;
  std::vector<FSBoundaryFace> bdryFacesFiltered;

  // Locally used, but external (i.e. located somewhere remote) FSDM bdry faces.
  std::vector<FSFace> extBdryFacesFSDM;
  // Separated local FSDM bdry faces, i.e. owner element somewhere remote.
  std::vector<FSBoundaryFace> sepBdryFacesFSDM;
  // the FSDM domain local faces (between two local elements, respectively)
  std::vector<FSFace> localFacesFSDM;
  // the FSDM halo-touching faces, i.e. owning element is local, but neighboring
  // element somewhere remote
  std::vector<FSFace> haloFacesFSDM;

  SeparateFaces(fsmesh, faceExtractor, localFacesFSDM, haloFacesFSDM,
                bdryFacesFSDM, extBdryFacesFSDM, sepBdryFacesFSDM);
  bdryFacesFiltered.reserve(bdryFacesFSDM.size());

  for(const auto& bdryFaces : bdryFacesFSDM) {
    if(marker == bdryFaces._marker)
      bdryFacesFiltered.push_back(bdryFaces);
  }
  return bdryFacesFiltered;
}
} // namespace FSFaceSeparator
_FS_END_NAMESPACE

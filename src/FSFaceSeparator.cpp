#include "FSClipping/FSFaceSeparator.h"

#include <map>
#include <tuple>
#include <vector>

_FS_BEGIN_NAMESPACE

namespace FSFaceSeparator {

namespace {

/// Attributes of a boundary surface cell needed to decide whether its face is clipped.
struct SurfaceCellAttributes {
  FS_intT marker;
  FS_intT fsdmGlobalNumber;
};

/// Reads the clipping-relevant attributes of a locally owned surface cell.
/**
 * @param meshData The local unstructured mesh data holding the surface cell pools.
 * @param cellType The type of the surface cell.
 * @param cell The index of the surface cell, which must be local to this proc.
 * @return Its FSDM boundary marker and global number, each -1 when the mesh does not carry that attribute.
 */
SurfaceCellAttributes ReadSurfaceCellAttributes(const FSUnstructMeshData& meshData, FSMeshEnums::CellType cellType,
                                                FS_intT cell)
{
  static const FSString attrCADGroupIDName = FSEnums::AttributeTypeToString(FSEnums::AT_CADGroupID);
  static const FSString attrGlobalNumberName = FSMeshEnums::AttributeTypeToString(FSMeshEnums::AT_GlobalNumber);

  SurfaceCellAttributes attributes{-1, -1};

  const FSCellPool* boundaryCellPool = meshData.GetCellPool(cellType);
  assert(boundaryCellPool);
  if(!boundaryCellPool)
    return attributes;

  if(boundaryCellPool->HasCellAttribute(attrCADGroupIDName)) {
    const FSIntArrayT& boundaryMarkers = boundaryCellPool->GetCellAttribute(attrCADGroupIDName);
    attributes.marker = boundaryMarkers(cell);
  }

  if(boundaryCellPool->HasCellAttribute(attrGlobalNumberName)) {
    const FSIntArrayT& globalBoundaryFacesIDs = boundaryCellPool->GetCellAttribute(attrGlobalNumberName);
    attributes.fsdmGlobalNumber = globalBoundaryFacesIDs(cell);
  }

  return attributes;
}

/// Orders a separated/external boundary face by the volume/surface cell pair it belongs to.
std::tuple<FS_intT, FS_intT, FS_intT, FS_intT> VolumeSidePairKey(const FSFaceConnectivity& face, bool ownerIsVolume)
{
  const FSCellFace& volume = ownerIsVolume ? face.mOwner : face.mNeighbor;
  const FSCellFace& surface = ownerIsVolume ? face.mNeighbor : face.mOwner;
  return std::make_tuple(volume.mCellProcID, volume.mCell, surface.mCellProcID, surface.mCell);
}

/// Sends the marker of the separated boundary faces to the procs owning their volume cell.
void ExchangeSeparatedBdryFaces(FSClac& clac, const std::vector<FSBoundaryFace>& sepBdryFaces,
                                const std::vector<FSFace>& extBdryFaces, const FS_intT markerBoundaryClipped,
                                std::vector<FSBoundaryFace>& bdryFaces)
{
  // send the marker of those local bdry faces that are needed somewhere remote
  std::map<FS_intT, FSClac::sizeT> proc2sendSize;
  for(const auto& sbf : sepBdryFaces) {
    assert(sbf._faceFSDM->mNeighbor.mCellProcID != -1);
    proc2sendSize[sbf._faceFSDM->mNeighbor.mCellProcID] += sbf._faceFSDM->GetBufSize(clac) + clac.GetBufSizeInt32(2);
  }

  STLIntSetT destProcs;
  for(const auto& [targetProcID, sendSize] : proc2sendSize) {
    clac.SelectSendBuffer(targetProcID);
    clac.InitSendBuffer(sendSize);
    for(const auto& sbf : sepBdryFaces) {
      if(sbf._faceFSDM->mNeighbor.mCellProcID != targetProcID)
        continue;
      FSFaceConnectivity face = *sbf._faceFSDM; // copy data into send buffer
      face.Pack(clac);
      FS_intT marker = sbf._marker;
      FS_intT fsdmFaceGlobalNumber = sbf._fsdmFaceGlobalNumber;
      clac.Pack(&marker, 1);
      clac.Pack(&fsdmFaceGlobalNumber, 1);
    }
    // now start sending the separated boundary faces
    clac.NBSendBuffer(targetProcID, false);
    destProcs.insert(targetProcID);
  }

  // receive all external boundary faces
  STLIntSetT sourceProcs;
  for(const auto& ebf : extBdryFaces)
    sourceProcs.insert(ebf._faceFSDM->mNeighbor.mCellProcID);

  for(const auto sourceProcID : sourceProcs) {
    clac.SelectReceiveBuffer(sourceProcID);
    clac.NBReceiveBuffer(sourceProcID, false);
  }
  clac.WaitAllReceives(sourceProcs, false);

  FS_intT sourceProcID = -1;
  for(const auto& ebf : extBdryFaces) {
    if(ebf._faceFSDM->mNeighbor.mCellProcID != sourceProcID) {
      // next/new neighbor process found to receive boundary faces from
      sourceProcID = ebf._faceFSDM->mNeighbor.mCellProcID;
      clac.SelectReceiveBuffer(sourceProcID);
    }

    FSFaceConnectivity sbf; // what we receive
    sbf.Unpack(clac);
    FS_intT marker = -1;
    FS_intT fsdmFaceGlobalNumber = -1;
    clac.Unpack(&marker, 1);
    clac.Unpack(&fsdmFaceGlobalNumber, 1);

    assert(sbf.mOwner == ebf._faceFSDM->mNeighbor); // check that the remote owner is the local neighbor
    assert(sbf.mNeighbor == ebf._faceFSDM->mOwner); // and vice versa

    if(marker == markerBoundaryClipped)
      bdryFaces.emplace_back(*ebf._faceFSDM, marker, fsdmFaceGlobalNumber);
  }

  // finally, finalize the send requests
  clac.WaitAllSends(destProcs, false);
}

} // namespace

void SeparateFaces(const FSMesh& fsmesh, const FSMeshFaceExtractor& faceExtractor,
                   std::vector<FSBoundaryFace>& bdryFaces, const FS_intT markerBoundaryClipped)
{
  assert(bdryFaces.empty());

  FS_intT procID = fsmesh.GetClac()->GetProcID();
  const FS_intT numAllFaces = faceExtractor.GetNLocalFaces();
  assert(fsmesh.GetMeshData());
  const FSUnstructMeshData& meshData = fsmesh.GetMeshData()->GetUnstructCells();

  /*
  The two halves of a volume/surface pair split by the partitioner. Neither is usable on its own: the separated
   side holds the marker, the external side holds the face the clipping needs. ExchangeSeparatedBdryFaces() matches
   them, so these lists never leave this function.
  */

  // Locally used, but external (i.e. located somewhere remote) FSDM bdry faces.
  std::vector<FSFace> extBdryFaces;
  // Separated local FSDM bdry faces, i.e. owner element somewhere remote.
  std::vector<FSBoundaryFace> sepBdryFaces;

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
      if(!FSMeshEnums::IsUnstructSurfaceCellType(face.mNeighbor.mCellType))
        continue; // neighbor is a volume cell, so this face is not on the boundary

      if(face.mNeighbor.mCellProcID == procID) // neighbor cell is local, too
      {
        const SurfaceCellAttributes attributes =
          ReadSurfaceCellAttributes(meshData, face.mNeighbor.mCellType, face.mNeighbor.mCell);
        if(attributes.marker == markerBoundaryClipped)
          bdryFaces.emplace_back(face, attributes.marker, attributes.fsdmGlobalNumber);
      } else // i.e. the surface cell closing this volume cell is located on a different domain
        extBdryFaces.emplace_back(face);
    } else // i.e. faceIdx < 0, indicating that ownerCell is not a volume cell
    {
      assert(FSMeshEnums::IsUnstructSurfaceCellType(face.mOwner.mCellType));
      const SurfaceCellAttributes attributes =
        ReadSurfaceCellAttributes(meshData, face.mOwner.mCellType, face.mOwner.mCell);
      sepBdryFaces.emplace_back(face, attributes.marker, attributes.fsdmGlobalNumber);
    }
  }

  // Order both halves on the same key, so that the n-th face a proc sends is the n-th face the other proc expects.
  std::sort(extBdryFaces.begin(), extBdryFaces.end(), [](const FSFace& a, const FSFace& b) {
    return VolumeSidePairKey(*a._faceFSDM, true) < VolumeSidePairKey(*b._faceFSDM, true);
  });
  std::sort(sepBdryFaces.begin(), sepBdryFaces.end(), [](const FSBoundaryFace& a, const FSBoundaryFace& b) {
    return VolumeSidePairKey(*a._faceFSDM, false) < VolumeSidePairKey(*b._faceFSDM, false);
  });

  if(fsmesh.GetClac()->GetNProcs() > 1) {
    // check there are as many local separated boundary faces as there are external boundary faces somewhere remote
    auto missingBdryFaces = static_cast<FS_int32T>(extBdryFaces.size() - sepBdryFaces.size()); // possibly != 0
    fsmesh.GetClac()->SumInPlace(&missingBdryFaces, 1);
    assert(missingBdryFaces == 0);

    ExchangeSeparatedBdryFaces(*fsmesh.GetClac(), sepBdryFaces, extBdryFaces, markerBoundaryClipped, bdryFaces);
  }

  // order the faces according to FSDM cell number of the owner and neighbor cell
  std::sort(bdryFaces.begin(), bdryFaces.end());
}

} // namespace FSFaceSeparator
_FS_END_NAMESPACE

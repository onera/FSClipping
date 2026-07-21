#include "FSClipping/FSClippingEngine.h"
#include "FSClipping/FSFaceExchange.h"
#include "FSClipping/FSMeshReconstruction.h"
#include <FSMeshData.h>

_FS_BEGIN_NAMESPACE

//-----------------------------------------------------------------------------
//
//  BuildTopology
//

FSTopologyData FSClippingEngine::BuildTopology(FSMesh& mesh, const BoundaryExtraction& be,
                                               const std::vector<FSFaceMatch>& matches, Mode mode)
{
  FSTopologyAssembler topologyAssembler(tol_);
  FSTopologyData surfaceTopology = topologyAssembler.BuildSurfaceTopo(matches, be.faceKeys);
  if(mode == Mode::Surface)
    return surfaceTopology;

  FSFloatArrayT oldCoords;
  FS_intT nodeOffset;
  FSQuantityDescArrayT coordDesc;
  mesh.GetMeshData()->GetUnstructCells().GetCoordinates3D(coordDesc, oldCoords, nodeOffset);

  FSTopologyData volumeTopology;
  // A proc whose partition has no boundary face of a given type has no entry
  // in be.volumeCells — BuildVolumeTopo with an empty pool passes all cells
  // of that type through unclipped.
  static const std::set<FS_intT> kEmptyPool;
  for(const auto& t : mesh.GetCellTypes()) {
    if(FSMeshEnums::IsUnstructVolumeCellType(t)) {
      const auto& cell2Node = mesh.GetCell2Node(t);
      const auto& cellPool = mesh.GetMeshData()->GetUnstructCells().GetCellPool(t);
      const auto it = be.volumeCells.find(t);
      const auto& bdryCellPool = (it != be.volumeCells.end()) ? it->second : kEmptyPool;

      topologyAssembler.BuildVolumeTopo(cell2Node, bdryCellPool, *cellPool, oldCoords, volumeTopology);
    }
  }

  topologyAssembler.AppendUnclippedSurfaces(mesh, be.surfaceCells, oldCoords, volumeTopology);

  volumeTopology.cellParent[FSMeshEnums::CellType::CT_Poly2D] =
    surfaceTopology.cellParent[FSMeshEnums::CellType::CT_Poly2D];
  volumeTopology.cellParentType[FSMeshEnums::CellType::CT_Poly2D] =
    surfaceTopology.cellParentType[FSMeshEnums::CellType::CT_Poly2D];
  volumeTopology.poly2DGlobalNumbers = surfaceTopology.poly2DGlobalNumbers;

  return volumeTopology;
}

//-----------------------------------------------------------------------------
//
//  Reconstruct
//

FSMesh FSClippingEngine::Reconstruct(FSClac& clac, FSMesh& meshOriginal, FSTopologyData& topo, bool copyAttributes)
{
  FSMeshReconstruction meshReconstruction(clac, tol_);
  FSMesh clippedMesh = meshReconstruction.Build(meshOriginal.GetMeshData()->GetUnstructCells(), topo);
  if(copyAttributes)
    meshReconstruction.CopyAttributes(meshOriginal, clippedMesh);
  return clippedMesh;
}

//-----------------------------------------------------------------------------
//
//  RunMatcherProc
//

void FSClippingEngine::RunMatcherProc(FSClac& globalClac, FSClac& localClac, FS_floatT tol)
{
  auto allGathered = FSFaceExchange::GatherReceiveAll(globalClac);

  FSFaceMatcher matcher(localClac, allGathered.meshA.faces.faces, allGathered.meshB.faces.faces, tol);
  std::vector<FSFaceMatch> matches;
  matcher.ComputeMatches(matches);

  // Each side gets its own global numbering (nodes + Poly2D cells), assigned
  // here once so that mesh procs sharing an interface node receive the same ID.
  FSFaceMatcher::AssignGlobalNodeIds(matches, tol);
  FSMatchExchange::ScatterSend(globalClac, matches, allGathered.meshA);

  matcher.ComputeInvertedMatches(matches);
  FSFaceMatcher::AssignGlobalNodeIds(matches, tol);
  FSMatchExchange::ScatterSend(globalClac, matches, allGathered.meshB);
}

_FS_END_NAMESPACE

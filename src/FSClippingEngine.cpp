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
  for(const auto& t : mesh.GetCellTypes()) {
    if(FSMeshEnums::IsUnstructVolumeCellType(t)) {
      const auto& cell2Node = mesh.GetCell2Node(t);
      const auto& cellPool = mesh.GetMeshData()->GetUnstructCells().GetCellPool(t);
      const auto& bdryCellPool = be.volumeCells.at(t);

      topologyAssembler.BuildVolumeTopo(cell2Node, bdryCellPool, *cellPool, oldCoords, volumeTopology);
    }
  }

  topologyAssembler.AppendUnclippedSurfaces(mesh, be.surfaceCells, oldCoords, volumeTopology);

  volumeTopology.cellParent[FSMeshEnums::CellType::CT_Poly2D] =
    surfaceTopology.cellParent[FSMeshEnums::CellType::CT_Poly2D];
  volumeTopology.cellParentType[FSMeshEnums::CellType::CT_Poly2D] =
    surfaceTopology.cellParentType[FSMeshEnums::CellType::CT_Poly2D];

  return volumeTopology;
}

//-----------------------------------------------------------------------------
//
//  Reconstruct
//

FSMesh FSClippingEngine::Reconstruct(FSClac& clac, FSMesh& meshOriginal, FSTopologyData& topo, bool copyAttributes)
{
  FSMeshReconstruction meshReconstruction(clac);
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
  auto subjectReceived = FSFaceExchange::Receive(globalClac, 0);
  auto clippedReceived = FSFaceExchange::Receive(globalClac, 1);

  FSFaceMatcher matcher(localClac, subjectReceived.faces, clippedReceived.faces, tol);
  std::vector<FSFaceMatch> matches;
  matcher.ComputeMatches(matches);

  FSMatchExchange::Send(globalClac, 0, matches);
  matcher.ComputeInvertedMatches(matches);
  FSMatchExchange::Send(globalClac, 1, matches);
}

_FS_END_NAMESPACE

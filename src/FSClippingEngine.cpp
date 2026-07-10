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
  auto allGathered = FSFaceExchange::GatherReceiveAll(globalClac);
  std::cout << "Nb faces receive mesh A : " << allGathered.meshA.faces.faces.size() << std::endl;
  std::cout << "Nb faces receive mesh B : " << allGathered.meshB.faces.faces.size() << std::endl;

  FSFaceMatcher matcher(localClac, allGathered.meshA.faces.faces, allGathered.meshB.faces.faces, tol);
  std::vector<FSFaceMatch> matches;
  matcher.ComputeMatches(matches);


  // for(const auto& m : matches) {
  //   std::cout << "Matches " << m.type << " own by elem1 " << m.elemOwner1 << " and elem2 :  " << m.elemOwner2
  //             << std::endl;
  // }
  //   ScatterSend routes original matches to mesh A procs and InvertSingle(match) to mesh B procs.
  //   InvertSingle swaps owners and clippedPoly3D <-> clippedPoly3D_face2, equivalent to InvertMatches.
  FSMatchExchange::ScatterSend(globalClac, matches, allGathered.meshA);
  matcher.ComputeInvertedMatches(matches);
  std::cout << "Nb matches extract final : " << matches.size() << std::endl;
  FSMatchExchange::ScatterSend(globalClac, matches, allGathered.meshB);
}

_FS_END_NAMESPACE

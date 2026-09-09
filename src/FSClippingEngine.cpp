#include "FSClipping/FSClippingEngine.h"
#include "FSClipping/FSFaceExchange.h"
#include "FSClipping/FSMeshReconstruction.h"
#include <FSMeshData.h>
#include <numeric>

_FS_BEGIN_NAMESPACE

//-----------------------------------------------------------------------------
//
//  BuildTopology
//

FSTopologyData FSClippingEngine::BuildTopology(FSMesh& mesh, const BoundaryExtraction& be,
                                               const std::vector<FSFaceMatch>& matches, Mode mode)
{
  // const FS_intT master = mesh.GetClac()->GetProcID();
  FSTopologyAssembler topologyAssembler(tol_);
  FSTopologyData surfaceTopology = topologyAssembler.BuildSurfaceTopo(matches, be.faceKeys);
  // FS_int32T nProc = static_cast<FS_int32T>(1);
  // std::vector<FS_int32T> countsProc(mesh.GetClac()->NProcs(), 0);
  // mesh.GetClac()->AllGather(&nProc, 1, countsProc.data(), 1);
  // FS_intT numProc = std::accumulate(countsProc.begin(), countsProc.end(), 0);
  // if(!master)
  //   std::cout << "Surfacique done" << numProc << std::endl;
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
      const auto it = be.volumeCells.find(t);
      const auto& bdryCellPool = (it != be.volumeCells.end()) ? it->second : std::set<FS_intT>();

      topologyAssembler.BuildVolumeTopo(cell2Node, bdryCellPool, *cellPool, oldCoords, volumeTopology);
    }
  }
  topologyAssembler.AppendUnclippedSurfaces(mesh, be.surfaceCells, oldCoords, volumeTopology);
  volumeTopology.cellParent[FSMeshEnums::CellType::CT_Poly2D] =
    surfaceTopology.cellParent[FSMeshEnums::CellType::CT_Poly2D];
  volumeTopology.cellParentType[FSMeshEnums::CellType::CT_Poly2D] =
    surfaceTopology.cellParentType[FSMeshEnums::CellType::CT_Poly2D];
  volumeTopology.poly2DGlobalNumbers = surfaceTopology.poly2DGlobalNumbers;
  volumeTopology.poly2DMarkers = surfaceTopology.poly2DMarkers;

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

bool FSClippingEngine::RunMatcherProc(FSClac& globalClac, FSClac& localClac, FS_floatT tol)
{
  auto allGathered = FSFaceExchange::GatherReceiveAll(globalClac);
  std::cout << "Nb faces mesh A :" << allGathered.meshA.faces.faces.size() << std::endl;
  std::cout << "Nb faces mesh B :" << allGathered.meshB.faces.faces.size() << std::endl;
  if(allGathered.meshA.faces.faces.empty() || allGathered.meshB.faces.faces.empty())
    return false;

  FSFaceMatcher matcher(localClac, allGathered.meshA.faces.faces, allGathered.meshB.faces.faces, tol);
  std::vector<FSFaceMatch> matches;
  matcher.ComputeMatches(matches);

  if(matches.empty()) {
    return false;
  }

  // Each side gets its own global numbering (nodes + Poly2D cells)
  FSFaceMatcher::AssignGlobalNodeIds(matches, tol);
  FSMatchExchange::ScatterSend(globalClac, matches, allGathered.meshA);

  matcher.ComputeInvertedMatches(matches);
  FSFaceMatcher::AssignGlobalNodeIds(matches, tol);
  FSMatchExchange::ScatterSend(globalClac, matches, allGathered.meshB);
  std::cout << "Nb matches clipper proc :" << matches.size() << std::endl;
  return true;
}

_FS_END_NAMESPACE

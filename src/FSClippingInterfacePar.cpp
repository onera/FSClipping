#include "FSClipping/FSClippingInterfacePar.h"
#include "FSClipping/FSFaceExchange.h"
#include "FSClipping/FSMeshReconstruction.h"
#include "FSClipping/FSTopologyAssembler.h"
#include "FSClipping/FSFaceMatcher.h"
#include <FSMeshData.h>

_FS_BEGIN_NAMESPACE

FSMesh FSClippingInterfacePar::BuildSurfaceInterface(FSMesh& mesh)
{

  if(globalClac_.GetNProcs() != 3)
    FSError.SetAndPrintAndExit("3 MPI processes are required");


  FS_intT meshId = globalClac_.GetProcID();


  if(meshId == 0 || meshId == 1) {

    FSMeshFaceExtractor extraction;
    auto boundaryExtraction = FSBoundaryFaceProvider::Extract(mesh, extraction, boundaryMarkerMesh_, tol_, false);

    FSFaceExchange::Send(globalClac_, 2, boundaryExtraction.faces);

    std::vector<FSFaceMatch> matches = FSMatchExchange::Receive(globalClac_, 2);

    FSTopologyAssembler topologyAssembler(tol_);
    FSTopologyData surfaceTopology = topologyAssembler.BuildSurfaceTopo(matches, boundaryExtraction.faceKeys);

    FSMeshReconstruction meshReconstruction(clac_);
    return meshReconstruction.Build(mesh.GetMeshData()->GetUnstructCells(), surfaceTopology);
  }

  else {
    auto subjectReceived = FSFaceExchange::Receive(globalClac_, 0);
    auto clippedReceived = FSFaceExchange::Receive(globalClac_, 1);

    // --- ComputeMatches (clac is 1-proc → BVH purely local) ---
    FSFaceMatcher matcher(clac_, subjectReceived.faces, clippedReceived.faces, tol_);
    std::vector<FSFaceMatch> matches;
    matcher.ComputeMatches(matches);

    FSMatchExchange::Send(globalClac_, 0, matches);
    matcher.InvertMatches(matches);
    FSMatchExchange::Send(globalClac_, 1, matches);
  }

  return FSMesh(&clac_);
}

FSMesh FSClippingInterfacePar::BuildVolumeInterface(FSMesh& mesh)
{

  if(globalClac_.GetNProcs() != 3)
    FSError.SetAndPrintAndExit("3 MPI processes are required");


  FS_intT meshId = globalClac_.GetProcID();


  if(meshId == 0 || meshId == 1) {

    FSMeshFaceExtractor extraction;
    auto boundaryExtraction = FSBoundaryFaceProvider::Extract(mesh, extraction, boundaryMarkerMesh_, tol_);

    FSFaceExchange::Send(globalClac_, 2, boundaryExtraction.faces);

    std::vector<FSFaceMatch> matches = FSMatchExchange::Receive(globalClac_, 2);

    FSTopologyAssembler topologyAssembler(tol_);
    FSTopologyData surfaceTopology = topologyAssembler.BuildSurfaceTopo(matches, boundaryExtraction.faceKeys);

    // Get the coordinate of the original mesh
    FSUnstructMeshData& meshData = mesh.GetMeshData()->GetUnstructCells();
    FSFloatArrayT oldCoords;
    FS_intT nodeOffset;
    FSQuantityDescArrayT coordDesc;
    meshData.GetCoordinates3D(coordDesc, oldCoords, nodeOffset);

    FSTopologyData volumeTopology;
    for(const auto& t : mesh.GetCellTypes()) {
      if(FSMeshEnums::IsUnstructVolumeCellType(t)) {
        const auto& cell2Node1 = mesh.GetCell2Node(t);
        const auto& cellPool1 = mesh.GetMeshData()->GetUnstructCells().GetCellPool(t);
        const auto& bdryCellPool1 = boundaryExtraction.volumeCells.at(t);

        topologyAssembler.BuildVolumeTopo(cell2Node1, bdryCellPool1, *cellPool1, oldCoords, volumeTopology);
      }
    }

    topologyAssembler.AppendUnclippedSurfaces(mesh, boundaryExtraction.surfaceCells, oldCoords, volumeTopology);

    volumeTopology.cellParent[FSMeshEnums::CellType::CT_Poly2D] =
      surfaceTopology.cellParent[FSMeshEnums::CellType::CT_Poly2D];
    volumeTopology.cellParentType[FSMeshEnums::CellType::CT_Poly2D] =
      surfaceTopology.cellParentType[FSMeshEnums::CellType::CT_Poly2D];
    FSMeshReconstruction meshReconstruction(clac_);
    return meshReconstruction.Build(mesh.GetMeshData()->GetUnstructCells(), volumeTopology);
  }

  else {
    auto subjectReceived = FSFaceExchange::Receive(globalClac_, 0);
    auto clippedReceived = FSFaceExchange::Receive(globalClac_, 1);

    // --- ComputeMatches (clac is 1-proc → BVH purely local) ---
    FSFaceMatcher matcher(clac_, subjectReceived.faces, clippedReceived.faces, tol_);
    std::vector<FSFaceMatch> matches;
    matcher.ComputeMatches(matches);

    FSMatchExchange::Send(globalClac_, 0, matches);
    matcher.InvertMatches(matches);
    FSMatchExchange::Send(globalClac_, 1, matches);
  }

  return FSMesh(&clac_);
}
_FS_END_NAMESPACE
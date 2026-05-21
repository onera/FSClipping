#include "FSClipping/FSClippingInterfacePar.h"
#include "FSClipping/FSFaceExchange.h"
#include "FSClipping/FSMeshReconstruction.h"
#include "FSClipping/FSTopologyAssembler.h"
#include <FSMeshData.h>

_FS_BEGIN_NAMESPACE

FSMesh FSClippingInterfacePar::BuildSurfaceInterface(FSMesh& mesh)
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
    if(meshId == 0) {
      FSMeshReconstruction meshReconstruction(clac_);
      meshReconstruction.Build(mesh.GetMeshData()->GetUnstructCells(), surfaceTopology);
    }
    return mesh;
  }

  else {
    auto subjectReceived = FSFaceExchange::Receive(globalClac_, 0);
    auto clippedReceived = FSFaceExchange::Receive(globalClac_, 1);

    // --- ComputeMatches (clac is 1-proc → BVH purely local) ---
    FSFaceMatcher matcher(clac_, subjectReceived.faces, clippedReceived.faces, tol_);
    std::vector<FSFaceMatch> matches;
    matcher.ComputeMatches(matches);

    FSMatchExchange::Send(globalClac_, 0, matches);
    FSMatchExchange::Send(globalClac_, 1, matches);
    // // - Build surface topology (faceKeys not needed for surface-only) ---
    // FSTopologyAssembler assembler(tol_);
    // const std::unordered_set<GeomFaceKey, GeomFaceKeyHash> emptyFaceKeys;
    // FSTopologyData surfaceTopo = assembler.BuildSurfaceTopo(matches, emptyFaceKeys);
  }



  return mesh;
}

_FS_END_NAMESPACE
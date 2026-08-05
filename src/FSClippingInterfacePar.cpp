#include "FSClipping/FSClippingInterfacePar.h"
#include "FSClipping/FSBoundaryFaceProvider.h"
#include "FSClipping/FSFaceExchange.h"

_FS_BEGIN_NAMESPACE

FSMesh FSClippingInterfacePar::BuildSurfaceInterface(FSMesh& mesh, FS_intT meshID)
{
  return BuildInterface(mesh, FSClippingEngine::Mode::Surface, true, meshID);
}

FSMesh FSClippingInterfacePar::BuildVolumeInterface(FSMesh& mesh, FS_intT meshID)
{
  return BuildInterface(mesh, FSClippingEngine::Mode::Volume, true, meshID);
}

FSMesh FSClippingInterfacePar::BuildInterface(FSMesh& mesh, FSClippingEngine::Mode mode, bool matchRemoteFaces,
                                              const FS_intT meshID)
{
  // proc 0: matcher role (clac_ is 1-proc → BVH purely local)
  if(meshID == 0) {
    FSClippingEngine::RunMatcherProc(globalClac_, clac_, tol_);
    return FSMesh(&clac_);
  }

  FSMeshFaceExtractor extractor;
  auto boundaryExtraction =
    FSBoundaryFaceProvider::Extract(mesh, extractor, boundaryMarkerMesh_, tol_, matchRemoteFaces);

  FSFaceExchange::GatherSend(globalClac_, 0, meshID, boundaryExtraction.faces);
  std::vector<FSFaceMatch> matches = FSMatchExchange::ScatterReceive(globalClac_, 0);
  FSClippingEngine engine(tol_);
  FSTopologyData topology = engine.BuildTopology(mesh, boundaryExtraction, matches, mode);
  return engine.Reconstruct(clac_, mesh, topology);
}

_FS_END_NAMESPACE

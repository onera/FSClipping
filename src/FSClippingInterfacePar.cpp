#include "FSClipping/FSClippingInterfacePar.h"
#include "FSClipping/FSBoundaryFaceProvider.h"
#include "FSClipping/FSFaceExchange.h"

_FS_BEGIN_NAMESPACE

FSMesh FSClippingInterfacePar::BuildSurfaceInterface(FSMesh& mesh)
{
  return BuildInterface(mesh, FSClippingEngine::Mode::Surface, false);
}

FSMesh FSClippingInterfacePar::BuildVolumeInterface(FSMesh& mesh)
{
  return BuildInterface(mesh, FSClippingEngine::Mode::Volume, true);
}

FSMesh FSClippingInterfacePar::BuildInterface(FSMesh& mesh, FSClippingEngine::Mode mode, bool matchRemoteFaces)
{
  if(globalClac_.GetNProcs() != 3)
    FSError.SetAndPrintAndExit("3 MPI processes are required");

  const FS_intT meshId = globalClac_.GetProcID();

  // proc 2: matcher role (clac_ is 1-proc → BVH purely local)
  if(meshId == 2) {
    FSClippingEngine::RunMatcherProc(globalClac_, clac_, tol_);
    return FSMesh(&clac_);
  }

  // procs 0 and 1: extract, exchange with the matcher, rebuild
  FSMeshFaceExtractor extractor;
  auto boundaryExtraction = FSBoundaryFaceProvider::Extract(mesh, extractor, boundaryMarkerMesh_, tol_, matchRemoteFaces);

  FSFaceExchange::Send(globalClac_, 2, boundaryExtraction.faces);
  std::vector<FSFaceMatch> matches = FSMatchExchange::Receive(globalClac_, 2);

  FSClippingEngine engine(tol_);
  FSTopologyData topology = engine.BuildTopology(mesh, boundaryExtraction, matches, mode);
  return engine.Reconstruct(clac_, mesh, topology);
}

_FS_END_NAMESPACE

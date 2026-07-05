#include "FSClipping/FSClippedMesh.h"
#include "FSClipping/FSClippedMeshParams.h"

#include "FSClipping/FSClippingEngine.h"
#include "FSClipping/FSFaceExchange.h"
#include "FSDataManagerData.h"
#include "FSTimer.h"


_FS_BEGIN_NAMESPACE

//-----------------------------------------------------------------------------
//
//  Create
//

FSDataManagerOp* FSClippedMesh::Create(FSClac* clac) { return new FSClippedMesh(clac); }

//-----------------------------------------------------------------------------
//
//  Constructor
//

FSClippedMesh::FSClippedMesh(FSClac* clac) : FSDataManagerOp(clac) {}


//-----------------------------------------------------------------------------
//
//  Destructor
//

FSClippedMesh::~FSClippedMesh() {}


//-----------------------------------------------------------------------------
//
//  GetClassName
//

FSString FSClippedMesh::GetClassName() const { return "FSClippedMesh"; }


//-----------------------------------------------------------------------------
//
//  DoOp
//

bool FSClippedMesh::DoOp(FSDataManagerData*& data, const FSDataManagerOpParams* params)
{
  const FSClippedMeshParams* clippedMeshParams = dynamic_cast<const FSClippedMeshParams*>(params);

  if(clippedMeshParams == nullptr) {
    FSError("FSClippedMesh: clipped mesh parameters NOT initialized");
    return false;
  }

  if(mClac->GetNProcs() != 3)
    FSError.SetAndPrintAndExit("FSClippedMesh: 3 MPI processes required");

  FSTimer timer(mClac, sTimerLevel);
  timer.Start();

  const FS_intT meshId = mClac->GetProcID();
  FSClac meshClac(MPI_COMM_SELF);

  mData = data;

  // proc 2 (matcher) has no DataManager data — only procs 0 and 1 require it
  if(meshId != 2 && mData == nullptr) {
    FSError("FSClippedMesh: data manager data NOT initialized.");
    return false;
  }

  mParams = *clippedMeshParams;

  if(!mParams.IsInitialized())
    return false;

  const bool okFlag = GenerateClippedMesh(meshId, meshClac);
  timer.Stop();
  timer.Print(0, "FSClippedMesh: created clipped mesh:");
  return okFlag;
}

//-----------------------------------------------------------------------------
//
//  GenerateClippedMesh
//

bool FSClippedMesh::GenerateClippedMesh(FS_intT meshId, FSClac& meshClac)
{
  if(meshId == 0 || meshId == 1) {

    const FS_intT marker = (meshId == 0) ? mParams.mMarker1 : mParams.mMarker2;
    const FSString& meshKey = (meshId == 0) ? mParams.mMeshKeyOriginal1 : mParams.mMeshKeyOriginal2;
    const FSString& clippedKey = (meshId == 0) ? mParams.mMeshKeyClipped1 : mParams.mMeshKeyClipped2;

    if(!mData->HasMesh(meshKey))
      return false;
    FSMesh* mesh = mData->GetMesh(meshKey, false);
    FSMesh* meshClipped = mData->GetMesh(clippedKey, mesh->GetClac(), true);

    BoundaryExtraction be;
    ExtractBoundaryFaces(*mesh, marker, be);

    FSFaceExchange::Send(*mClac, 2, be.faces);
    const std::vector<FSFaceMatch> matches = FSMatchExchange::Receive(*mClac, 2);

    FSClippingEngine engine(mParams.mTol);
    FSTopologyData meshClippedTopo = engine.BuildTopology(*mesh, be, matches, FSClippingEngine::Mode::Volume);
    *meshClipped = engine.Reconstruct(*mesh->GetClac(), *mesh, meshClippedTopo, true);
#ifdef FS_SAFETYCHECKS
    if(!meshClipped->Check()) {
      FSError("FSClippedMesh: resulting mesh of sub-elements is invalid.");
      return false;
    }
#endif

  } else { // meshId == 2
    FSClippingEngine::RunMatcherProc(*mClac, meshClac, mParams.mTol);
  }

  return true;
}

//-----------------------------------------------------------------------------
//
//  ExtractBoundaryFaces
//

bool FSClippedMesh::ExtractBoundaryFaces(FSMesh& mesh, FS_intT marker, BoundaryExtraction& be)
{
  FSMeshFaceExtractor& fex = (mClac->GetProcID() == 0) ? mFaceExtractor1 : mFaceExtractor2;
  be = FSBoundaryFaceProvider::Extract(mesh, fex, marker, mParams.mTol, false);

  if(be.faces.empty()) {
    FSLog("No faces were extracted from the mesh");
    return false;
  }
  return true;
}

_FS_END_NAMESPACE

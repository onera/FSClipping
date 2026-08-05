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

  FSTimer timer(mClac, sTimerLevel);
  timer.Start();

  mData = data;
  mParams = *clippedMeshParams;
  if(!mParams.IsInitialized())
    return false;

  // Determine role from which meshes are locally present in the data manager.
  const bool hasA = mData && mData->HasMesh(mParams.mMeshKeyOriginal1);
  const bool hasB = mData && mData->HasMesh(mParams.mMeshKeyOriginal2);
  const FS_intT meshId = hasA ? 1 : (hasB ? 2 : 0);

  mParams = *clippedMeshParams;

  if(!mParams.IsInitialized())
    return false;

  const bool okFlag = GenerateClippedMesh(meshId);
  timer.Stop();
  timer.Print(0, "FSClippedMesh: created clipped mesh:");
  return okFlag;
}

//-----------------------------------------------------------------------------
//
//  GenerateClippedMesh
//

bool FSClippedMesh::GenerateClippedMesh(FS_intT meshId)
{
  constexpr FS_intT clipperProc = 0;

  if(meshId == 1 || meshId == 2) {

    const FS_intT marker = (meshId == 1) ? mParams.mMarker1 : mParams.mMarker2;
    const FSString& meshKey = (meshId == 1) ? mParams.mMeshKeyOriginal1 : mParams.mMeshKeyOriginal2;
    const FSString& clippedKey = (meshId == 1) ? mParams.mMeshKeyClipped1 : mParams.mMeshKeyClipped2;

    FSMesh* mesh = mData->GetMesh(meshKey, false);
    FSMesh* meshClipped = mData->GetMesh(clippedKey, mesh->GetClac(), true);

    BoundaryExtraction be;
    ExtractBoundaryFaces(*mesh, marker, be);

    FSFaceExchange::GatherSend(*mClac, clipperProc, meshId, be.faces);
    const std::vector<FSFaceMatch> matches = FSMatchExchange::ScatterReceive(*mClac, clipperProc);

    FSClippingEngine engine(mParams.mTol);
    FSTopologyData meshClippedTopo = engine.BuildTopology(*mesh, be, matches, FSClippingEngine::Mode::Volume);
    *meshClipped = engine.Reconstruct(*mesh->GetClac(), *mesh, meshClippedTopo, true);
#ifdef FS_SAFETYCHECKS
    if(!meshClipped->Check()) {
      FSError("FSClippedMesh: resulting mesh of sub-elements is invalid.");
      return false;
    }
#endif

  } else { // meshId == 0
    FSClac selfClac(FSClac::sSelfComm);
    FSClippingEngine::RunMatcherProc(*mClac, selfClac, mParams.mTol);
  }

  return true;
}

//-----------------------------------------------------------------------------
//
//  ExtractBoundaryFaces
//

bool FSClippedMesh::ExtractBoundaryFaces(FSMesh& mesh, FS_intT marker, BoundaryExtraction& be)
{
  FSMeshFaceExtractor& fex = (mClac->GetProcID() == 1) ? mFaceExtractor1 : mFaceExtractor2;
  be = FSBoundaryFaceProvider::Extract(mesh, fex, marker, mParams.mTol, true);

  return true;
}

_FS_END_NAMESPACE

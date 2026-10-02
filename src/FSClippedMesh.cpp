#include "FSClipping/FSClippedMesh.h"
#include "FSClipping/FSClippedMeshParams.h"

#include "FSClipping/FSClippingEngine.h"
#include "FSClipping/FSFaceExchange.h"
#include "FSDataManagerData.h"
#include "FSTimer.h"
#include <numeric>


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

  // --- check ---
  if(clippedMeshParams == nullptr) {
    FSError("FSClippedMesh: clipped mesh parameters NOT initialized");
    return false;
  }

  // --- timer ---
  FSTimer timer(mClac, sTimerLevel);
  timer.Start();

  // --- copy data ---
  mData = data;
  // if(mData == nullptr) {
  //   FSError("FSClippedMesh: data manager data NOT initialized.");
  //   return false;
  // }

  // --- set parameters ---
  mParams = *clippedMeshParams;

  if(!mParams.IsInitialized()) {
    FSError("FSClippedMesh: clipped mesh parameters NOT initialized.");
    return false;
  }

  // Determine role from which meshes are locally present in the data manager.
  const bool hasA = mData && mData->HasMesh(mParams.mMeshKeyOriginal1);
  const bool hasB = mData && mData->HasMesh(mParams.mMeshKeyOriginal2);

  if(mClac->WorldProcID()) { // clipper proc is 0
    if(!hasA && !hasB) {
      FSError("FSClippedMesh: original mesh NOT initialized");
      return false;
    }
  }

  const FS_intT meshId = hasA ? 1 : (hasB ? 2 : 0);

  bool okFlag = GenerateClippedMesh(meshId);

  const intT nProcs = FSCLAC_NPROCS(mClac);
  if(nProcs > 1)
    mClac->AgreeOnSuccess(okFlag);

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

    // --- check input mesh
    FSMesh* mesh = mData->GetMesh(meshKey, false);
    FSMesh* meshClipped = mData->GetMesh(clippedKey, mesh->GetClac(), true);

    assert(mesh != nullptr);
    assert(meshClipped != nullptr);

    if((!mesh->IsInitialized()) || (!mesh->IsUnstructured())) {
      FSError("FSClippedMesh: original mesh NOT initialized or NOT unstructured");
      return false;
    }

    // 1 - Extract the boundaries faces with the marker
    FSMeshFaceExtractor fex;
    BoundaryExtraction be = FSBoundaryFaceProvider::Extract(*mesh, fex, marker, mParams.mTol);

    FS_int32T nFaces32 = static_cast<FS_int32T>(be.faces.size());
    std::vector<FS_int32T> countsFaces(mesh->GetClac()->NProcs(), 0);
    mesh->GetClac()->AllGather(&nFaces32, 1, countsFaces.data(), 1);
    if(!(std::accumulate(countsFaces.begin(), countsFaces.end(), 0) > 0)) {
      FSError("FSClippedMesh: no faces were extract, check your boundaries marker");
      return false;
    }

    // 2 - Send the faces extrated and compute the matches
    FSFaceExchange::GatherSend(*mClac, clipperProc, meshId, be.faces);
    const std::vector<FSFaceMatch> matches = FSMatchExchange::ScatterReceive(*mClac, clipperProc);

    // 4 - Reconstruct the topologie of the entire mesh in parrallel
    FSClippingEngine engine(mParams.mTol);
    FSTopologyData meshClippedTopo = engine.BuildTopology(*mesh, be, matches);

    // 5 - Construct the new clipped mesh
    *meshClipped = engine.Reconstruct(*mesh->GetClac(), *mesh, meshClippedTopo, true);

#ifdef FS_SAFETYCHECKS
    if(!meshClipped->Check()) {
      FSError("FSClippedMesh: resulting of clipped mesh is invalid.");
      return false;
    }
#endif

  } else {
    // 3 - Clipper proc (meshID = 0) receive the extrated faces and compute the matches
    FSClac selfClac(FSClac::sSelfComm);
    bool ok = FSClippingEngine::RunMatcherProc(*mClac, selfClac, mParams.mTol);
    if(!ok) {
      FSError("FSClippedMesh: They were an error with the clipper proc, any matches was computed, check your marker or "
              "the tolerance.");
      return false;
    }
  }

  return true;
}

_FS_END_NAMESPACE

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

    assert(mesh != nullptr);
    assert(meshClipped != nullptr);

    if((!mesh->IsInitialized()) || (!mesh->IsUnstructured())) {
      FSString msg = "FSClippedMesh: original mesh " + meshKey + "NOT initialized";
      FSError(msg);
      return false;
    }

    // 1 - Extract the boundaries faces with the marker
    FSMeshFaceExtractor fex;
    BoundaryExtraction be = FSBoundaryFaceProvider::Extract(*mesh, fex, marker, mParams.mTol);
    // condition sur les bdry sinon return false, operation atomique pour recuperer le nb total de face extraite ?

    // 2 - Send the faces extrated and compute the matches
    FSFaceExchange::GatherSend(*mClac, clipperProc, meshId, be.faces);
    const std::vector<FSFaceMatch> matches = FSMatchExchange::ScatterReceive(*mClac, clipperProc);
    // condition sur les matches -> Si tous les matches sont nul, alors on a un pbm

    // 4 - Reconstruct the topologie of the entire mesh in parrallel
    FSClippingEngine engine(mParams.mTol);
    FSTopologyData meshClippedTopo = engine.BuildTopology(*mesh, be, matches, FSClippingEngine::Mode::Volume);

    // 5 - Construct the new clipped mesh
    *meshClipped = engine.Reconstruct(*mesh->GetClac(), *mesh, meshClippedTopo, true);

#ifdef FS_SAFETYCHECKS
    if(!meshClipped->Check()) {
      FSError("FSClippedMesh: resulting mesh of sub-elements is invalid.");
      return false;
    }
#endif

  } else {
    // 3 - Clipper proc (meshID = 0) receive the extrated faces and compute the matches
    FSClac selfClac(FSClac::sSelfComm);
    FSClippingEngine::RunMatcherProc(*mClac, selfClac, mParams.mTol);
  }

  return true;
}

_FS_END_NAMESPACE

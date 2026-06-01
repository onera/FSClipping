#include "FSClipping/FSClippedMesh.h"
#include "FSClipping/FSClippedMeshParams.h"

#include "FSClipping/FSFaceExchange.h"
#include "FSClipping/FSFaceMatcher.h"
#include "FSClipping/FSMeshReconstruction.h"
#include "FSDataManagerData.h"
#include "FSTimer.h"
#include <FSMeshData.h>


_FS_BEGIN_NAMESPACE

//-----------------------------------------------------------------------------
//
//  Create
//

FSDataManagerOp* FSClippedMesh::Create(FSClac* clac)
{
  return new FSClippedMesh(clac);
}

//-----------------------------------------------------------------------------
//
//  Constructor
//

FSClippedMesh::FSClippedMesh(FSClac* clac)
  : FSDataManagerOp(clac)
{
}


//-----------------------------------------------------------------------------
//
//  Destructor
//

FSClippedMesh::~FSClippedMesh()
{
}


//-----------------------------------------------------------------------------
//
//  GetClassName
//

FSString FSClippedMesh::GetClassName() const
{
  return "FSClippedMesh";
}


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

    FSTopologyData meshClippedTopo;
    if(!GenerateMeshClippedTopo(*mesh, be, matches, meshClippedTopo))
      return false;
    if(!GenerateMesh(*mesh, meshClippedTopo, *meshClipped))
      return false;
#ifdef FS_SAFETYCHECKS
    if(!meshClipped->Check()) {
      FSError("FSClippedMesh: resulting mesh of sub-elements is invalid.");
      return false;
    }
#endif

  } else { // meshId == 2

    auto subjectReceived = FSFaceExchange::Receive(*mClac, 0);
    auto clippedReceived = FSFaceExchange::Receive(*mClac, 1);

    FSFaceMatcher matcher(meshClac, subjectReceived.faces, clippedReceived.faces, mParams.mTol);
    std::vector<FSFaceMatch> matches;
    matcher.ComputeMatches(matches);
    FSMatchExchange::Send(*mClac, 0, matches);
    matcher.InvertMatches(matches);
    FSMatchExchange::Send(*mClac, 1, matches);
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

//-----------------------------------------------------------------------------
//
//  GenerateMeshClippedTopo
//

bool FSClippedMesh::GenerateMeshClippedTopo(FSMesh& mesh, const BoundaryExtraction& be1,
                                            const std::vector<FSFaceMatch>& matches,
                                            FSTopologyData& meshClippedTopo)
{
  FSTopologyAssembler topologyAssembler(mParams.mTol);
  FSTopologyData surfaceClippedTopo = topologyAssembler.BuildSurfaceTopo(matches, be1.faceKeys);

  FSUnstructMeshData& meshData = mesh.GetMeshData()->GetUnstructCells();
  FSFloatArrayT oldCoords;
  FS_intT nodeOffset;
  FSQuantityDescArrayT coordDesc;
  meshData.GetCoordinates3D(coordDesc, oldCoords, nodeOffset);

  for(const auto& t : mesh.GetCellTypes()) {
    if(FSMeshEnums::IsUnstructVolumeCellType(t)) {
      const auto& cell2Node = mesh.GetCell2Node(t);
      const auto& cellPool = meshData.GetCellPool(t);
      const auto& bdryPool = be1.volumeCells.at(t);
      topologyAssembler.BuildVolumeTopo(cell2Node, bdryPool, *cellPool, oldCoords, meshClippedTopo);
    }
  }

  topologyAssembler.AppendUnclippedSurfaces(mesh, be1.surfaceCells, oldCoords, meshClippedTopo);

  meshClippedTopo.cellParent[FSMeshEnums::CellType::CT_Poly2D] = surfaceClippedTopo.cellParent[FSMeshEnums::CellType::CT_Poly2D];
  meshClippedTopo.cellParentType[FSMeshEnums::CellType::CT_Poly2D] = surfaceClippedTopo.cellParentType[FSMeshEnums::CellType::CT_Poly2D];

  return true;
}

bool FSClippedMesh::GenerateMesh(FSMesh& meshOriginal1, FSTopologyData& meshClippedTopo, FSMesh& clippedMesh)
{
  FSClac* originalClac = meshOriginal1.GetClac();
  FSUnstructMeshData& meshDataOrig1 = meshOriginal1.GetMeshData()->GetUnstructCells();
  FSMeshReconstruction meshReconstruction(*originalClac);
  clippedMesh = meshReconstruction.Build(meshDataOrig1, meshClippedTopo);
  // clippedMesh.PrintInfo();
  return true;
}
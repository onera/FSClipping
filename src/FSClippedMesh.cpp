#include "FSClipping/FSClippedMesh.h"
#include "FSClipping/FSClippedMeshParams.h"

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
  bool okFlag = true;
  const FSClippedMeshParams* clippedMeshParams = dynamic_cast<const FSClippedMeshParams*>(params);

  if(clippedMeshParams == nullptr) {
    FSError("FSClippedMesh: clipped mesh parameters NOT initialized");
    return false;
  }

  // --- timer ---
  FSTimer timer(mClac, sTimerLevel);
  timer.Start();

  mData = data; // copy pointer

  if(mData == nullptr) {
    FSError("FSLinearSubElementsMesh: data manager data NOT initialized.");
    return false;
  }

  /// --- set parameters ---
  mParams = *clippedMeshParams;

  if(!mParams.IsInitialized())
    return false;

  if(!mData->HasMesh(mParams.mMeshKeyOriginal1) || !mData->HasMesh(mParams.mMeshKeyOriginal2))
    return false;

  // --- check input mesh ---
  FSMesh* meshOriginal1 = mData->GetMesh(mParams.mMeshKeyOriginal1, false);
  FSMesh* meshOriginal2 = mData->GetMesh(mParams.mMeshKeyOriginal2, false);
  FSMesh* meshClipped = mData->GetMesh(mParams.mMeshKeyClipped, false);

  if((!meshOriginal1->IsInitialized()) || (!meshOriginal1->IsUnstructured()) ||
     (!meshOriginal2->IsInitialized()) || (!meshOriginal2->IsUnstructured())) {
    FSError("FSClippedMesh: original mesh NOT initialized.");
    return false;
  }

  okFlag = GenerateClippedMesh(*meshOriginal1, *meshOriginal2, *meshClipped);

#ifdef FS_SAFETYCHECKS
  if(!meshClipped->Check()) {
    FSError("FSLinearSubElementsMesh: resulting mesh of sub-elements is invalid.");
    return false;
  }
#endif
  // meshClipped->Check();
  // meshClipped->PrintInfo();
  timer.Stop();
  timer.Print(0, "FSClippedMesh: created clipped mesh:");

  return okFlag;
}

//-----------------------------------------------------------------------------
//
//  GenerateSubElementsMesh
//

bool FSClippedMesh::GenerateClippedMesh(FSMesh& meshOriginal1, FSMesh& meshOriginal2, FSMesh& meshClipped)
{

  bool okFlag = true;

  // --- check input mesh ---
  if((!meshOriginal1.IsInitialized()) || (!meshOriginal1.IsUnstructured()) ||
     (!meshOriginal2.IsInitialized()) || (!meshOriginal2.IsUnstructured()))
    okFlag = false;

  BoundaryExtraction boundaryExtraction1, boundaryExtraction2;
  FSTopologyData meshClippedTopo;

  FSUnstructMeshData& meshDataOrig1 = meshOriginal1.GetMeshData()->GetUnstructCells();
  FSUnstructMeshData& meshDataOrig2 = meshOriginal2.GetMeshData()->GetUnstructCells();

  const bool originalHasLocalNumbering1 = meshDataOrig1.HasLocalNumbering();
  if(originalHasLocalNumbering1)
    meshDataOrig1.CreateGlobalNumbering();

  const bool originalHasLocalNumbering2 = meshDataOrig2.HasLocalNumbering();
  if(originalHasLocalNumbering2)
    meshDataOrig2.CreateGlobalNumbering();

  if(okFlag) {

    okFlag = ExtractBoundaryFaces(meshOriginal1, meshOriginal2, boundaryExtraction1, boundaryExtraction2);

    if(okFlag)
      okFlag = GenerateMeshClippedTopo(meshOriginal1, meshOriginal2, boundaryExtraction1, boundaryExtraction2, meshClippedTopo);

    if(okFlag)
      okFlag = GenerateMesh(meshOriginal1, meshClippedTopo, meshClipped);

    if(!okFlag) {
      FSLog("Could not create mesh of sub-elements");
      return false;
    }
  }
  return okFlag;
}

bool FSClippedMesh::ExtractBoundaryFaces(FSMesh& meshOriginal1, FSMesh& meshOriginal2, BoundaryExtraction& boundaryExtraction1, BoundaryExtraction& boundaryExtraction2)
{
  boundaryExtraction1 = FSBoundaryFaceProvider::Extract(meshOriginal1, mFaceExtractor1, mParams.mMarker1, mParams.mTol);
  boundaryExtraction2 = FSBoundaryFaceProvider::Extract(meshOriginal2, mFaceExtractor2, mParams.mMarker2, mParams.mTol);

  if(boundaryExtraction1.faces.empty() || boundaryExtraction2.faces.empty()) {
    FSLog("No faces were extract from the mesh");
    return false;
  }
  return true;
}

bool FSClippedMesh::GenerateMeshClippedTopo(FSMesh& meshOriginal1, FSMesh& meshOriginal2, BoundaryExtraction& boundaryExtraction1, BoundaryExtraction& boundaryExtraction2, FSTopologyData& meshClippedTopo)
{
  auto subjectFaces = boundaryExtraction1.faces;
  auto clippedFaces = boundaryExtraction2.faces;
  auto clippedClac = meshOriginal2.GetClac();

  // Compute the matches and the corresponding intersections
  FSFaceMatcher matcher(*clippedClac, subjectFaces, clippedFaces, mParams.mTol);
  std::vector<FSFaceMatch> matches;
  matcher.ComputeMatches(matches);

  // Compute the new surface topology of the border
  FSTopologyAssembler topologyAssembler(mParams.mTol);
  FSTopologyData surfaceClippedTopo = topologyAssembler.BuildSurfaceTopo(matches, boundaryExtraction1.faceKeys);

  // Get the coordinate of the original mesh
  FSUnstructMeshData& meshDataOrig1 = meshOriginal1.GetMeshData()->GetUnstructCells();
  FSFloatArrayT oldCoords;
  FS_intT nodeOffset;
  FSQuantityDescArrayT coordDesc;
  meshDataOrig1.GetCoordinates3D(coordDesc, oldCoords, nodeOffset);

  const auto& cellType = meshOriginal1.GetCellTypes();
  for(const auto& t : cellType) {
    if(FSMeshEnums::IsUnstructVolumeCellType(t)) {
      const auto& cell2Node1 = meshOriginal1.GetCell2Node(t);
      const auto& cellPool1 = meshOriginal1.GetMeshData()->GetUnstructCells().GetCellPool(t);
      const auto& bdryCellPool1 = boundaryExtraction1.volumeCells.at(t);

      topologyAssembler.BuildVolumeTopo(cell2Node1, bdryCellPool1, *cellPool1, oldCoords, meshClippedTopo);
      // volumeTopology = topolyBuilder.BuildVolumeTopo(cell2Node, surfaceTopology); in the future
    }
  }

  topologyAssembler.AppendUnclippedSurfaces(meshOriginal1, boundaryExtraction1.surfaceCells, oldCoords, meshClippedTopo);

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

  return true;
}
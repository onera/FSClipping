/*
Attributes of the clipped Poly2D cells, run at exactly 9 MPI procs. It is the smallest count at which RCB splits a
volume/surface pair on these meshes, which is what this test is about.

A Poly2D is born from a match, and its parent is the surface cell the matched face closes. When RCB splits a
volume/surface pair, that parent is owned by another proc: its index means nothing in the local pool, so the
CADGroupID cannot be read back locally and has to travel with the match. Without that, the copy either reads an
unrelated local cell or runs past the end of the pool.
*/

#include "FSClipping/FSClippedMesh.h"
#include "FSClipping/FSClippedMeshParams.h"
#include "TestUtilsParallel.hpp"
#include "gtest/gtest.h"
#include <FSDataManagerData.h>
#include <FSMeshImportParamsHDF5.h>

_FS_BEGIN_NAMESPACE

namespace {

/// Counts the Poly2D cells of `mesh` whose CADGroupID differs from `expectedMarker`.
FS_int32T CountWrongPoly2DMarkers(const FSMesh& mesh, FS_intT expectedMarker)
{
  static const FSString attrCADGroupIDName = FSEnums::AttributeTypeToString(FSEnums::AT_CADGroupID);
  const FSUnstructMeshData& meshData = mesh.GetMeshData()->GetUnstructCells();

  const FSCellPool* cellPool = meshData.GetCellPool(FSMeshEnums::CT_Poly2D);
  if(!cellPool || meshData.GetCellAttributes(FSMeshEnums::CT_Poly2D).Size() == 0)
    return 0;

  const FSIntArrayT& markers = meshData.GetCellAttribute(attrCADGroupIDName, FSMeshEnums::CT_Poly2D);

  FS_int32T wrong = 0;
  for(FS_intT i = 0; i < markers.Size(); ++i)
    if(markers[i] != expectedMarker)
      ++wrong;
  return wrong;
}

/// Counts the Quad4 cells of `mesh` still carrying `marker`: each one should have been replaced by Poly2D.
FS_int32T CountQuad4WithMarker(const FSMesh& mesh, FS_intT marker)
{
  static const FSString attrCADGroupIDName = FSEnums::AttributeTypeToString(FSEnums::AT_CADGroupID);
  const FSUnstructMeshData& meshData = mesh.GetMeshData()->GetUnstructCells();
  if(!meshData.GetCellPool(FSMeshEnums::CT_Quad4) || !meshData.HasCellAttribute(attrCADGroupIDName, FSMeshEnums::CT_Quad4))
    return 0;

  const FSIntArrayT& markers = meshData.GetCellAttribute(attrCADGroupIDName, FSMeshEnums::CT_Quad4);
  FS_int32T n = 0;
  for(FS_intT i = 0; i < markers.Size(); ++i)
    if(markers[i] == marker)
      ++n;
  return n;
}

/// Same as polyMeshExtractFaces, but returns the result instead of asserting, so it can be reduced over the procs.
bool ExtractFaces(FSUnstructMeshData& unstructMeshData)
{
  unstructMeshData.CreateLocalNumbering(); // face nodal data needs it
  FSMeshFaceExtractor fex;
  FSIntArrayT faces;
  const bool success = fex.PrepareFaceConnectivity(unstructMeshData, false, false) &&
                       fex.PrepareFaceNodeCoordinates(FSQuantityDescArrayT()) &&
                       (faces.Resize(fex.GetNLocalFaces(), 6), fex.ExtractFaces(unstructMeshData, faces, false));
  if(!success)
    FSError.Print();
  return success;
}

} // namespace

TEST(FSClippingTestPoly2DAttributesPar, ClippedFacesKeepTheirCADGroupID)
{
  FSClac globalClac(MPI_COMM_WORLD);
  SKIP_UNLESS_EXACT_PROCS(globalClac, 9);

  const FS_intT procId = globalClac.GetProcID();
  const FS_intT meshID = SymMeshID(procId);

  FSClac clac;
  globalClac.DivideIntoGroups(meshID, clac);

  FSDataManagerData* data = new FSDataManagerData(&globalClac);

  FSClippedMeshParams clippedParams;
  clippedParams.mMeshKeyOriginal1 = "original1";
  clippedParams.mMeshKeyOriginal2 = "original2";
  clippedParams.mMeshKeyClipped1 = "clippedMesh1";
  clippedParams.mMeshKeyClipped2 = "clippedMesh2";
  clippedParams.mMarker1 = 6;
  clippedParams.mMarker2 = 5;
  clippedParams.mTol = 1e-8;

  if(meshID == 1) {
    FSMesh* meshPtr = data->GetMesh("original1", &clac);
    LoadMeshWithClac(meshPtr, MeshPath("input/cube_hexa_coarse_par.grid"));
    polyMeshRepartition(&clac, meshPtr->GetMeshData());
    meshPtr->GetMeshData()->GetUnstructCells().CreateLocalNumbering();
  }

  if(meshID == 2) {
    FSMesh* meshPtr = data->GetMesh("original2", &clac);
    LoadMeshWithClac(meshPtr, MeshPath("input/cube_hexa_fine_par.grid"));
    polyMeshRepartition(&clac, meshPtr->GetMeshData());
    meshPtr->GetMeshData()->GetUnstructCells().CreateLocalNumbering();
  }

  FSClippedMesh clippedMesh(&globalClac);
  ASSERT_TRUE(clippedMesh.DoOp(data, &clippedParams));

  if(meshID == 0)
    return; // clipper proc holds no clipped mesh

  const FSString clippedKey = (meshID == 1) ? "clippedMesh1" : "clippedMesh2";
  const FS_intT expectedMarker = (meshID == 1) ? clippedParams.mMarker1 : clippedParams.mMarker2;

  FS_int32T localWrong = CountWrongPoly2DMarkers(*data->GetMesh(clippedKey, &clac), expectedMarker);
  FS_int32T totalWrong = 0;
  clac.Sum(&localWrong, &totalWrong, 1);

  // Every clipped face lies on the clipped boundary, so all Poly2D must carry that marker
  EXPECT_EQ(totalWrong, 0);
}

/*
R5: when RCB separates a marked Quad4 from its hexa, the Quad4 must still be dropped. It used to be kept next to the
Poly2D created on the hexa side, which made the face extraction find two neighbours for the same face.
*/
TEST(FSClippingTestPoly2DAttributesPar, SeparatedSurfaceCellIsReplaced)
{
  FSClac globalClac(MPI_COMM_WORLD);
  SKIP_UNLESS_EXACT_PROCS(globalClac, 9);

  const FS_intT meshID = SymMeshID(globalClac.GetProcID());
  FSClac clac;
  globalClac.DivideIntoGroups(meshID, clac);

  FSDataManagerData* data = new FSDataManagerData(&globalClac);

  FSClippedMeshParams clippedParams;
  clippedParams.mMeshKeyOriginal1 = "original1";
  clippedParams.mMeshKeyOriginal2 = "original2";
  clippedParams.mMeshKeyClipped1 = "clippedMesh1";
  clippedParams.mMeshKeyClipped2 = "clippedMesh2";
  clippedParams.mMarker1 = 1;
  clippedParams.mMarker2 = 2;
  clippedParams.mTol = 1e-6;

  if(meshID == 1) {
    FSMesh* meshPtr = data->GetMesh("original1", &clac);
    FSMeshImportParamsHDF5 params;
    params.mMeshFilename = MeshPath("input/cube_hexa_translate.h5");
    ASSERT_TRUE(meshPtr->ImportMesh(&params));
    polyMeshRepartition(&clac, meshPtr->GetMeshData());
    meshPtr->GetMeshData()->GetUnstructCells().CreateLocalNumbering();
  }

  if(meshID == 2) {
    FSMesh* meshPtr = data->GetMesh("original2", &clac);
    LoadMeshWithClac(meshPtr, MeshPath("input/cube_hexa_fine.grid"));
    polyMeshRepartition(&clac, meshPtr->GetMeshData());
    meshPtr->GetMeshData()->GetUnstructCells().CreateLocalNumbering();
  }

  FSClippedMesh clippedMesh(&globalClac);
  ASSERT_TRUE(clippedMesh.DoOp(data, &clippedParams));

  // Reduced on every proc, clipper included, since only rank 0 reports
  FS_int32T local[2] = {0, 0}, total[2] = {0, 0}; // marked Quad4 left, failed face extractions
  if(meshID != 0) {
    FSMesh* clipped = data->GetMesh((meshID == 1) ? "clippedMesh1" : "clippedMesh2", &clac);
    local[0] = CountQuad4WithMarker(*clipped, (meshID == 1) ? clippedParams.mMarker1 : clippedParams.mMarker2);
    local[1] = ExtractFaces(clipped->GetMeshData()->GetUnstructCells()) ? 0 : 1;
  }
  globalClac.Sum(local, total, 2);
  EXPECT_EQ(total[0], 0) << "marked Quad4 kept next to its Poly2D";
  EXPECT_EQ(total[1], 0) << "face extraction failed";
}

_FS_END_NAMESPACE
